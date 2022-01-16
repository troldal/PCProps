/*

8888888b.   .d8888b.  8888888b.
888   Y88b d88P  Y88b 888   Y88b
888    888 888    888 888    888
888   d88P 888        888   d88P 888d888 .d88b.  88888b.  .d8888b
8888888P"  888        8888888P"  888P"  d88""88b 888 "88b 88K
888        888    888 888        888    888  888 888  888 "Y8888b.
888        Y88b  d88P 888        888    Y88..88P 888 d88P      X88
888         "Y8888P"  888        888     "Y88P"  88888P"   88888P'
                                                 888
                                                 888
                                                 888

Copyright (c) 2020 Kenneth Troldal Balslev

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#ifndef PCPROPS_THOMSON_HPP
#define PCPROPS_THOMSON_HPP

#include <cmath>
#include <functional>

namespace PCProps::LiquidVolume
{
    /**
     * @brief The CLVThomson class encapsulates the Thomson method for estimating compressed liquid volume.
     * @details The Thomson method looks as follows:
     * \f[ V = V_{0} \cdot \left[ 1 - C \cdot ln \ \frac{B + P}{B + P_{0}} \right] \f]
     * \f[ \frac{B}{P_{c}} = -1 +
     * a \cdot \left( 1 - \frac{T}{T_{c}} \right)^{1/3} +
     * b \cdot \left( 1 - \frac{T}{T_{c}} \right)^{2/3} +
     * c \cdot \left( 1 - \frac{T}{T_{c}} \right) +
     * d \cdot \left( 1 - \frac{T}{T_{c}} \right)^{4/3} \f]
     * \f[ C = 0.0861488 + 0.0344483 \cdot \omega \f]
     * \f[ where \f]
     * \f[ a = -9.070217 \f]
     * \f[ b = 62.45326 \f]
     * \f[ c = -135.1102 \f]
     * \f[ d = e^{4.79594 + 0.250047 \cdot \omega + 1.14188 \cdot \omega^2} \f]
     *
     * In the equations above, \f$ V_{0}\f$ and \f$ P_{0} \f$ refers to a known reference volume and pressure. Often
     * the saturation volume (\f$ V_{s} \f$) and pressure (\f$ P_{s} \f$) is used.
     * @warning This method may yield NaN as the result when T is close to Tc.
     */
    class Thomson final
    {
        double m_criticalTemperature {}; /**< The critical temperature [K]. */
        double m_criticalPressure {};    /**< The critical pressure [Pa]. */
        double m_acentricFactor {};      /**< The acentric factor [-]. */

    public:

        /**
         * @brief Constructor, taking critical properties and function objects for
         * saturated liquid volume and vapor pressure as arguments.
         * @param criticalTemperature The critical temperature [K]
         * @param criticalPressure The critical pressure [Pa]
         * @param acentricFactor The acentric factor [-]
         */
        Thomson(
            double                               criticalTemperature,
            double                               criticalPressure,
            double                               acentricFactor)
            : m_criticalTemperature(criticalTemperature),
              m_criticalPressure(criticalPressure),
              m_acentricFactor(acentricFactor)
        {}

        /**
         *
         * @tparam PC
         * @param pureComponent
         */
        template<typename PC>
        explicit Thomson(const PC& pureComponent)
            : Thomson(pureComponent.property("CriticalTemperature"),
                      pureComponent.property("CriticalPressure"),
                      pureComponent.property("AcentricFactor"))
        {}

        /**
         * @brief Copy constructor
         */
        Thomson(const Thomson& other) = default;

        /**
         * @brief Move constructor
         */
        Thomson(Thomson&& other) noexcept = default;

        /**
         * @brief Destructor
         */
        ~Thomson() = default;

        /**
         * @brief Copy assignment operator
         */
        Thomson& operator=(const Thomson& other) = default;

        /**
         * @brief Move assignment operator
         */
        Thomson& operator=(Thomson&& other) noexcept = default;

        /**
         *
         * @param params
         * @return
         */
        double operator()(std::vector<double> params) const {
            switch (params.size()) {
                case 4:
                    return operator()(params[0], params[1], params[2], params[3]);

                default:
                    return 0.0;
            }
        }

        /**
         *
         * @param temperature
         * @param pressure
         * @param satPressure
         * @param satVolume
         * @return
         */
        double operator()(double temperature, double pressure, double satPressure, double satVolume) const {

            using std::exp;
            using std::log;
            using std::pow;
            using std::min;

            double tr = min(temperature / m_criticalTemperature, 1.0);
            double C  = 0.0861488 + 0.0344483 * m_acentricFactor;
            double B  = m_criticalPressure * (-1 - 9.070217 * pow(1 - tr, 1.0 / 3.0) + 62.45326 * pow(1 - tr, 2.0 / 3.0) - 135.1102 * (1 - tr) +
                exp(4.79594 + 0.250047 * m_acentricFactor + 1.14188 * m_acentricFactor * m_acentricFactor) * pow(1 - tr, 4.0 / 3.0));

            double result = satVolume * (1 - C * log((B + pressure) / (B + satPressure)));
            return (std::isnan(result) || std::isinf(result)) ? satVolume : result;
        }
    };

}    // namespace PCProps::LiquidVolume
#endif    // PCPROPS_THOMSON_HPP
