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
        std::function<double(double)> m_saturatedVolumeFunction {}; /** Function object for calculation of saturated liquid volume. */
        std::function<double(double)> m_vaporPressureFunction {};   /** Function object for calculation of vapor pressure. */

        double m_criticalTemperature {}; /** The critical temperature [K]. */
        double m_criticalPressure {};    /** The critical pressure [Pa]. */
        double m_acentricFactor {};      /** The acentric factor [-]. */

        /**
         * @brief Constructor, taking critical properties and function objects for
         * saturated liquid volume and vapor pressure as arguments.
         * @param criticalTemperature The critical temperature [K]
         * @param criticalPressure The critical pressure [Pa]
         * @param acentricFactor The acentric factor [-]
         * @param satVolumeFunction Function object for calculating saturated liquid volume [m3/mol] as function of temperature [K]
         * @param vaporPressureFunction Function object for calculating vapor pressure [Pa] as function of temperature [K]
         */
        Thomson(
            double                               criticalTemperature,
            double                               criticalPressure,
            double                               acentricFactor,
            const std::function<double(double)>& satVolumeFunction,
            const std::function<double(double)>& vaporPressureFunction);

    public:
        /**
         * @brief Copy constructor
         */
        Thomson(const Thomson& other);

        /**
         * @brief Move constructor
         */
        Thomson(Thomson&& other) noexcept;

        /**
         * @brief Destructor
         */
        ~Thomson();

        /**
         * @brief Copy assignment operator
         */
        Thomson& operator=(const Thomson& other);

        /**
         * @brief Move assignment operator
         */
        Thomson& operator=(Thomson&& other) noexcept;

        /**
         * @brief Function call operator, taking temperature [K] and pressure [Pa] as arguments
         * and returns compressed liquid volume [m3/mol]
         * @param temperature The temperature [K]
         * @param pressure The pressure [Pa]
         * @return The compressed liquid volume [m3/mol]
         * @warning This may yield NaN as the result when T is close to Tc.
         */
        double operator()(double temperature, double pressure);

        /**
         * @brief Static factory function for creating an CLVThomson object.
         * @param criticalTemperature The critical temperature [K]
         * @param criticalPressure The critical pressure [Pa]
         * @param acentricFactor The acentric factor [-]
         * @param satVolumeFunction Function object for calculating saturated liquid volume [m3/mol] as function of temperature [K]
         * @param vaporPressureFunction Function object for calculating vapor pressure [Pa] as function of temperature [K]
         * @return An CLVThomson object
         */
        static Thomson create(
            double                               criticalTemperature,
            double                               criticalPressure,
            double                               acentricFactor,
            const std::function<double(double)>& satVolumeFunction,
            const std::function<double(double)>& vaporPressureFunction);
    };

}    // namespace PCProps::LiquidVolume
#endif    // PCPROPS_THOMSON_HPP
