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

#ifndef PCPROPS_AALTO_HPP
#define PCPROPS_AALTO_HPP

#include <cmath>
#include <functional>

namespace PCProps::LiquidVolume
{
    /**
     * @brief The CLVAalto class encapsulates the Aalto method for estimating compressed liquid volume.
     * @details The Aalto method looks as follows:
     * \f[ V = V_{s} \cdot
     * \frac{A \cdot P_{c} + C^{(D-T_{r})^B} \cdot (P - P_{s})}
     * {A \cdot P_{c} + C \cdot (P - P_{s})} \f]
     * \f[ A = a_{0} + a_{1} \cdot T_{r} + a_{2} \cdot T_{r}^3 + a_{3} \cdot T_{r}^6 + \frac{a_{4}}{T_{r}} \f]
     * \f[ B = b_{0} + \omega_{SRK} \cdot b_{1} \f]
     * \f[ where \f]
     * \f[ a_{0} = -170.335 \f]
     * \f[ a_{1} = -28.578 \f]
     * \f[ a_{2} = 124.809 \f]
     * \f[ a_{3} = -55.5393 \f]
     * \f[ a_{4} = 130.01 \f]
     * \f[ b_{0} = 0.164813 \f]
     * \f[ b_{1} = -0.0914427 \f]
     * \f[ C = e^1 \f]
     * \f[ D = 1.00588 \f]
     */
    class Aalto final
    {
        std::function<double(double)> m_saturatedVolumeFunction {}; /** Function object for calculation of saturated liquid volume. */
        std::function<double(double)> m_vaporPressureFunction {};   /** Function object for calculation of vapor pressure. */

        double m_criticalTemperature {}; /** The critical temperature [K]. */
        double m_criticalPressure {};    /** The critical pressure [Pa]. */

        double D { 1.00588 };
        double C { 2.718281828 };
        double B {};

    public:

        /**
         * @brief
         */
        Aalto() = default;

        /**
         * @brief Constructor, taking critical properties and function objects for
         * saturated liquid volume and vapor pressure as arguments.
         * @param criticalTemperature The critical temperature [K]
         * @param criticalPressure The critical pressure [Pa]
         * @param acentricFactor The acentric factor [-]
         * @param satVolumeFunction Function object for calculating saturated liquid volume [m3/mol] as function of temperature [K]
         * @param vaporPressureFunction Function object for calculating vapor pressure [Pa] as function of temperature [K]
         */
        Aalto(
            double                               criticalTemperature,
            double                               criticalPressure,
            double                               acentricFactor,
            const std::function<double(double)>& satVolumeFunction,
            const std::function<double(double)>& vaporPressureFunction)
            : m_saturatedVolumeFunction(satVolumeFunction),
              m_vaporPressureFunction(vaporPressureFunction),
              m_criticalTemperature(criticalTemperature),
              m_criticalPressure(criticalPressure),
              B { 0.164813 - 0.0914427 * acentricFactor }
        {}

        /**
         * @brief Copy constructor
         */
        Aalto(const Aalto& other) = default;

        /**
         * @brief Move constructor
         */
        Aalto(Aalto&& other) noexcept = default;

        /**
         * @brief Destructor
         */
        ~Aalto() = default;

        /**
         * @brief Copy assignment operator
         */
        Aalto& operator=(const Aalto& other) = default;

        /**
         * @brief Move assignment operator
         */
        Aalto& operator=(Aalto&& other) noexcept = default;

        /**
         * @brief Function call operator, taking temperature [K] and pressure [Pa] as arguments
         * and returns compressed liquid volume [m3/mol]
         * @param temperature The temperature [K]
         * @param pressure The pressure [Pa]
         * @return The compressed liquid volume [m3/mol]
         */
        double operator()(double temperature, double pressure) const
        {
            using std::pow;
            double tr = temperature / m_criticalTemperature;

            double A = -170.335 - 28.578 * tr + 124.809 * pow(tr, 3) - 55.5393 * pow(tr, 6) + 130.01 / tr;

            return m_saturatedVolumeFunction(temperature) * ((A * m_criticalPressure + pow(C, pow(D - tr, B)) * (pressure - m_vaporPressureFunction(temperature))) /
                (A * m_criticalPressure + C * (pressure - m_vaporPressureFunction(temperature))));
        }
    };
}    // namespace PCProps::LiquidVolume

#endif    // PCPROPS_AALTO_HPP
