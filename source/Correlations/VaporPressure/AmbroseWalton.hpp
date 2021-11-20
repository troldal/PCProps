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

#ifndef PCPROPS_AMBROSEWALTON_HPP
#define PCPROPS_AMBROSEWALTON_HPP

#include <cmath>

namespace PCProps::VaporPressure
{
    /**
     * @brief
     */
    class AmbroseWalton
    {
        double m_criticalTemperature = 0.0; /**< The critical temperature for the vapor pressure estimation. */
        double m_criticalPressure    = 0.0; /**< The critical pressure for the vapor pressure estimation. */
        double m_acentricFactor      = 0.0; /**< The acentric factor for the vapor pressure estimation. */

    public:
        /**
         * @brief Constructor, default
         * @note A default constructed VPAmbroseWalton object has all constants set to zero. Because the constants cannot
         * be modifed after object construction, a default constructed VPAmbroseWalton object is in itself of little use; the
         * main purpose is to serve as a placeholder for a correctly constructed object later on.
         */
        AmbroseWalton() = default;

        /**
         * @brief Constructor, taking critical properties and acentric factor as arguments.
         * @param criticalTemperature The critical temperature [K].
         * @param criticalPressure The critical pressure [Pa]
         * @param acentricFactor The acentric factor [-]
         */
        AmbroseWalton(double criticalTemperature, double criticalPressure, double acentricFactor)
            : m_criticalTemperature { criticalTemperature },
              m_criticalPressure { criticalPressure },
              m_acentricFactor { acentricFactor }
        {}

        /**
         * @brief Copy constructor
         */
        AmbroseWalton(const AmbroseWalton& other) = default;

        /**
         * @brief Move constructor
         */
        AmbroseWalton(AmbroseWalton&& other) noexcept = default;

        /**
         * @brief Destructor
         */
        ~AmbroseWalton() = default;

        /**
         * @brief Copy assignment operator
         */
        AmbroseWalton& operator=(const AmbroseWalton& other) = default;

        /**
         * @brief Move assignment operator
         */
        AmbroseWalton& operator=(AmbroseWalton&& other) noexcept = default;

        /**
         * @brief Function call operator, taking the temperature as an argument, and returns the vapor
         * pressure [Pa] at that temperature.
         * @param temperature The temperature [K] at which to calculate the vapor pressure.
         * @return The vapor pressure [Pa]
         * @warning If the object is default constructed only, operator() will yield NaN as the result.
         */
        double operator()(double temperature) const
        {
            using std::exp;
            using std::pow;
            auto tau = 1 - (temperature / m_criticalTemperature);
            auto f0  = (-5.97616 * tau + 1.29874 * pow(tau, 1.5) - 0.60394 * pow(tau, 2.5) - 1.06841 * pow(tau, 5)) / (1 - tau);
            auto f1  = (-5.03365 * tau + 1.11505 * pow(tau, 1.5) - 5.41217 * pow(tau, 2.5) - 7.46628 * pow(tau, 5)) / (1 - tau);
            auto f2  = (-0.64771 * tau + 2.41539 * pow(tau, 1.5) - 4.26979 * pow(tau, 2.5) + 3.25259 * pow(tau, 5)) / (1 - tau);

            return exp(f0 + m_acentricFactor * f1 + pow(m_acentricFactor, 2) * f2) * m_criticalPressure;
        }
    };
} // namespace PCProps::VaporPressure

#endif    // PCPROPS_AMBROSEWALTON_HPP
