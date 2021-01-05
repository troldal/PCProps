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
        AmbroseWalton();

        /**
         * @brief Constructor, taking critical properties and acentric factor as arguments.
         * @param criticalTemperature The critical temperature [K].
         * @param criticalPressure The critical pressure [Pa]
         * @param acentricFactor The acentric factor [-]
         */
        AmbroseWalton(double criticalTemperature, double criticalPressure, double acentricFactor);

        /**
         * @brief Copy constructor
         */
        AmbroseWalton(const AmbroseWalton& other);

        /**
         * @brief Move constructor
         */
        AmbroseWalton(AmbroseWalton&& other) noexcept;

        /**
         * @brief Destructor
         */
        ~AmbroseWalton();

        /**
         * @brief Copy assignment operator
         */
        AmbroseWalton& operator=(const AmbroseWalton& other);

        /**
         * @brief Move assignment operator
         */
        AmbroseWalton& operator=(AmbroseWalton&& other) noexcept;

        /**
         * @brief Function call operator, taking the temperature as an argument, and returns the vapor
         * pressure [Pa] at that temperature.
         * @param temperature The temperature [K] at which to calculate the vapor pressure.
         * @return The vapor pressure [Pa]
         * @warning If the object is default constructed only, operator() will yield NaN as the result.
         */
        double operator()(double temperature) const;

        /**
         * @brief Getter for the critical temperature used for the vapor pressure estimation.
         * @return The critical temperature [K]
         */
        double criticalTemperature() const;

        /**
         * @brief Getter for the critical pressure used for the vapor pressure estimation.
         * @return The critical pressure [Pa]
         */
        double criticalPressure() const;

        /**
         * @brief Getter for the critical pressure used for the vapor pressure estimation.
         * @return The acentric factor [-]
         */
        double acentricFactor() const;

    };
} // namespace PCProps::VaporPressure

#endif    // PCPROPS_AMBROSEWALTON_HPP
