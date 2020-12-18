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

#ifndef PCPROPS_VPRIEDEL_HPP
#define PCPROPS_VPRIEDEL_HPP

#include <string>
#include <array>

namespace PCProps::VaporPressure
{

    enum class VPRiedelType { Organic, Acid, Alcohol};

    /**
     * @brief
     */
    class VPRiedel
    {

        double m_criticalTemperature = 0.0;
        double m_criticalPressure = 0.0;

        std::array<double, 4> m_coefficients {0.0, 0.0, 0.0, 0.0};

    public:

        /**
         * @brief Constructor, default, with no parameters. Calling the operator() on a default constructed object
         * will yield NaN as the result. To turn it into a valid object, use the move or copy assignment operator.
         */
        VPRiedel();

        /**
         * @brief Constructor, taking the normal boiling temperature [K], the critical temperature [K]
         * and the critical pressure [Pa] as arguments
         * @param boilingTemperature The normal boiling temperature [K] (i.e. at 101325 pascals).
         * @param criticalTemperature The critical temperature of the fluid [K].
         * @param criticalPressure The critical pressure of the fluid [Pa].
         * @param type The type of component; can be acid, alcohol or general organic compound.
         */
        VPRiedel(double boilingTemperature, double criticalTemperature, double criticalPressure, VPRiedelType type = VPRiedelType::Organic);

        /**
         * @brief
         * @param coeffA
         * @param coeffB
         * @param coeffC
         * @param coeffD
         */
        VPRiedel(double criticalTemperature, double criticalPressure, double coeffA, double coeffB, double coeffC, double coeffD);

        /**
         * @brief Destructor
         */
        ~VPRiedel();

        /**
         * @brief Copy constructor
         */
        VPRiedel(const VPRiedel& other);

        /**
         * @brief Move constructor
         */
        VPRiedel(VPRiedel&& other) noexcept;

        /**
         * @brief Copy assignment operator
         */
        VPRiedel& operator=(const VPRiedel& other);

        /**
         * @brief Move assignment operator
         */
        VPRiedel& operator=(VPRiedel&& other) noexcept;

        /**
         * @brief operator(), yielding the saturation pressure at the requested temperature for the fluid.
         * @param temperature The temperature [K] at which to get the vapor pressure.
         * @return The vapor pressure [Pa]
         * @warning If the object is default constructed only, operator() will yield NaN as the result.
         * To make the object valid, call the reset member function.
         */
        double operator()(double temperature) const;

        /**
         * @brief Get the critical temperature used for the vapor pressure estimation.
         * @return The critical temperature [K]
         */
        double criticalTemperature() const;

        /**
         * @brief Get the critical pressure used for the vapor pressure estimation.
         * @return The critical pressure [Pa]
         */
        double criticalPressure() const;

        /**
         * @brief Get the Riedel equation coefficients.
         * @return A std::array with the four Riedel coefficients.
         */
        std::array<double, 4> coefficients() const;

    };
} // namespace PCProps::VaporPressure

#endif    // PCPROPS_VPRIEDEL_HPP
