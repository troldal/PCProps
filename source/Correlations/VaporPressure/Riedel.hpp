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

#ifndef PCPROPS_RIEDEL_HPP
#define PCPROPS_RIEDEL_HPP

#include <array>
#include <cmath>
#include <string>

namespace PCProps::VaporPressure
{
    enum class VPRiedelType { Organic, Acid, Alcohol };

    /**
     * @brief
     */
    class Riedel
    {
        double m_criticalTemperature = 0.0;
        double m_criticalPressure    = 0.0;

        std::array<double, 4> m_coefficients { 0.0, 0.0, 0.0, 0.0 };

    public:
        /**
         * @brief Constructor, default, with no parameters. Calling the operator() on a default constructed object
         * will yield NaN as the result. To turn it into a valid object, use the move or copy assignment operator.
         */
        Riedel() = default;

        /**
         * @brief Constructor, taking the normal boiling temperature [K], the critical temperature [K]
         * and the critical pressure [Pa] as arguments
         * @param boilingTemperature The normal boiling temperature [K] (i.e. at 101325 pascals).
         * @param criticalTemperature The critical temperature of the fluid [K].
         * @param criticalPressure The critical pressure of the fluid [Pa].
         * @param type The type of component; can be acid, alcohol or general organic compound.
         */
        Riedel(double boilingTemperature, double criticalTemperature, double criticalPressure, VPRiedelType type = VPRiedelType::Organic)
            : m_criticalTemperature { criticalTemperature },
              m_criticalPressure { criticalPressure }
        {
            double tbr = boilingTemperature / criticalTemperature;
            double h = tbr * std::log(criticalPressure/101325.0)/(1.0 - tbr);
            double K = [&]() {
                   switch (type) {
                       case VPRiedelType::Organic:
                           return 0.0838;

                       case VPRiedelType::Acid:
                           return -0.12 + 0.025 * h;

                       case VPRiedelType::Alcohol:
                           return 0.373 - 0.030 * h;

                       default:
                           return 0.0;
                   }
            }();

            double psi     = -35.0 + 36.0 / tbr + 42.0 * std::log(tbr) - std::pow(tbr, 6);
            double alpha_c = (3.758 * K * psi + std::log(criticalPressure / 101325.0)) / (K * psi - std::log(tbr));

            m_coefficients[3] = K * (alpha_c - 3.758);
            m_coefficients[2] = alpha_c - 42.0 * m_coefficients[3];
            m_coefficients[1] = -36.0 * m_coefficients[3];
            m_coefficients[0] = 35.0 * m_coefficients[3];
        }

        /**
         * @brief
         * @param criticalTemperature
         * @param criticalPressure
         * @param coeffA
         * @param coeffB
         * @param coeffC
         * @param coeffD
         */
        Riedel(double criticalTemperature, double criticalPressure, double coeffA, double coeffB, double coeffC, double coeffD)
            : m_criticalTemperature { criticalTemperature },
              m_criticalPressure { criticalPressure },
              m_coefficients { coeffA, coeffB, coeffC, coeffD }
        {}

        /**
         * @brief Destructor
         */
        ~Riedel() = default;

        /**
         * @brief Copy constructor
         */
        Riedel(const Riedel& other) = default;

        /**
         * @brief Move constructor
         */
        Riedel(Riedel&& other) noexcept = default;

        /**
         * @brief Copy assignment operator
         */
        Riedel& operator=(const Riedel& other) = default;

        /**
         * @brief Move assignment operator
         */
        Riedel& operator=(Riedel&& other) noexcept = default;

        /**
         * @brief operator(), yielding the saturation pressure at the requested temperature for the fluid.
         * @param temperature The temperature [K] at which to get the vapor pressure.
         * @return The vapor pressure [Pa]
         * @warning If the object is default constructed only, operator() will yield NaN as the result.
         * To make the object valid, call the reset member function.
         */
        double operator()(double temperature) const
        {
            using std::log;
            using std::pow;
            auto tr = temperature / m_criticalTemperature;
            return exp(m_coefficients[0] + m_coefficients[1] / tr + m_coefficients[2] * log(tr) + m_coefficients[3] * pow(tr, 6)) * m_criticalPressure;
        }
    };
} // namespace PCProps::SaturationPressure

#endif    // PCPROPS_RIEDEL_HPP
