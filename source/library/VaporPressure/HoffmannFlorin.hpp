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

#ifndef PCPROPS_HOFFMANNFLORIN_HPP
#define PCPROPS_HOFFMANNFLORIN_HPP

#include <array>
#include <cmath>

namespace PCProps::VaporPressure::detail
{
    double hfFunc(double temperature)
    {
        using std::log10;
        return (1.0 / temperature) - 7.9151E-3 + 2.6726E-3 * log10(temperature) - 0.8625E-6 * temperature;
    }
}    // namespace PCProps::VaporPressure::detail

namespace PCProps::VaporPressure
{
    /**
     * @brief
     */
    class HoffmannFlorin
    {
        std::array<double, 2> m_coefficients { 0.0, 0.0 };

    public:
        /**
         * @brief Constructor, default
         */
        HoffmannFlorin() = default;

        /**
         * @brief Constructor, taking temperature and vapor pressure for two reference points.
         * @param ref1Temp Reference 1 temperature [K]
         * @param ref1Psat Reference 1 vapor pressure [Pa]
         * @param ref2Temp Reference 2 temperature [K]
         * @param ref2Psat Reference 2 vapor pressure [Pa]
         */
        HoffmannFlorin(double ref1Temp, double ref1Psat, double ref2Temp, double ref2Psat)
            : m_coefficients { std::log(ref1Psat) - std::log(ref1Psat / ref2Psat) * detail::hfFunc(ref1Temp) / (detail::hfFunc(ref1Temp) - detail::hfFunc(ref2Temp)),
                               std::log(ref1Psat / ref2Psat) / (detail::hfFunc(ref1Temp) - detail::hfFunc(ref2Temp)) }
        {}

        /**
         * @brief Constructor, taking the two Hoffmann-Florin coefficients as arguments.
         * @param coeffA Coefficient A.
         * @param coeffB Coefficient B.
         */
        HoffmannFlorin(double coeffA, double coeffB) : m_coefficients { coeffA, coeffB } {}

        /**
         * @brief Copy constructor
         */
        HoffmannFlorin(const HoffmannFlorin& other) = default;

        /**
         * @brief Move constructor
         */
        HoffmannFlorin(HoffmannFlorin&& other) noexcept = default;

        /**
         * @brief Destructor
         */
        ~HoffmannFlorin() = default;

        /**
         * @brief Copy assignment operator
         */
        HoffmannFlorin& operator=(const HoffmannFlorin& other) = default;

        /**
         * @brief Move assignment operator
         */
        HoffmannFlorin& operator=(HoffmannFlorin&& other) noexcept = default;

        /**
         * @brief Function call operator, taking temperature [K] as an argument and returns the vapor pressure [Pa]
         * @param temperature The temperature at which to calculate the capor pressure [K]
         * @return The vapor pressure [Pa]
         */
        double operator()(double temperature) const
        {
            using std::exp;
            return exp(m_coefficients[0] + m_coefficients[1] * detail::hfFunc(temperature));
        }

        /**
         * @brief Getter for the Hoffman-Florin coefficients.
         * @return A std::array with the two coefficients.
         */
        std::array<double, 2> coefficients() const
        {
            return m_coefficients;
        }
    };

}    // namespace PCProps::VaporPressure
#endif    // PCPROPS_HOFFMANNFLORIN_HPP
