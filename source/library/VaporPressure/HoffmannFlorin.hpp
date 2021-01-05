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
        HoffmannFlorin();

        /**
         * @brief Constructor, taking temperature and vapor pressure for two reference points.
         * @param ref1Temp Reference 1 temperature [K]
         * @param ref1Psat Reference 1 vapor pressure [Pa]
         * @param ref2Temp Reference 2 temperature [K]
         * @param ref2Psat Reference 2 vapor pressure [Pa]
         */
        HoffmannFlorin(double ref1Temp, double ref1Psat, double ref2Temp, double ref2Psat);

        /**
         * @brief Constructor, taking the two Hoffmann-Florin coefficients as arguments.
         * @param coeffA Coefficient A.
         * @param coeffB Coefficient B.
         */
        HoffmannFlorin(double coeffA, double coeffB);

        /**
         * @brief Copy constructor
         */
        HoffmannFlorin(const HoffmannFlorin& other);

        /**
         * @brief Move constructor
         */
        HoffmannFlorin(HoffmannFlorin&& other) noexcept;

        /**
         * @brief Destructor
         */
        ~HoffmannFlorin();

        /**
         * @brief Copy assignment operator
         */
        HoffmannFlorin& operator=(const HoffmannFlorin& other);

        /**
         * @brief Move assignment operator
         */
        HoffmannFlorin& operator=(HoffmannFlorin&& other) noexcept;

        /**
         * @brief Function call operator, taking temperature [K] as an argument and returns the vapor pressure [Pa]
         * @param temperature The temperature at which to calculate the capor pressure [K]
         * @return The vapor pressure [Pa]
         */
        double operator()(double temperature) const;

        /**
         * @brief Getter for the Hoffman-Florin coefficients.
         * @return A std::array with the two coefficients.
         */
        std::array<double, 2> coefficients() const;
    };

}    // namespace PCProps::VaporPressure
#endif    // PCPROPS_HOFFMANNFLORIN_HPP
