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

#ifndef PCPROPS_ALYLEE_HPP
#define PCPROPS_ALYLEE_HPP

#include <cmath>

namespace PCProps::HeatCapacity
{
    /**
     * @brief
     */
    class AlyLee
    {
        double m_A { 0.0 };
        double m_B { 0.0 };
        double m_C { 0.0 };
        double m_D { 0.0 };
        double m_E { 0.0 };

    public:

        /**
         * @brief
         */
        struct CreateFromDIPPR
        {
            double A;
            double B;
            double C;
            double D;
            double E;
        };

        /**
         * @brief
         */
        AlyLee() = default;

        /**
         * @brief
         * @param coeffA
         * @param coeffB
         * @param coeffC
         * @param coeffD
         * @param coeffE
         */
        AlyLee(double coeffA, double coeffB, double coeffC, double coeffD, double coeffE) : m_A(coeffA), m_B(coeffB), m_C(coeffC), m_D(coeffD), m_E(coeffE) {}

        /**
         * @brief
         * @param coefficients
         */
        explicit AlyLee(const CreateFromDIPPR& c) : AlyLee(c.A / 1000, c.B / 1000, c.C, c.D / 1000, c.E) {}

        /**
         * @brief
         * @param other
         */
        AlyLee(const AlyLee& other) =default;

        /**
         * @brief
         * @param other
         */
        AlyLee(AlyLee&& other) noexcept = default;

        /**
         * @brief
         */
        ~AlyLee() = default;

        /**
         * @brief
         * @param other
         * @return
         */
        AlyLee& operator=(const AlyLee& other) = default;

        /**
         * @brief
         * @param other
         * @return
         */
        AlyLee& operator=(AlyLee&& other) noexcept = default;

        /**
         * @brief
         * @param temperature
         * @return
         */
        double operator()(double temperature) const
        {
            using std::cosh;
            using std::pow;
            using std::sinh;

            return m_A + m_B * pow((m_C / temperature) / sinh(m_C / temperature), 2) + m_D * pow((m_E / temperature) / cosh(m_E / temperature), 2);
        }
    };

}    // namespace PCProps::HeatCapacity

#endif    // PCPROPS_ALYLEE_HPP
