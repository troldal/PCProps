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

#ifndef PCPROPS_ANTOINEEXTENDED_HPP
#define PCPROPS_ANTOINEEXTENDED_HPP

#include <array>
#include <cmath>

namespace PCProps::VaporPressure
{
    /**
     * @brief The VPAntoineExt class encapsulates the extended Antoine equation for calculating vapor pressures.
     * @details The VPAntoineExt class encapsulates the extended Antoine equation for calculating vapor pressures:
     * \f[ ln \ P_{sat} = A + \frac{B}{T+C}+D \cdot T + E \cdot ln \ T + F \cdot T^G \f]
     * Essentially, this equation is just a collection of useful terms. The first two terms represents the
     * original Antoine equation. Hence, the original Antoine equation can be achieved by setting the coefficients
     * D-G to zero. Similarly, other vapor pressure equations, such as DIPPR, uses a form similar to the extended Antoine
     * equation.
     *
     * The coefficients can only be set in the constructor; if coefficients needs to be changed after construction,
     * a new object has to be created.
     *
     * Note that operator() returns Psat, the vapor pressure [Pa], not \f$ln \ P_{sat} \f$.
     */
    class AntoineExtended
    {
        // ===== Private Data Members

        double m_coeffA {};
        double m_coeffB {};
        double m_coeffC {};
        double m_coeffD {};
        double m_coeffE {};
        double m_coeffF {};
        double m_coeffG {};

    public:
        struct CreateFromDIPPR
        {
            double A;
            double B;
            double C;
            double D;
            double E;
        };

        struct CreateFromYaws
        {
            double A;
            double B;
            double C;
            double D;
            double E;
        };

        /**
         * @brief Default constructor.
         * @note A default constructed VPAntoineExt object has all coefficients set to zero. Because the coefficients cannot
         * be modifed after object construction, a default constructed VPAntoineExt object is in itself of little use; the
         * main purpose is to serve as a placeholder for a correctly constructed object later on.
         */
        AntoineExtended() = default;

        /**
         * @brief Constructor, taking three or more coefficients as input.
         * @param A Coefficient A
         * @param B Coefficient B
         * @param C Coefficient C
         * @param D Coefficient D
         * @param E Coefficient E
         * @param F Coefficient F
         * @param G Coefficient G
         * @note Please note that the VPAntoineExt class assumes temperatures in Kelvin, pressures in Pascals and that
         * the equation uses natural logarithms (rather than base-10 logarithm). If coefficients exists with a different
         * basis, they will have to be converted first.
         */
        AntoineExtended(double A, double B, double C, double D = 0.0, double E = 0.0, double F = 0.0, double G = 0.0)
            : m_coeffA { A },
              m_coeffB { B },
              m_coeffC { C },
              m_coeffD { D },
              m_coeffE { E },
              m_coeffF { F },
              m_coeffG { G }
        {}

        /**
         * @brief Constructor, taking a CreateFromDIPPR struct, for creating an object using DIPPR coefficients (eg. from Perry)
         * @param coefficients A CreateFromDIPPR with the DIPPR coefficients.
         */
        explicit AntoineExtended(const CreateFromDIPPR& coefficients)
            : AntoineExtended { coefficients.A, coefficients.B, 0.0, 0.0, coefficients.C, coefficients.D, coefficients.E }
        {}

        /**
         * @brief Constructor, taking a CreateFromYaws struct, for creating an object using coefficients from Yaws handbooks.
         * @param coefficients A CreateFromYaws with the Yaws coefficients.
         */
        explicit AntoineExtended(const CreateFromYaws& c)
            : AntoineExtended(log(133.322368) + c.A * log(10), c.B * log(10), 0.0, c.D * log(10), c.C, c.E * log(10), 2)
        {}


        /**
         * @brief Copy constructor.
         */
        AntoineExtended(const AntoineExtended& other) = default;

        /**
         * @brief Move constructor.
         */
        AntoineExtended(AntoineExtended&& other) noexcept = default;

        /**
         * @brief Destructor.
         */
        ~AntoineExtended() = default;

        /**
         * @brief Copy assignment operator.
         */
        AntoineExtended& operator=(const AntoineExtended& other) = default;

        /**
         * @brief Move assignment operator.
         */
        AntoineExtended& operator=(AntoineExtended&& other) noexcept = default;

        /**
         * @brief Function call operator, yielding the saturation pressure [Pa] at the requested temperature [K] for the component.
         * @param temperature The temperature [K] at which to get the vapor pressure.
         * @return The vapor pressure [Pa].
         * @warning If the object is default constructed only, operator() will yield 1.0 as the result.
         */
        double operator()(double temperature) const
        {
            using std::exp;
            using std::log;
            using std::pow;
            return exp(m_coeffA + m_coeffB / (temperature + m_coeffC) + m_coeffD * temperature + m_coeffE * log(temperature) + m_coeffF * pow(temperature, static_cast<int>(m_coeffG)));
        }

    };

}    // namespace PCProps::VaporPressure

#endif    // PCPROPS_ANTOINEEXTENDED_HPP
