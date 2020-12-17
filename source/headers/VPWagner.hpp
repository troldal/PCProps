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

#ifndef PCPROPS_VPWAGNER_HPP
#define PCPROPS_VPWAGNER_HPP

#include <array>

namespace PCProps::VaporPressure
{
    enum class VPWagnerForm {Form25, Form36};

    /**
     * @brief The VPWagner class encapsulates the Wagner equation(s) for calculating vapor pressures.
     * @details The VPWagner class encapsulates the Wagner equation(s) for calculating vapor pressures:
     * \f[ ln \ \frac {P_{s}}{P_{c}} = \frac {1}{T_{r}}\left [ A\cdot(1-T_{r}) + B\cdot(1-T_{r})^{1.5}+C\cdot(1-T_{r})^{3}+D\cdot(1-T_{r})^6 \right] \f]
     * and:
     * \f[ ln \ \frac {P_{s}}{P_{c}} = \frac {1}{T_{r}}\left [ A\cdot(1-T_{r}) + B\cdot(1-T_{r})^{1.5}+C\cdot(1-T_{r})^{2.5}+D\cdot(1-T_{r})^5 \right] \f]
     *
     * The two forms are called the 3-6 form and the 2.5-5 form, respectively. Apparently, the 2.5-5 form is slightly more accurate,
     * and some authors prefers this. The VPWagner class supports both forms. The default is the 2.5-5 form, but the 3-6 form
     * can be chosen by using the VPWagnerForm::Form36 flag in the constructor.
     *
     * The coefficients and critical properties can only be set in the constructor;
     * if coefficients needs to be changed after construction, a new object has to be created.
     *
     * Note that operator() returns Psat, the vapor pressure [Pa], not \f$ln \ \frac {P_{sat}}{P_{c}} \f$.
     *
     * All pressures are in Pascals and all temperatures in Kelvin.
     */
    class VPWagner
    {
        // ===== Private Data Members

        double m_criticalTemperature = 0.0;
        double m_criticalPressure = 0.0;
        std::array<double, 4> m_coefficients { 0.0, 0.0, 0.0, 0.0};

        double m_expC = 2.5;
        double m_expD = 5.0;

    public:

        /**
         * @brief Default constructor.
         * @note A default constructed VPWagner object has all coefficients and critical properties set to zero.
         * Because the coefficients cannot be modifed after object construction, a default constructed VPWagner object
         * is in itself of little use; indeed, using the function call operator on a default constructed object will
         * result in division-by-zero and will return NaN. The main purpose is to serve as a placeholder for a correctly
         * constructed object later on.
         */
        VPWagner();

        /**
         * @brief Constructor, taking three or more coefficients as input.
         * @param criticalTemperature The critical temperature of the component [K].
         * @param criticalPressure The critical pressure of the component [Pa].
         * @param A Coefficient A
         * @param B Coefficient B
         * @param C Coefficient C
         * @param D Coefficient D
         * @note Please note that the VPWagner class assumes temperatures in Kelvin, pressures in Pascals and that
         * the equation uses natural logarithms (rather than base-10 logarithm). If coefficients exists with a different
         * basis, they will have to be converted first.
         */
        VPWagner(double criticalTemperature, double criticalPressure, double A, double B, double C, double D, VPWagnerForm form = VPWagnerForm::Form25);

        /**
         * @brief Copy constructor.
         */
        VPWagner(const VPWagner& other);

        /**
         * @brief Move constructor.
         */
        VPWagner(VPWagner&& other) noexcept;

        /**
         * @brief Destructor.
         */
        ~VPWagner();

        /**
         * @brief Copy assignment operator.
         */
        VPWagner& operator=(const VPWagner& other);

        /**
         * @brief Move assignment operator.
         */
        VPWagner& operator=(VPWagner&& other) noexcept;

        /**
         * @brief Function call operator, yielding the saturation pressure [Pa] at the requested temperature [K] for the component.
         * @param temperature The temperature [K] at which to get the vapor pressure.
         * @return The vapor pressure [Pa].
         * @warning If the object is default constructed only, operator() will return NaN.
         */
        double operator()(double temperature) const;

        /**
         * @brief Get the critical temperature used in the correlation.
         * @return The critical temperature [K].
         */
        double criticalTemperature() const;

        /**
         * @brief Get the critical pressure used in the correlation.
         * @return The critical temperature [Pa].
         */
        double criticalPressure() const;

        /**
         * @brief Get the coefficients of the current object.
         * @return An array with the coefficients A-D.
         */
        std::array<double, 4> coefficients() const;

        /**
         * @brief Get the form of the Wagner equation (2.5-5 or 3-6)
         * @return A VPWagnerForm object, representing the Wagner form.
         */
        VPWagnerForm form() const;
    };

}    // namespace PCProps::VaporPressure

#endif    // PCPROPS_VPWAGNER_HPP
