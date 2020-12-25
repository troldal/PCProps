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

#ifndef PCPROPS_SLVYENWOODS_HPP
#define PCPROPS_SLVYENWOODS_HPP

namespace PCProps::LiquidVolume
{
    /**
     * @brief The SLVYenWoods class encapsulates the Yen-Woods equation for calculating the saturated molar volume for liquids.
     * @details The SLVYenWoods class encapsulates the Yen-Woods equation for calculating the saturated molar volume for liquids:
     * \f[ \frac{V_{c}}{V_{s}} = 1
     * + A \cdot \left(1 - \frac{T}{T_{c}} \right)^{\frac{1}{3}}
     * + B \cdot \left(1 - \frac{T}{T_{c}} \right)^{\frac{2}{3}}
     * + C \cdot \left(1 - \frac{T}{T_{c}} \right)
     * + D \cdot \left(1 - \frac{T}{T_{c}} \right)^{\frac{4}{3}} \f]
     *
     * The above equation is the original form of the Yen-Woods equation. A modified form of this equation uses an exponent
     * of 0.35 in the second term:
     * \f[ \frac{V_{c}}{V_{s}} = 1
     * + A \cdot \left(1 - \frac{T}{T_{c}} \right)^{0.35}
     * + B \cdot \left(1 - \frac{T}{T_{c}} \right)^{\frac{2}{3}}
     * + C \cdot \left(1 - \frac{T}{T_{c}} \right)
     * + D \cdot \left(1 - \frac{T}{T_{c}} \right)^{\frac{4}{3}} \f]
     *
     * PPDS and DIPPR (equation 116) is based on the modified Yen-Woods equation,
     * but uses a form that looks slightly different:
     * \f[ \rho_{s} = \rho_{c}
     * + A \cdot \left(1 - \frac{T}{T_{c}} \right)^{0.35}
     * + B \cdot \left(1 - \frac{T}{T_{c}} \right)^{\frac{2}{3}}
     * + C \cdot \left(1 - \frac{T}{T_{c}} \right)
     * + D \cdot \left(1 - \frac{T}{T_{c}} \right)^{\frac{4}{3}} \f]
     *
     * This equation is explicit in density rather than molar volume. Other than that, the first term in the modified
     * Yen-Woods equation (a constant of 1.0) has been replaced by the critical density. This implies that the four
     * coefficients have a unit, rather than being dimensionless. Therefore, to transform this form to the modified
     * Yen-Woods equation, the coefficients will have to be converted.
     * Assuming that the density is in terms of molar density [mol/m3], the equation can be rearranged to
     * follow the form of the modified Yen-Woods equation:
     *
     * \f[ \frac{V_{c}}{V_{s}} = 1 +
     * V_{c} \cdot A \cdot \left(1 - \frac{T}{T_{c}} \right)^{0.35} +
     * V_{c} \cdot B \cdot \left(1 - \frac{T}{T_{c}} \right)^{\frac{2}{3}} +
     * V_{c} \cdot C \cdot \left(1 - \frac{T}{T_{c}} \right) +
     * V_{c} \cdot D \cdot \left(1 - \frac{T}{T_{c}} \right)^{\frac{4}{3}} \f]
     *
     * This shows that the coefficients will have to multiplied by the critical volume in order to be compatible
     * with the modified Yen-Woods equation. If the density unit is different (e.g. kg/m3 or mol/dm3), then the coefficients
     * have to be converted accordingly.
     *
     * To simplify the creation of SLVYenWoods objects from coefficients from PPDS or DIPPR, static factory functions
     * are provided, that will perform the coefficient conversion during object construction.
     *
     * The return value of the function call operator will be the saturated molar volume in m3/mol. The Yen-Woods
     * equations in it's various forms traditionally calculates density, but for the sake of consistency, the
     * SLVYenWoods will calculate the molar volume instead.
     *
     * An SLVYenWoods object can be created by supplying the four Yen-Woods coefficients to the constructor. The coefficients
     * have to correspond to the original form of the equation and must be based on molar volume [m3/mol]. It is recommended
     * to use the static factory functions to convert coefficients from known data sources, or estimate them from
     * critical properties.
     *
     * Note that the coefficients can only be set in the constructor or using the factory functions;
     * if coefficients needs to be changed after construction, a new object has to be created.
     *
     */
    class SLVYenWoods final
    {
    public:
        /**
         * @brief A nested enum class, defining which form of the Yen-Woods equation to use.
         */
        enum class Form { Original, Modified };

    private:
        Form m_type = Form::Modified; /**< The form of the equation used (original or modified). */

        double m_criticalTemperature = 0.0; /**< The critical temperature [K] of the component. */
        double m_criticalVolume      = 0.0; /**< The critical volume [m3/mol] of the component. */
        double m_A                   = 0.0; /**< The Yen-Woods coefficient A */
        double m_B                   = 0.0; /**< The Yen-Woods coefficient B */
        double m_C                   = 0.0; /**< The Yen-Woods coefficient C */
        double m_D                   = 0.0; /**< The Yen-Woods coefficient D */

        /**
         * @brief Constructor, taking critical properties and coefficients A-D as arguments.
         * @details This constructor is declared private, as the static factory functions should be used for object creation.
         * @param criticalTemperature The critical temperature [K]
         * @param criticalVolume The critical volume [m3/mol]
         * @param coeffA The A coefficient of the original Yen-Woods equation.
         * @param coeffB The B coefficient of the original Yen-Woods equation.
         * @param coeffC The C coefficient of the original Yen-Woods equation.
         * @param coeffD The D coefficient of the original Yen-Woods equation.
         * @param type The form of the equation used (original or modified).
         */
        SLVYenWoods(
            double            criticalTemperature,
            double            criticalVolume,
            double            coeffA,
            double            coeffB,
            double            coeffC,
            double            coeffD,
            SLVYenWoods::Form type = Form::Modified);

    public:
        /**
         * @brief Copy constructor
         */
        SLVYenWoods(const SLVYenWoods& other);

        /**
         * @brief Move constructor
         */
        SLVYenWoods(SLVYenWoods&& other) noexcept;

        /**
         * @brief Destructor
         */
        ~SLVYenWoods();

        /**
         * @brief Copy assignment operator
         */
        SLVYenWoods& operator=(const SLVYenWoods& other);

        /**
         * @brief Move assignment operator
         */
        SLVYenWoods& operator=(SLVYenWoods&& other) noexcept;

        /**
         * @brief Function call operator, taking temperature [K] as an argument and returns the liquid molar volume [m3/mol]
         * @param temperature The temperature [K]
         * @return The saturated liquid molar volume [m3/mol]
         */
        double operator()(double temperature);

        /**
         * @brief Factory function for creating an SLVYenWoods object using Original Yen-Woods coefficients.
         * @details This function creates an SLVYenWoods object using coefficients corresponding to the
         * original Yen-Woods equation:
         * \f[ \frac{V_{c}}{V_{s}} = 1
         * + A \cdot \left(1 - \frac{T}{T_{c}} \right)^{\frac{1}{3}}
         * + B \cdot \left(1 - \frac{T}{T_{c}} \right)^{\frac{2}{3}}
         * + C \cdot \left(1 - \frac{T}{T_{c}} \right)
         * + D \cdot \left(1 - \frac{T}{T_{c}} \right)^{\frac{4}{3}} \f]
         * The coefficients used should correspond to molar volume with the unit of m3/mol
         * @param criticalTemperature The critical temperature [K]
         * @param criticalVolume The critical volume [m3/mol]
         * @param coeffA The Original Yen-Woods A coefficient.
         * @param coeffB The Original Yen-Woods B coefficient.
         * @param coeffC The Original Yen-Woods C coefficient.
         * @param coeffD The Original Yen-Woods D coefficient.
         * @return A SLVYenWoods object constructed from Original Yen-Woods coefficients.
         * @note This function will create an object corresponding to the original form of the Yen-Woods equation.
         */
        static SLVYenWoods createFromOriginalYenWoodsCoefficients(
            double criticalTemperature,
            double criticalVolume,
            double coeffA,
            double coeffB,
            double coeffC,
            double coeffD);

        /**
         * @brief Factory function for creating an SLVYenWoods object using Modified Yen-Woods coefficients.
         * @details This function creates an SLVYenWoods object using coefficients corresponding to the
         * modified Yen-Woods equation:
         * \f[ \frac{V_{c}}{V_{s}} = 1
         * + A \cdot \left(1 - \frac{T}{T_{c}} \right)^{0.35}
         * + B \cdot \left(1 - \frac{T}{T_{c}} \right)^{\frac{2}{3}}
         * + C \cdot \left(1 - \frac{T}{T_{c}} \right)
         * + D \cdot \left(1 - \frac{T}{T_{c}} \right)^{\frac{4}{3}} \f]
         * The coefficients used should correspond to molar volume with the unit of m3/mol
         * @param criticalTemperature The critical temperature [K]
         * @param criticalVolume The critical volume [m3/mol]
         * @param coeffA The Modified Yen-Woods A coefficient.
         * @param coeffB The Modified Yen-Woods B coefficient.
         * @param coeffC The Modified Yen-Woods C coefficient.
         * @param coeffD The Modified Yen-Woods D coefficient.
         * @return A SLVYenWoods object constructed from Modified Yen-Woods coefficients.
         * @note This function will create an object corresponding to the modified form of the Yen-Woods equation.
         */
        static SLVYenWoods createFromModifiedYenWoodsCoefficients(
            double criticalTemperature,
            double criticalVolume,
            double coeffA,
            double coeffB,
            double coeffC,
            double coeffD);

        /**
         * @brief Factory function for creating an SLVYenWoods object using PPDS coefficients.
         * @details The PPDS equation follows the same form as the (modified) Yen-Woods equation, but is used to calculate
         * saturated liquid density [kg/m3], rather than molar volume [m3/mol]. Coefficients for the PPDS equation can be
         * found in a number of references, notably the VDI Heat Atlas. This factory function converts the PPDS coefficients to
         * the Yen-Woods form. This is done simply by multiplying the coefficients by the term \f$ \frac {V_{c} \cdot 1000}{M_{w}} \f$
         * @param criticalTemperature The critical temperature [K]
         * @param criticalVolume The critical volume [m3/mol]
         * @param molecularWeight The molecular weight [kg/kmol]
         * @param coeffA The PPDS A coefficient.
         * @param coeffB The PPDS B coefficient.
         * @param coeffC The PPDS C coefficient.
         * @param coeffD The PPDS D coefficient.
         * @return A SLVYenWoods object constructed from PPDS coefficients.
         * @note This function will create an object corresponding to the modified form of the Yen-Woods equation.
         */
        static SLVYenWoods createFromPPDSCoefficients(
            double criticalTemperature,
            double criticalVolume,
            double molecularWeight,
            double coeffA,
            double coeffB,
            double coeffC,
            double coeffD);

        /**
         * @brief Factory function for creating an SLVYenWoods object using DIPPR 116 coefficients.
         * @details The DIPPR 116 equation follows the same form as the (modified) Yen-Woods equation, but is used to calculate
         * saturated liquid density [mol/dm3], rather than molar volume [m3/mol]. Coefficients for the DIPPR 116 equation is
         * only available for very few substances, notably for water. This factory function converts the DIPPR 116 coefficients to
         * the Yen-Woods form. This is done simply by multiplying the coefficients by the term \f$ {V_{c}} \cdot 1000 \f$
         * @param criticalTemperature The critical temperature [K]
         * @param criticalVolume The critical volume [m3/mol]
         * @param coeffA The DIPPR 116 A coefficient.
         * @param coeffB The DIPPR 116 B coefficient.
         * @param coeffC The DIPPR 116 C coefficient.
         * @param coeffD The DIPPR 116 D coefficient.
         * @return A SLVYenWoods object constructed from DIPPR 116 coefficients.
         * @note This function will create an object corresponding to the modified form of the Yen-Woods equation.
         */
        static SLVYenWoods createFromDIPPR116Coefficients(
            double criticalTemperature,
            double criticalVolume,
            double coeffA,
            double coeffB,
            double coeffC,
            double coeffD);

        /**
         * @brief Factory function for creating an SLVYenWoods object using the original correlations in the paper
         * by Yen and Woods. This should only be used if coefficients are not directly available.
         * @details In the original paper by Yen and Woods, the authors proposes a set of simple correlations, that can
         * be used to estimate the four coefficients from the critical compressibility, \f$ {Z_{c}} \f$, in case
         * coefficients fitted to experimental data are not available. The coefficients are estimated as follows:
         * \f[ A = 17.4425 - 214.578 \cdot Z_{z} + 989.625 \cdot Z_{c}^2 - 1522.06 \cdot Z_{c}^3 \f]
         * \f[ B = -3.28257 + 13.6377 \cdot Z_{z} + 107.4844 \cdot Z_{c}^2 - 384.211 \cdot Z_{c}^3 \ \textrm{(if} \ Z_{c} \leq \textrm{
         * 0.26)} \f] \f[ B = 60.2091 - 402.063 \cdot Z_{z} + 501.0 \cdot Z_{c}^2 + 641.0 \cdot Z_{c}^3 \ \textrm{(if} \ Z_{c} > \textrm{
         * 0.26)} \f] \f[ C = 0 \f] \f[ D = 0.93 - B \f]
         * @param criticalTemperature The critical temperature [K]
         * @param criticalVolume The critical volume [m3/mol]
         * @param criticalCompressibility The critical compressibility [-]
         * @return A SLVYenWoods object constructed from Yen and Woods estimation procedure.
         * @note This function will create an object corresponding to the original form of the Yen-Woods equation.
         */
        static SLVYenWoods createFromYenWoodsEstimation(double criticalTemperature, double criticalVolume, double criticalCompressibility);
    };

}    // namespace PCProps::LiquidVolume

#endif    // PCPROPS_SLVYENWOODS_HPP
