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

#ifndef PCPROPS_RACKETT_HPP
#define PCPROPS_RACKETT_HPP

namespace PCProps::LiquidVolume
{
    /**
     * @brief The SLVRackett class encapsulates the Rackett equation for calculating the saturated molar volume for liquids.
     * @details The SLVRackett class encapsulates the Rackett equation for calculating the saturated molar volume for liquids:
     * \f[ V_{s} = A \cdot {B^{\left [1+(1-T/C)^D\right]}} \f]
     *
     * The return value of the function call operator will be the saturated molar volume in m3/mol. Some references
     * expresses the Rackett equation in terms of the molar density [mol/m3] as the return value, which is simply the
     * reciprocal of the expression used in this class.
     *
     * An SLVRackett object can be created by supplying the four Rackett coefficients to the constructor. Coefficients
     * can be found in various sources, e.g. Perry's Chemical Engineer's Handbook, which tabulates coefficients from
     * the DIPPR 801 database (for molar density rather than molar volume).
     * Remember, however, that the units may need to be converted; this can be done simply by applying a conversion
     * factor to the 'A' coefficient.
     *
     * If coefficients are not available, they can be estimated from critical properties, acentric factor or a known
     * reference point. To construct an SLVRackett object from estimated coefficients, a number of static factory
     * functions can be used. Please refer to the documentation of those functions for further details.
     *
     * Note that the coefficients can only be set in the constructor or using the factory functions;
     * if coefficients needs to be changed after construction, a new object has to be created.
     *
     */
    class Rackett
    {
        double m_A = 0.0;
        double m_B = 0.0;
        double m_C = 0.0;
        double m_D = 0.0;

    public:
        // ===== Constructors & Destructors ===== //

        struct CreateFromDIPPR
        {
            double A;
            double B;
            double C;
            double D;
        };

        struct CreateFromYaws
        {
            double A;
            double B;
            double Tc;
            double n;
            double molecularWeight;
        };

        /**
         * @brief Constructor, default. Sets all coefficients to zero.
         * @warning Calling the function call operator on a default constructed SLVRackett object will result in
         * division by zero and return NaN.
         */
        Rackett();

        /**
         * @brief Constructor, taking the four Rackett coefficients as arguments.
         * @param coeffA Rackett coefficient A.
         * @param coeffB Rackett coefficient B.
         * @param coeffC Rackett coefficient C.
         * @param coeffD Rackett coefficient D.
         * @note The coefficients must be based on an temperature input in Kelvin, and a density return value
         * in m3/mol.
         */
        Rackett(double coeffA, double coeffB, double coeffC, double coeffD);

        /**
         * @brief
         * @param coefficients
         */
        Rackett(const CreateFromDIPPR& coefficients);

        /**
         * @brief
         * @param coefficients
         */
        Rackett(const CreateFromYaws& coefficients);

        /**
         * @brief Copy constructor.
         */
        Rackett(const Rackett& other);

        /**
         * @brief Move constructor.
         */
        Rackett(Rackett&& other) noexcept;

        /**
         * @brief Destructor.
         */
        ~Rackett();

        // ===== Manipulators ===== //

        /**
         * @brief Copy assignment operator.
         */
        Rackett& operator=(const Rackett& other);

        /**
         * @brief Move assignment operator.
         */
        Rackett& operator=(Rackett&& other) noexcept;

        // ===== Accessors ===== //

        /**
         * @brief Function call operator. This function is used to calculate the saturated liquid volume [m3/mol] at a
         * given temperature [K]
         * @param temperature The temperature [K] at which to evaluate the saturated liquid density.
         * @return The saturated liquid volume [m3/mol]
         * @warning Using the function call operator on a default constructed SLVRackett object will return zero.
         */
        double operator()(double temperature) const;

        // ===== Static Functions ===== //

        /**
         * @brief Factory function, for creating an SLVRackett object from the four Rackett coefficients. This is
         * identical to calling the constructor.
         * @param coeffA Rackett coefficient A.
         * @param coeffB Rackett coefficient B.
         * @param coeffC Rackett coefficient C.
         * @param coeffD Rackett coefficient D.
         * @return An SLVRackett object created from four Rackett coefficients.
         * @note The coefficients must be based on an temperature input in Kelvin, and a molar volume return value
         * in m3/mol.
         */
        static Rackett createFromCoefficients(double coeffA, double coeffB, double coeffC, double coeffD);

        /**
         * @brief Factory function, for creating an SLVRackett object from critical properties.
         * @details The saturated liquid volume can be estimated from critical properties:
         * \f[ V_{s} = {\frac{R \cdot T_{c}}{P_{c}}} \cdot {Z_{c}^{\left [1+(1-\frac{T}{T_{c}})^{\frac{2}{7}}\right]}} \f]
         * By using this static factory function and providing critical temperature [K], critical pressure [Pa] and
         * critical compressibility [-], the Rackett coefficients can be estimated.
         * The compressibility factor can be replaced by an adjustable parameter, \f$ Z_{RA} \f$. If available, it will
         * provide liquid volume estimates that are more accurate away from the critical point. The molar volume at the
         * critical point, however, will only be correct if the critical compressibility is used.
         * @param criticalTemperature The pure component critical temperature [K]
         * @param criticalPressure The pure component critical pressure [Pa]
         * @param criticalCompressibility The pure component critical compressibility [-] or the Rackett compressibility, if available.
         * @return An SLVRackett object created from Rackett coefficients estimated from critical properties.
         */
        static Rackett createFromCriticalProperties(double criticalTemperature, double criticalPressure, double criticalCompressibility);

        /**
         * @brief Factory function for creating an SLVRackett object using the Yamada-Gunn relation.
         * @details Yamada and Gunn in 1973 proposed that \f$ Z_{c} \f$ can be correlated with the acentric factor:
         * \f[ V_{s} = V_{c} \cdot (0.29056 - 0.08775 \cdot \omega)^{(1-\frac{T}{T_{c}})^{\frac{2}{7}}} \f]
         * This expression can be rearranged to the Rackett form used in this class:
         * \f[ V_{s} = {\frac{R \cdot T_{c}}{P_{c}}} \cdot {(0.29056 - 0.08775 \cdot \omega)^{\left [1+
         * (1-\frac{T}{T_{c}})^{\frac{2}{7}}\right]}} \f]
         * From this expression, the four Rackett coefficients are easily calculated. Using this static factory function
         * by providing critical temperature [K], critical volume [m3/mol] and acentric factor [-], the Rackett
         * coefficients can be estimated.
         * @param criticalTemperature The pure component critical temperature [K]
         * @param criticalPressure The pure component critical volume [m3/mol]
         * @param acentricFactor The pure component acentric factor [-]
         * @return An SLVRackett object created from Rackett coefficients estimated using the Yamada-Gunn relation.
         */
        static Rackett createFromAcentricFactor(double criticalTemperature, double criticalPressure, double acentricFactor);

        /**
         * @brief Factory function for creating an SLVRackett object using a known reference point and the Yamada-Gunn relation.
         * @details Using a reference point, with a know temperature and liquid density, the four Rackett coefficients
         * can be estimated using the following equation:
         * \f[ V_{s} = V_{c} \cdot (0.29056 - 0.08775)^{\left(1 - \frac{T}{T_{c}}\right)^{(2/7)}-\left(1 -
         * \frac{T^{ref}}{T_{c}}\right)^{(2/7)}} \f] This expression can be rearranged to the Rackett form used in this class: \f[ V_{s} =
         * {\frac {V_{s}^R}{(0.29056 - 0.08775 \cdot \omega)^{1+\left(1 - \frac{T^{ref}}{T_{c}}\right)^{ (2/7)}}}} \cdot {(0.29056 - 0.08775
         * \cdot \omega)^{\left [1+(1-\frac{T}{T_{c}})^{\frac{2}{7}}\right]}} \f] From this expression, the four Rackett coefficients are
         * easily calculated. Using this static factory function by providing critical temperature [K], experimental temperature [K] and
         * volume [m3/mol] and acentric factor [-], the Rackett coefficients can be estimated.
         * @param criticalTemperature The pure component critical temperature [K]
         * @param experimentalTemperature The temperature of the reference point [K]
         * @param experimentalVolume The molar volume of the reference point [m3/mol]
         * @param acentricFactor The pure component acentric factor [-]
         * @return An SLVRackett object created from Rackett coefficients estimated using a known reference point
         * and the Yamada-Gunn relation.
         */
        static Rackett createFromReferencePointA(double criticalTemperature, double experimentalTemperature, double experimentalVolume, double acentricFactor);

        /**
         * @brief Factory function for creating an SLVRackett object using a known reference point and the critical compressibility.
         * @details Using a reference point, with a know temperature and liquid density, the four Rackett coefficients
         * can be estimated using the following equation:
         * \f[ V_{s} = V_{c} \cdot Z_{c}^{\left(1 - \frac{T}{T_{c}}\right)^{(2/7)}-\left(1 - \frac{T^{ref}}{T_{c}}\right)^{(2/7)}} \f]
         * This expression can be rearranged to the Rackett form used in this class:
         * \f[ V_{s} = {\frac {V_{s}^R}{Z_{c}^{1+\left(1 - \frac{T^{ref}}{T_{c}}\right)^{(2/7)}}}} \cdot {Z_{c}^{\left [1+
         * (1-\frac{T}{T_{c}})^{\frac{2}{7}}\right]}} \f]
         * From this expression, the four Rackett coefficients are easily calculated. Using this static factory function
         * by providing critical temperature [K], experimental temperature [K] and volume [m3/mol] and critical compressibility [-], the
         * Rackett coefficients can be estimated.
         * @param criticalTemperature The pure component critical temperature [K]
         * @param experimentalTemperature The temperature of the reference point [K]
         * @param experimentalVolume The molar volume of the reference point [m3/mol]
         * @param criticalCompressibility The pure component critical compressibility [-]
         * @return An SLVRackett object created from Rackett coefficients estimated using a known reference point
         * and the critical compressibility.
         */
        static Rackett createFromReferencePointB(double criticalTemperature, double experimentalTemperature, double experimentalVolume, double criticalCompressibility);
    };

}    // namespace PCProps::LiquidVolume

#endif    // PCPROPS_RACKETT_HPP
