//
// Created by Kenneth Balslev on 11/01/2021.
//

#ifndef PCPROPS_KIRCHHOFFEXTENDED_HPP
#define PCPROPS_KIRCHHOFFEXTENDED_HPP

namespace PCProps::Viscosity
{
    class KirchhoffExtended
    {
        double m_A {};
        double m_B {};
        double m_C {};
        double m_D {};
        double m_E {};

    public:

        struct CreateFromDIPPR
        {
            double A;
            double B;
            double C;
            double D;
            double E;
        };

        /**
* @brief Constructor, default. Sets all coefficients to zero.
* @warning Calling the function call operator on a default constructed SLVRackett object will result in
* division by zero and return NaN.
*/
        KirchhoffExtended();

        /**
         * @brief Constructor, taking the four Rackett coefficients as arguments.
         * @param coeffA Rackett coefficient A.
         * @param coeffB Rackett coefficient B.
         * @param coeffC Rackett coefficient C.
         * @param coeffD Rackett coefficient D.
         * @note The coefficients must be based on an temperature input in Kelvin, and a density return value
         * in m3/mol.
         */
        KirchhoffExtended(double coeffA, double coeffB, double coeffC = 0.0, double coeffD = 0.0, double coeffE = 0.0);

        /**
         * @brief
         * @param coefficients
         */
        KirchhoffExtended(const CreateFromDIPPR& coefficients);

        /**
         * @brief Copy constructor.
         */
        KirchhoffExtended(const KirchhoffExtended& other);

        /**
         * @brief Move constructor.
         */
        KirchhoffExtended(KirchhoffExtended&& other) noexcept;

        /**
         * @brief Destructor.
         */
        ~KirchhoffExtended();

        // ===== Manipulators ===== //

        /**
         * @brief Copy assignment operator.
         */
        KirchhoffExtended& operator=(const KirchhoffExtended& other);

        /**
         * @brief Move assignment operator.
         */
        KirchhoffExtended& operator=(KirchhoffExtended&& other) noexcept;

        // ===== Accessors ===== //

        /**
         * @brief Function call operator. This function is used to calculate the saturated liquid volume [m3/mol] at a
         * given temperature [K]
         * @param temperature The temperature [K] at which to evaluate the saturated liquid density.
         * @return The saturated liquid volume [m3/mol]
         * @warning Using the function call operator on a default constructed SLVRackett object will return zero.
         */
        double operator()(double temperature) const;


    };
}    // namespace PCProps::Viscosity

#endif    // PCPROPS_KIRCHHOFFEXTENDED_HPP
