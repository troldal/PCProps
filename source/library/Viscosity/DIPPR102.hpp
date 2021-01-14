//
// Created by Kenneth Balslev on 11/01/2021.
//

#ifndef PCPROPS_DIPPR102_HPP
#define PCPROPS_DIPPR102_HPP

namespace PCProps::Viscosity
{
    class DIPPR102
    {
        double m_A = 0.0;
        double m_B = 0.0;
        double m_C = 0.0;
        double m_D = 0.0;

    public:
        /**
 * @brief Constructor, default. Sets all coefficients to zero.
 * @warning Calling the function call operator on a default constructed SLVRackett object will result in
 * division by zero and return NaN.
 */
        DIPPR102();

        /**
         * @brief Constructor, taking the four Rackett coefficients as arguments.
         * @param coeffA Rackett coefficient A.
         * @param coeffB Rackett coefficient B.
         * @param coeffC Rackett coefficient C.
         * @param coeffD Rackett coefficient D.
         * @note The coefficients must be based on an temperature input in Kelvin, and a density return value
         * in m3/mol.
         */
        DIPPR102(double coeffA, double coeffB, double coeffC = 0.0, double coeffD = 0.0);

        /**
         * @brief Copy constructor.
         */
        DIPPR102(const DIPPR102& other);

        /**
         * @brief Move constructor.
         */
        DIPPR102(DIPPR102&& other) noexcept;

        /**
         * @brief Destructor.
         */
        ~DIPPR102();

        // ===== Manipulators ===== //

        /**
         * @brief Copy assignment operator.
         */
        DIPPR102& operator=(const DIPPR102& other);

        /**
         * @brief Move assignment operator.
         */
        DIPPR102& operator=(DIPPR102&& other) noexcept;

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

#endif    // PCPROPS_DIPPR102_HPP
