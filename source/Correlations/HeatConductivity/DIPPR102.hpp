//
// Created by Kenneth Balslev on 11/01/2021.
//

#ifndef PCPROPS_DIPPR102_HEATCOND_HPP
#define PCPROPS_DIPPR102_HEATCOND_HPP

#include <cmath>

namespace PCProps::HeatConductivity
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
        DIPPR102() = default;

        /**
         * @brief Constructor, taking the four coefficients as arguments.
         * @param coeffA Rackett coefficient A.
         * @param coeffB Rackett coefficient B.
         * @param coeffC Rackett coefficient C.
         * @param coeffD Rackett coefficient D.
         * @note The coefficients must be based on an temperature input in Kelvin, and a heat conductivity return value
         * in W/m-K.
         */
        DIPPR102(double coeffA, double coeffB, double coeffC = 0.0, double coeffD = 0.0)
            : m_A(coeffA),
              m_B(coeffB),
              m_C(coeffC),
              m_D(coeffD)
        {}

        /**
         * @brief Copy constructor.
         */
        DIPPR102(const DIPPR102& other) = default;

        /**
         * @brief Move constructor.
         */
        DIPPR102(DIPPR102&& other) noexcept = default;

        /**
         * @brief Destructor.
         */
        ~DIPPR102() = default;

        // ===== Manipulators ===== //

        /**
         * @brief Copy assignment operator.
         */
        DIPPR102& operator=(const DIPPR102& other) = default;

        /**
         * @brief Move assignment operator.
         */
        DIPPR102& operator=(DIPPR102&& other) noexcept = default;

        // ===== Accessors ===== //

        /**
         * @brief Function call operator. This function is used to calculate the saturated liquid volume [m3/mol] at a
         * given temperature [K]
         * @param temperature The temperature [K] at which to evaluate the saturated liquid density.
         * @return The saturated liquid volume [m3/mol]
         * @warning Using the function call operator on a default constructed SLVRackett object will return zero.
         */
        double operator()(double temperature) const
        {
            using std::pow;
            return m_A * pow(temperature, m_B) / (1 + m_C/temperature + m_D/pow(temperature, 2));
        }


    };
}    // namespace PCProps::Viscosity

#endif    // PCPROPS_DIPPR102_HEATCOND_HPP
