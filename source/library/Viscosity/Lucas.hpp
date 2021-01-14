//
// Created by Kenneth Balslev on 11/01/2021.
//

#ifndef PCPROPS_LUCAS_HPP
#define PCPROPS_LUCAS_HPP

#include <functional>

namespace PCProps::Viscosity
{
    class Lucas
    {
        double m_criticalTemperature {};
        double m_criticalPressure {};
        double m_criticalCompressibility {};
        double m_molarWeight {};
        double m_dipoleMoment {};

    public:

        /**
         * @brief Constructor, default. Sets all coefficients to zero.
         * @warning Calling the function call operator on a default constructed SLVRackett object will result in
         * division by zero and return NaN.
         */
        Lucas();

        /**
         * @brief
         * @param critTemperature
         * @param critPressure
         * @param critCompressibility
         * @param molarWeight
         * @param dipoleMoment
         * @param vaporPressureFunction
         * @param idealGasViscosityFunction
         */
        Lucas(
            double                               critTemperature,
            double                               critPressure,
            double                               critCompressibility,
            double                               molarWeight,
            double                               dipoleMoment);

        /**
         * @brief Copy constructor.
         */
        Lucas(const Lucas& other);

        /**
         * @brief Move constructor.
         */
        Lucas(Lucas&& other) noexcept;

        /**
         * @brief Destructor.
         */
        ~Lucas();

        // ===== Manipulators ===== //

        /**
         * @brief Copy assignment operator.
         */
        Lucas& operator=(const Lucas& other);

        /**
         * @brief Move assignment operator.
         */
        Lucas& operator=(Lucas&& other) noexcept;

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

#endif    // PCPROPS_LUCAS_HPP
