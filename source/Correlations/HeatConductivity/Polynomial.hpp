//
// Created by Kenneth Balslev on 04/01/2021.
//

#ifndef PCPROPS_POLYNOMIAL_HEATCOND_HPP
#define PCPROPS_POLYNOMIAL_HEATCOND_HPP

#include <cmath>

namespace PCProps::HeatConductivity
{
    class Polynomial
    {
        double m_A {};
        double m_B {};
        double m_C {};
        double m_D {};
        double m_E {};
        double m_F {};
        double m_G {};

    public:
        struct CreateFromDIPPR
        {
            double A {};
            double B {};
            double C {};
            double D {};
            double E {};
        };

        Polynomial() = default;

        Polynomial(const CreateFromDIPPR& c) : m_A(c.A), m_B(c.B), m_C(c.C), m_D(c.D), m_E(c.E), m_F(0.0), m_G(0.0) {}

        Polynomial(const Polynomial& other) = default;

        Polynomial(Polynomial&&) noexcept = default;

        ~Polynomial() = default;

        Polynomial& operator=(const Polynomial& other) = default;

        Polynomial& operator=(Polynomial&& other) noexcept = default;

        double operator()(double temperature) const
        {
            using std::pow;
            return m_A + m_B * temperature + m_C * pow(temperature, 2) + m_D * pow(temperature, 3) + m_E * pow(temperature, 4) + m_F * pow(temperature, 5) + m_G * pow(temperature, 6);
        }
    };
}    // namespace PCProps::HeatCapacity

#endif    // PCPROPS_POLYNOMIAL_HEATCOND_HPP
