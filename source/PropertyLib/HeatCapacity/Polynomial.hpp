//
// Created by Kenneth Balslev on 04/01/2021.
//

#ifndef PCPROPS_POLYNOMIAL_HPP
#define PCPROPS_POLYNOMIAL_HPP

#include <cmath>

namespace PCProps::HeatCapacity
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

        struct CreateFromYaws
        {
            double A {};
            double B {};
            double C {};
            double D {};
        };

        Polynomial() = default;

        Polynomial(const CreateFromDIPPR& c) : m_A(c.A / 1000), m_B(c.B / 1000), m_C(c.C / 1000), m_D(c.D / 1000), m_E(c.E / 1000) {}

        Polynomial(const CreateFromYaws& c) : m_A(c.A), m_B(c.B), m_C(c.C), m_D(c.D) {}

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

#endif    // PCPROPS_POLYNOMIAL_HPP
