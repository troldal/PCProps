//
// Created by Kenneth Balslev on 04/01/2021.
//

#include <cmath>

#include "Polynomial.hpp"

namespace PCProps::HeatCapacity
{
    Polynomial::Polynomial() = default;

    Polynomial::Polynomial(const Polynomial::CreateFromDIPPR& c) : m_A(c.A / 1000), m_B(c.B / 1000), m_C(c.C / 1000), m_D(c.D / 1000), m_E(c.E / 1000) {}

    Polynomial::Polynomial(const Polynomial::CreateFromYaws& c) : m_A(c.A), m_B(c.B), m_C(c.C), m_D(c.D) {}

    Polynomial::Polynomial(const Polynomial& other) = default;

    Polynomial::Polynomial(Polynomial&&) noexcept = default;

    Polynomial::~Polynomial() = default;

    Polynomial& Polynomial::operator=(const Polynomial& other) = default;

    Polynomial& Polynomial::operator=(Polynomial&& other) noexcept = default;

    double Polynomial::operator()(double temperature) const
    {
        return evaluateCp(temperature);
    }

    double Polynomial::evaluateCp(double temperature) const
    {
        using std::pow;
        return m_A + m_B * temperature + m_C * pow(temperature, 2) + m_D * pow(temperature, 3) + m_E * pow(temperature, 4) + m_F * pow(temperature, 5) + m_G * pow(temperature, 6);
    }

    double Polynomial::derivativeOfCp(double temperature) const
    {
        using std::pow;
        return m_B + 2 * m_C * temperature + 3 * m_D * pow(temperature, 2) + 4 * m_E * pow(temperature, 3) + 5 * m_F * pow(temperature, 4) + 6 * m_G * pow(temperature, 5);
    }

    double Polynomial::integralOfCp(double temperature) const
    {
        using std::pow;
        return m_A * temperature + m_B * pow(temperature, 2) / 2 + m_C * pow(temperature, 3) / 3 + m_D * pow(temperature, 4) / 4 + m_E * pow(temperature, 5) / 5 +
               m_F * pow(temperature, 6) / 6 + m_G * pow(temperature, 7) / 7;
    }

    double Polynomial::integralOfCpOverT(double temperature) const
    {
        using std::log;
        using std::pow;
        return m_A * log(temperature) + m_B * temperature + m_C * pow(temperature, 2) / 2 + m_D * pow(temperature, 3) / 3 + m_E * pow(temperature, 4) / 4 +
               m_F * pow(temperature, 5) / 5 + m_G * pow(temperature, 6) / 6;
    }

}    // namespace PCProps::HeatCapacity