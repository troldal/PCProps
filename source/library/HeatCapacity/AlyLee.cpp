//
// Created by Kenneth Balslev on 31/12/2020.
//

#include <cmath>

#include "AlyLee.hpp"

namespace PCProps::HeatCapacity
{
    AlyLee::AlyLee() = default;

    AlyLee::AlyLee(double coeffA, double coeffB, double coeffC, double coeffD, double coeffE) : m_A(coeffA), m_B(coeffB), m_C(coeffC), m_D(coeffD), m_E(coeffE) {}

    AlyLee::AlyLee(const AlyLee::CreateFromDIPPR& c) : AlyLee(c.A / 1000, c.B / 1000, c.C, c.D / 1000, c.E) {}

    AlyLee::AlyLee(const AlyLee& other) = default;

    AlyLee::AlyLee(AlyLee&& other) noexcept = default;

    AlyLee::~AlyLee() = default;

    AlyLee& AlyLee::operator=(const AlyLee& other) = default;

    AlyLee& AlyLee::operator=(AlyLee&& other) noexcept = default;

    double AlyLee::operator()(double temperature)
    {
        return evaluateCp(temperature);
    }

    double AlyLee::evaluateCp(double temperature) const
    {
        using std::cosh;
        using std::pow;
        using std::sinh;

        return m_A + m_B * pow((m_C / temperature) / sinh(m_C / temperature), 2) + m_D * pow((m_E / temperature) / cosh(m_E / temperature), 2);
    }

    double AlyLee::derivativeOfCp(double temperature) const
    {
        using std::cosh;
        using std::pow;
        using std::sinh;

        return (2 * m_B * pow(m_C, 3) * cosh(m_C / temperature) / pow(temperature, 4) / pow(sinh(m_C / temperature), 3) -
                2 * m_B * pow(m_C, 2) / pow(temperature, 3) / pow(sinh(m_C / temperature), 2) +
                2 * m_D * pow(m_E, 3) * sinh(m_E / temperature) / pow(temperature, 4) / pow(cosh(m_E / temperature), 3) -
                2 * m_D * pow(m_E, 2) / pow(temperature, 3) / pow(cosh(m_E / temperature), 2)) /
               1000;
    }

    double AlyLee::integralOfCp(double temperature) const
    {
        using std::tanh;

        return (m_A * temperature + m_B * m_C / tanh(m_C / temperature) - m_D * m_E * tanh(m_E / temperature)) / 1000;
    }

    double AlyLee::integralOfCpOverT(double temperature) const
    {
        using std::cosh;
        using std::pow;
        using std::sinh;
        using std::tanh;

        return m_A * log(temperature) + m_B * m_C / (temperature * tanh(m_C / temperature)) - m_B * log(sinh(m_C / temperature)) -
               m_D * m_E / temperature * tanh(m_E / temperature) + m_D * log(cosh(m_E / temperature));
    }

}    // namespace PCProps::HeatCapacity
