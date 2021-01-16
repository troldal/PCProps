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
        using std::cosh;
        using std::pow;
        using std::sinh;

        return m_A + m_B * pow((m_C / temperature) / sinh(m_C / temperature), 2) + m_D * pow((m_E / temperature) / cosh(m_E / temperature), 2);
    }
}    // namespace PCProps::HeatCapacity
