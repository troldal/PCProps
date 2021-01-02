//
// Created by Kenneth Balslev on 31/12/2020.
//

#include <cmath>

#include "IGAlyLee.hpp"

namespace PCProps::HeatCapacity
{
    IGAlyLee::IGAlyLee() = default;

    IGAlyLee::IGAlyLee(double coeffA, double coeffB, double coeffC, double coeffD, double coeffE)
        : m_A(coeffA),
          m_B(coeffB),
          m_C(coeffC),
          m_D(coeffD),
          m_E(coeffE)
    {}

    IGAlyLee::IGAlyLee(const IGAlyLee& other) = default;

    IGAlyLee::IGAlyLee(IGAlyLee&& other) noexcept = default;

    IGAlyLee::~IGAlyLee() = default;

    IGAlyLee& IGAlyLee::operator=(const IGAlyLee& other) = default;

    IGAlyLee& IGAlyLee::operator=(IGAlyLee&& other) noexcept = default;

    double IGAlyLee::operator()(double temperature)
    {
        return evaluate(temperature);
    }

    double IGAlyLee::evaluate(double temperature)
    {
        using std::cosh;
        using std::pow;
        using std::sinh;

        return (m_A + m_B * pow((m_C / temperature) / sinh(m_C / temperature), 2) +
                m_D * pow((m_E / temperature) / cosh(m_E / temperature), 2)) /
               1000;
    }

    double IGAlyLee::differential(double temperature)
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

    double IGAlyLee::integral(double temperature)
    {
        using std::tanh;

        return (m_A * temperature + m_B * m_C / tanh(m_C / temperature) - m_D * m_E * tanh(m_E / temperature)) / 1000;
    }

}    // namespace PCProps::HeatCapacity
