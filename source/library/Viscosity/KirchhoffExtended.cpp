//
// Created by Kenneth Balslev on 11/01/2021.
//

#include <cmath>

#include "KirchhoffExtended.hpp"

namespace PCProps::Viscosity {
    KirchhoffExtended::KirchhoffExtended() = default;

    KirchhoffExtended::KirchhoffExtended(double coeffA, double coeffB, double coeffC, double coeffD, double coeffE)
        : m_A(coeffA),
          m_B(coeffB),
          m_C(coeffC),
          m_D(coeffD),
          m_E(coeffE)
    {}

    KirchhoffExtended::KirchhoffExtended(const KirchhoffExtended::CreateFromDIPPR& c)
        : KirchhoffExtended(c.A, c.B, c.C, c.D, c.E)
    {}

    KirchhoffExtended::KirchhoffExtended(const KirchhoffExtended& other) = default;

    KirchhoffExtended::KirchhoffExtended(KirchhoffExtended&& other) noexcept = default;

    KirchhoffExtended::~KirchhoffExtended() = default;

    KirchhoffExtended& KirchhoffExtended::operator=(const KirchhoffExtended& other) = default;

    KirchhoffExtended& KirchhoffExtended::operator=(KirchhoffExtended&& other) noexcept = default;

    double KirchhoffExtended::operator()(double temperature) const
    {
        using std::exp;
        using std::pow;
        using std::log;
        return exp(m_A + m_B/temperature + m_C * log(temperature) + m_D * pow(temperature, m_E));
    }

} // namespace PCProps::Viscosity