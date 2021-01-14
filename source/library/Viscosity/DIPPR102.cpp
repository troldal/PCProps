//
// Created by Kenneth Balslev on 11/01/2021.
//

#include <cmath>

#include "DIPPR102.hpp"

namespace PCProps::Viscosity {
    DIPPR102::DIPPR102() = default;

    DIPPR102::DIPPR102(double coeffA, double coeffB, double coeffC, double coeffD)
        : m_A(coeffA),
          m_B(coeffB),
          m_C(coeffC),
          m_D(coeffD)
    {}

    DIPPR102::DIPPR102(const DIPPR102& other) = default;

    DIPPR102::DIPPR102(DIPPR102&& other) noexcept =default;

    DIPPR102::~DIPPR102() = default;

    DIPPR102& DIPPR102::operator=(const DIPPR102& other) = default;

    DIPPR102& DIPPR102::operator=(DIPPR102&& other) noexcept = default;

    double DIPPR102::operator()(double temperature) const
    {
        using std::pow;
        return m_A * pow(temperature, m_B) / (1 + m_C/temperature + m_D/pow(temperature, 2));
    }

} // namespace PCProps::Viscosity