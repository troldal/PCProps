//
// Created by Kenneth Balslev on 04/01/2021.
//

#include <cmath>

#include "PPDSLiquid.hpp"

using std::log;
using std::pow;

namespace PCProps::HeatCapacity
{
    PPDSLiquid::PPDSLiquid() = default;

    PPDSLiquid::PPDSLiquid(const PPDSLiquid::CreateFromPPDS& c)
        : m_A(c.A * 8.31446261815324),
          m_B(c.B * 8.31446261815324),
          m_C(c.C * 8.31446261815324),
          m_D(c.D * 8.31446261815324),
          m_E(c.E * 8.31446261815324),
          m_F(c.F * 8.31446261815324),
          m_criticalTemperature(c.criticalTemperature)
    {}

    PPDSLiquid::PPDSLiquid(const PPDSLiquid::CreateFromDIPPR& c)
        : m_A(pow(c.A, 2) / 1000),
          m_B(c.B / 1000),
          m_C(-2.0 * c.A * c.C / 1000),
          m_D(-c.A * c.D / 1000),
          m_E(-pow(c.C, 2) / 3000),
          m_F(-c.C * c.D / 2000),
          m_G(-pow(c.D, 2) / 5000),
          m_criticalTemperature(c.criticalTemperature)
    {}

    PPDSLiquid::PPDSLiquid(const PPDSLiquid& other) = default;

    PPDSLiquid::PPDSLiquid(PPDSLiquid&& other) noexcept = default;

    PPDSLiquid::~PPDSLiquid() = default;

    PPDSLiquid& PPDSLiquid::operator=(const PPDSLiquid& other) = default;

    PPDSLiquid& PPDSLiquid::operator=(PPDSLiquid&& other) noexcept = default;

    double PPDSLiquid::operator()(double temperature) const
    {
        auto tau = 1 - (temperature / m_criticalTemperature);
        return m_A / tau + m_B + m_C * tau + m_D * pow(tau, 2) + m_E * pow(tau, 3) + m_F * pow(tau, 4) + m_G * pow(tau, 5);
    }
}    // namespace PCProps::HeatCapacity