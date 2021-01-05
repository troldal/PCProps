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
        return evaluateCp(temperature);
    }

    double PPDSLiquid::evaluateCp(double temperature) const
    {
        auto tau = 1 - (temperature / m_criticalTemperature);

        return m_A / tau + m_B + m_C * tau + m_D * pow(tau, 2) + m_E * pow(tau, 3) + m_F * pow(tau, 4) + m_G * pow(tau, 5);
    }

    double PPDSLiquid::derivativeOfCp(double temperature) const
    {
        auto tau = 1 - (temperature / m_criticalTemperature);
        return m_A / (m_criticalTemperature * pow(tau, 2)) - m_C / m_criticalTemperature - m_D * 2 * tau / m_criticalTemperature - m_E * 3 * pow(tau, 2) / m_criticalTemperature -
               m_F * 4 * pow(tau, 3) / m_criticalTemperature - m_G * 5 * pow(tau, 4) / m_criticalTemperature;
    }

    double PPDSLiquid::integralOfCp(double temperature) const
    {
        auto tau = 1 - (temperature / m_criticalTemperature);
        return m_B * temperature + m_C * temperature - m_C * pow(temperature, 2) / (2 * m_criticalTemperature) - m_A * m_criticalTemperature * log(tau) -
               m_D * m_criticalTemperature * pow(tau, 3) / 3 - m_E * m_criticalTemperature * pow(tau, 4) / 4 - m_F * m_criticalTemperature * pow(tau, 5) / 5 -
               m_G * m_criticalTemperature * pow(tau, 6) / 6;
    }

    double PPDSLiquid::integralOfCpOverT(double temperature) const
    {
        auto tau = 1 - (temperature / m_criticalTemperature);
        return log(temperature) * (m_A + m_B + m_C + m_D + m_E + m_F + m_G) - m_A * log(tau) - m_C * temperature / m_criticalTemperature -
               2 * m_D * temperature / m_criticalTemperature + m_D * pow(temperature, 2) / (2 * pow(m_criticalTemperature, 2)) - 3 * m_E * temperature / m_criticalTemperature -
               m_E * pow(temperature, 3) / (3 * pow(m_criticalTemperature, 3)) + 3 * m_E * pow(temperature, 2) / (2 * pow(m_criticalTemperature, 2)) -
               4 * m_F * temperature / m_criticalTemperature - 4 * m_F * pow(temperature, 3) / (3 * pow(m_criticalTemperature, 3)) +
               m_F * pow(temperature, 4) / (4 * pow(m_criticalTemperature, 4)) + 3 * m_F * pow(temperature, 2) / pow(m_criticalTemperature, 2) -
               5 * m_G * temperature / m_criticalTemperature - 10 * m_G * pow(temperature, 3) / (3 * pow(m_criticalTemperature, 3)) -
               m_G * pow(temperature, 5) / (5 * pow(m_criticalTemperature, 5)) + 5 * m_G * pow(temperature, 4) / (4 * pow(m_criticalTemperature, 4)) +
               5 * m_G * pow(temperature, 2) / pow(m_criticalTemperature, 2);
    }

}    // namespace PCProps::HeatCapacity