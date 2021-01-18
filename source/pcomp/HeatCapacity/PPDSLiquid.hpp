//
// Created by Kenneth Balslev on 04/01/2021.
//

#ifndef PCPROPS_PPDSLIQUID_HPP
#define PCPROPS_PPDSLIQUID_HPP

#include <cmath>

namespace PCProps::HeatCapacity
{
    class PPDSLiquid
    {
        double m_A {};
        double m_B {};
        double m_C {};
        double m_D {};
        double m_E {};
        double m_F {};
        double m_G {};
        double m_criticalTemperature {};

    public:
        struct CreateFromPPDS
        {
            double A {};
            double B {};
            double C {};
            double D {};
            double E {};
            double F {};
            double criticalTemperature {};
            double molecularWeight {};
        };

        struct CreateFromDIPPR
        {
            double A {};
            double B {};
            double C {};
            double D {};
            double criticalTemperature {};
        };

        PPDSLiquid() = default;

        explicit PPDSLiquid(const CreateFromPPDS& c)
            : m_A(c.A * 8.31446261815324),
              m_B(c.B * 8.31446261815324),
              m_C(c.C * 8.31446261815324),
              m_D(c.D * 8.31446261815324),
              m_E(c.E * 8.31446261815324),
              m_F(c.F * 8.31446261815324),
              m_criticalTemperature(c.criticalTemperature)
        {}

        explicit PPDSLiquid(const CreateFromDIPPR& c)
            : m_A(pow(c.A, 2) / 1000),
              m_B(c.B / 1000),
              m_C(-2.0 * c.A * c.C / 1000),
              m_D(-c.A * c.D / 1000),
              m_E(-pow(c.C, 2) / 3000),
              m_F(-c.C * c.D / 2000),
              m_G(-pow(c.D, 2) / 5000),
              m_criticalTemperature(c.criticalTemperature)
        {}

        PPDSLiquid(const PPDSLiquid& other) = default;

        PPDSLiquid(PPDSLiquid&& other) noexcept = default;

        ~PPDSLiquid() = default;

        PPDSLiquid& operator=(const PPDSLiquid& other) = default;

        PPDSLiquid& operator=(PPDSLiquid&& other) noexcept = default;

        double operator()(double temperature) const
        {
            auto tau = 1 - (temperature / m_criticalTemperature);
            return m_A / tau + m_B + m_C * tau + m_D * pow(tau, 2) + m_E * pow(tau, 3) + m_F * pow(tau, 4) + m_G * pow(tau, 5);
        }
    };
}    // namespace PCProps::HeatCapacity

#endif    // PCPROPS_PPDSLIQUID_HPP
