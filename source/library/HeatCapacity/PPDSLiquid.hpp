//
// Created by Kenneth Balslev on 04/01/2021.
//

#ifndef PCPROPS_PPDSLIQUID_HPP
#define PCPROPS_PPDSLIQUID_HPP

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

        PPDSLiquid();

        explicit PPDSLiquid(const CreateFromPPDS& coefficients);

        explicit PPDSLiquid(const CreateFromDIPPR& coefficients);

        PPDSLiquid(const PPDSLiquid& other);

        PPDSLiquid(PPDSLiquid&& other) noexcept;

        ~PPDSLiquid();

        PPDSLiquid& operator=(const PPDSLiquid& other);

        PPDSLiquid& operator=(PPDSLiquid&& other) noexcept;

        double operator()(double temperature) const;

        double evaluateCp(double temperature) const;

        double derivativeOfCp(double temperature) const;

        double integralOfCp(double temperature) const;

        double integralOfCpOverT(double temperature) const;
    };
}    // namespace PCProps::HeatCapacity

#endif    // PCPROPS_PPDSLIQUID_HPP
