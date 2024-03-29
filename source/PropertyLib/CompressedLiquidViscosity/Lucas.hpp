//
// Created by Kenneth Balslev on 19/01/2021.
//

#ifndef PCPROPS_COMPRESSEDLIQUID_LUCAS_HPP
#define PCPROPS_COMPRESSEDLIQUID_LUCAS_HPP

#include <cmath>
#include <functional>

namespace PCProps::CompressedLiquidViscosity
{
    class Lucas
    {

        double m_criticalPressure{};
        double m_criticalTemperature{};
        double m_acentricFactor{};

        std::function<double(double)> m_satLiquidViscosity {};
        std::function<double(double)> m_vaporPressure {};

    public:

        Lucas() = default;

        Lucas(double criticalTemperature,
              double criticalPressure,
              double acentricFactor,
              const std::function<double(double)>& satLiquidViscosity = {},
              const std::function<double(double)>& vaporPressureFunction = {})
            : m_criticalTemperature(criticalTemperature),
              m_criticalPressure(criticalPressure),
              m_acentricFactor(acentricFactor),
              m_satLiquidViscosity{satLiquidViscosity},
              m_vaporPressure{vaporPressureFunction}
        {}

        Lucas(const Lucas& other) = default;

        Lucas(Lucas&& other) noexcept = default;

        ~Lucas() = default;

        Lucas& operator=(const Lucas& other) = default;

        Lucas& operator=(Lucas&& other) noexcept = default;

        double operator()(double temperature, double pressure, double satPressure, double satViscosity) const {

            using std::max;
            using std::pow;

            // ===== Calculate Tr and delta Pr. If the pressure is lower than the vapor pressure, the pressure is assumed
            // ===== equal to the vapor pressure. If the pressure is lower than the vapor pressure, the fluid would be in
            // ===== the gas phase rather than the liquid phase. However, at saturation conditions, calculations may show
            // ===== that the liquid pressure is less than the vapor pressure, for numeric reasons.
            double tr  = temperature / m_criticalTemperature;
            double dpr = (max(0.0, pressure - satPressure)) / m_criticalPressure;

            double A = 0.9991 - (4.674E-4 / (1.0523 * pow(tr, -0.03877) - 1.0513));
            double D = (0.3257 / pow(1.0039 - pow(tr, 2.573), 0.2906)) - 0.2086;
            double C =
                       -0.07921 + 2.1616 * tr -
                           13.4040 * pow(tr, 2) +
                           44.1706 * pow(tr, 3) -
                           84.8291 * pow(tr, 4) +
                           96.1209 * pow(tr, 5) -
                           59.8127 * pow(tr, 6) +
                           15.6719 * pow(tr, 7);

            return (1.0 + D * pow(dpr / 2.118, A)) / (1.0 + C * m_acentricFactor * dpr) * satViscosity;
        }


        double operator()(double temperature, double pressure) const
        {

            return operator()(temperature, pressure, m_vaporPressure(temperature), m_satLiquidViscosity(temperature));
//            using std::max;
//            using std::pow;
//
//            // ===== Calculate Tr and delta Pr. If the pressure is lower than the vapor pressure, the pressure is assumed
//            // ===== equal to the vapor pressure. If the pressure is lower than the vapor pressure, the fluid would be in
//            // ===== the gas phase rather than the liquid phase. However, at saturation conditions, calculations may show
//            // ===== that the liquid pressure is less than the vapor pressure, for numeric reasons.
//            double tr  = temperature / m_criticalTemperature;
//            double dpr = (max(0.0, pressure - m_vaporPressure(temperature))) / m_criticalPressure;
//
//            double A = 0.9991 - (4.674E-4 / (1.0523 * pow(tr, -0.03877) - 1.0513));
//            double D = (0.3257 / pow(1.0039 - pow(tr, 2.573), 0.2906)) - 0.2086;
//            double C =
//                -0.07921 + 2.1616 * tr -
//                13.4040 * pow(tr, 2) +
//                44.1706 * pow(tr, 3) -
//                84.8291 * pow(tr, 4) +
//                96.1209 * pow(tr, 5) -
//                59.8127 * pow(tr, 6) +
//                15.6719 * pow(tr, 7);
//
//            return (1.0 + D * pow(dpr / 2.118, A)) / (1.0 + C * m_acentricFactor * dpr) * m_satLiquidViscosity(temperature);
        }
    };
}    // namespace PCProps::CompressedLiquidViscosity

#endif    // PCPROPS_COMPRESSEDLIQUID_LUCAS_HPP
