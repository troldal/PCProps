//
// Created by Kenneth Balslev on 19/01/2021.
//

#ifndef PCPROPS_COMPRESSEDLIQUID_LUCAS_HPP
#define PCPROPS_COMPRESSEDLIQUID_LUCAS_HPP

#include <cmath>
#include <functional>

namespace PCProps::CompressedLiquidViscosity
{
    /**
     *
     */
    class Lucas
    {
        double m_criticalPressure{};
        double m_criticalTemperature{};
        double m_acentricFactor{};

    public:

        /**
         *
         */
        Lucas() = default;

        /**
         *
         * @param criticalTemperature
         * @param criticalPressure
         * @param acentricFactor
         */
        Lucas(double criticalTemperature,
              double criticalPressure,
              double acentricFactor)
            : m_criticalTemperature(criticalTemperature),
              m_criticalPressure(criticalPressure),
              m_acentricFactor(acentricFactor)
        {}

        /**
         *
         * @tparam PC
         * @param pureComponent
         */
        template<typename PC>
        explicit Lucas(const PC& pureComponent) : Lucas(pureComponent.property("CriticalTemperature"),
                                               pureComponent.property("CriticalPressure"),
                                               pureComponent.property("AcentricFactor"))
        {}

        /**
         *
         * @param other
         */
        Lucas(const Lucas& other) = default;

        /**
         *
         * @param other
         */
        Lucas(Lucas&& other) noexcept = default;

        /**
         *
         */
        ~Lucas() = default;

        /**
         *
         * @param other
         * @return
         */
        Lucas& operator=(const Lucas& other) = default;

        /**
         *
         * @param other
         * @return
         */
        Lucas& operator=(Lucas&& other) noexcept = default;

        /**
         *
         * @param params
         * @return
         */
        double operator()(std::vector<double> params) const {
            switch (params.size()) {
                case 4:
                    return operator()(params[0], params[1], params[2], params[3]);

                default:
                    return 0.0;
            }
        }

        /**
         *
         * @param temperature
         * @param pressure
         * @param satPressure
         * @param satViscosity
         * @return
         */
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
    };
}    // namespace PCProps::CompressedLiquidViscosity

#endif    // PCPROPS_COMPRESSEDLIQUID_LUCAS_HPP
