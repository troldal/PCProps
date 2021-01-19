//
// Created by Kenneth Balslev on 19/01/2021.
//

#ifndef PCPROPS_COMPRESSEDVAPOR_LUCAS_HPP
#define PCPROPS_COMPRESSEDVAPOR_LUCAS_HPP

#include <stdexcept>

namespace PCProps::CompressedVaporViscosity
{
    class Lucas
    {

        double m_criticalTemperature{};
        double m_criticalPressure{};
        double m_criticalCompressibility{};
        double m_molarWeight {};
        double m_dipoleMoment{};

        std::function<double(double)> m_satVaporViscosity {};
        std::function<double(double)> m_vaporPressure {};

    public:

        Lucas() = default;

        Lucas(double criticalTemperature,
              double criticalPressure,
              double criticalCompressibility,
              double molarWeight,
              double dipoleMoment,
              const std::function<double(double)>& satVaporViscosity,
              const std::function<double(double)>& vaporPressureFunction)
            : m_criticalTemperature(criticalTemperature),
              m_criticalPressure(criticalPressure),
              m_criticalCompressibility(criticalCompressibility),
              m_molarWeight(molarWeight),
              m_dipoleMoment(dipoleMoment),
              m_satVaporViscosity{satVaporViscosity},
              m_vaporPressure{vaporPressureFunction}
        {}

        Lucas(const Lucas& other) = default;

        Lucas(Lucas&& other) noexcept = default;

        ~Lucas() = default;

        Lucas& operator=(const Lucas& other) = default;

        Lucas& operator=(Lucas&& other) noexcept = default;

        double operator()(double temperature, double pressure) const
        {
            using std::abs;
            using std::pow;

            double mu_r = 52.46 * pow(m_dipoleMoment, 2) * (m_criticalPressure / 1E5) * pow(m_criticalTemperature, -2);
            double tr   = temperature / m_criticalTemperature;
            double pr   = pressure / m_criticalPressure;

            double Fp_ig = [&]() {
                if (mu_r >= 0.0 && mu_r <= 0.022) return 1.0;
                if (mu_r > 0.22 && mu_r <= 0.075) return 1.0 + 30.55 * pow(0.292 - m_criticalCompressibility, 1.72);
                if (mu_r > 0.075) return 1.0 + (30.55 * pow(0.292 - m_criticalCompressibility, 1.72) * abs(0.96 + 0.1 * (tr - 0.7)));
                throw std::invalid_argument("Invalid dipole moment value.");
            }();

            if (tr <= 1.0 && pressure <= m_vaporPressure(temperature)) {
                double A      = 3.262 + 14.98 * pow(pr, 5.508);
                double B      = 1.39 + 5.746 * pr;
                double Z2     = 0.6 + 0.76 * pow(pr, A) + (6.99 * pow(pr, B) - 0.6) * (1.0 - tr);
                double ksi    = 0.176 * pow(m_criticalTemperature, 1.0 / 6) * pow(m_molarWeight, -0.5) * pow(m_criticalPressure / 1E5, -2.0 / 3);
                double eta_ig = m_satVaporViscosity(temperature);
                double Fp     = (1.0 + (Fp_ig - 1.0) * pow(Z2 / (ksi * eta_ig / 1.0E-7), -3)) / Fp_ig;
                return Z2 * Fp * 1E-7 / ksi;
            }

            if (tr <= 40.0 && pr <= 100.0) {
                double A      = 0.001245 / tr * exp(5.1726 * pow(tr, -0.3286));
                double B      = A * (1.6553 * tr - 1.2723);
                double C      = 0.4489 / tr * exp(3.0578 * pow(tr, -37.7332));
                double D      = 1.7368 / tr * exp(2.2310 * pow(tr, -7.6351));
                double E      = 1.3088;
                double F      = 0.9425 * exp(-0.1853 * pow(tr, 0.4489));
                double Z2     = 1.0 + (A * pow(pr, E)) / (B * pow(pr, F) + pow(1.0 + C * pow(pr, D), -1));
                double eta_ig = m_satVaporViscosity(temperature);
                double Fp     = (1.0 + (Fp_ig - 1.0) * pow(Z2, -3)) / Fp_ig;
                return eta_ig * Z2 * Fp;
            }

            throw std::invalid_argument("Lucas Viscosity Estimation Error: Invalid temperature/pressure range.");
        }
    };
}    // namespace PCProps::CompressedVaporViscosity

#endif    // PCPROPS_COMPRESSEDVAPOR_LUCAS_HPP
