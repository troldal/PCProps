//
// Created by Kenneth Balslev on 11/01/2021.
//

#include <cmath>
#include <string>

#include "Lucas.hpp"
#include <library/PCPropsException.hpp>

namespace PCProps::Viscosity {
    Lucas::Lucas() = default;

    Lucas::Lucas(double critTemperature,
                 double critPressure,
                 double critCompressibility,
                 double molarWeight,
                 double dipoleMoment)
        : m_criticalTemperature(critTemperature),
          m_criticalPressure(critPressure),
          m_criticalCompressibility(critCompressibility),
          m_molarWeight(molarWeight),
          m_dipoleMoment(dipoleMoment)
    {}

    Lucas::Lucas(const Lucas& other) = default;

    Lucas::Lucas(Lucas&& other) noexcept = default;

    Lucas::~Lucas() = default;

    Lucas& Lucas::operator=(const Lucas& other) = default;

    Lucas& Lucas::operator=(Lucas&& other) noexcept = default;

    double Lucas::operator()(double temperature) const
    {
        using std::pow;
        using std::abs;

        double mu_r = 52.46 * pow(m_dipoleMoment, 2) * (m_criticalPressure / 1E5) * pow(m_criticalTemperature, -2);
        double tr = temperature/m_criticalTemperature;

        double Fp_ig = [&]()
        {
               if (mu_r >= 0.0 && mu_r <= 0.022) return 1.0;
               if (mu_r > 0.22 && mu_r <= 0.075) return 1.0 + 30.55 * pow(0.292 - m_criticalCompressibility, 1.72);
               if (mu_r > 0.075) return 1.0 + (30.55 * pow(0.292 - m_criticalCompressibility, 1.72) * abs(0.96 + 0.1 * (tr - 0.7)));
               throw PCProps::PCPropsException("Invalid dipole moment value.");
        }();

        double ksi = 0.176 * pow(m_criticalTemperature, 1.0/6) * pow(m_molarWeight, -0.5) * pow(m_criticalPressure / 1E5, -2.0/3);

        return Fp_ig/ksi * (0.807 * pow(tr, 0.618) - 0.357 * exp(-0.449*tr) + 0.34 * exp(-4.058*tr) + 0.018) * 1E-7;
    }

} // namespace PCProps::Viscosity