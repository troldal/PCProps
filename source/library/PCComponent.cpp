/*

8888888b.   .d8888b.  8888888b.
888   Y88b d88P  Y88b 888   Y88b
888    888 888    888 888    888
888   d88P 888        888   d88P 888d888 .d88b.  88888b.  .d8888b
8888888P"  888        8888888P"  888P"  d88""88b 888 "88b 88K
888        888    888 888        888    888  888 888  888 "Y8888b.
888        Y88b  d88P 888        888    Y88..88P 888 d88P      X88
888         "Y8888P"  888        888     "Y88P"  88888P"   88888P'
                                                 888
                                                 888
                                                 888

Copyright (c) 2020 Kenneth Troldal Balslev

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#include <algorithm>
#include <cmath>

#include "PCComponent.hpp"

#include <library/PCPropsException.hpp>

namespace
{
    enum PhaseType { Liquid, Vapor, Dense, Undefined };

    PhaseType determinePhaseType(double pressure, double criticalPressure, double saturationPressure, double temperature, double criticalTemperature)
    {
        PhaseType phase = PhaseType::Undefined;

        if (temperature > criticalTemperature && pressure > criticalPressure)
            phase = PhaseType::Dense;
        else if ((temperature > criticalTemperature && pressure <= criticalPressure) || (temperature <= criticalTemperature && pressure <= saturationPressure))
            phase = PhaseType::Vapor;
        else if (temperature <= criticalTemperature && pressure > saturationPressure)
            phase = PhaseType::Liquid;

        return phase;
    }

    double computeCompressedLiquidViscosity(const PCProps::PCComponentData& data, double temperature, double pressure) {

        using std::pow;
        using std::max;

        // ===== Calculate Tr and delta Pr. If the pressure is lower than the vapor pressure, the pressure is assumed
        // ===== equal to the vapor pressure. If the pressure is lower than the vapor pressure, the fluid would be in
        // ===== the gas phase rather than the liquid phase. However, at saturation conditions, calculations may show
        // ===== that the liquid pressure is less than the vapor pressure, for numeric reasons.
        double tr = temperature / data.criticalTemperature.value();
        double dpr = (max(0.0, pressure - data.equationOfState.saturationPressure(temperature))) / data.criticalPressure.value();

        double A = 0.9991 - (4.674E-4 / (1.0523 * pow(tr, -0.03877) - 1.0513));
        double D = (0.3257 / pow(1.0039 - pow(tr, 2.573), 0.2906)) - 0.2086;
        double C = -0.07921 +
            2.1616 * tr -
            13.4040 * pow(tr, 2) +
            44.1706 * pow(tr, 3) -
            84.8291 * pow(tr, 4) +
            96.1209 * pow(tr, 5) -
            59.8127 * pow(tr, 6) +
            15.6719 * pow(tr, 7);

        return (1.0 + D * pow(dpr/2.118, A)) / (1.0 + C * data.acentricFactor.value() * dpr) * data.saturatedLiquidViscosityCorrelation(temperature);
    }

    double computeCompressedVaporViscosity(const PCProps::PCComponentData& data, double temperature, double pressure) {

        using std::pow;
        using std::abs;

        double mu_r = 52.46 * pow(data.dipoleMoment.value(), 2) * (data.criticalPressure.value() / 1E5) * pow(data.criticalTemperature.value(), -2);
        double tr = temperature/data.criticalTemperature.value();
        double pr = pressure/data.criticalPressure.value();

        double Fp_ig = [&]()
        {
            if (mu_r >= 0.0 && mu_r <= 0.022) return 1.0;
            if (mu_r > 0.22 && mu_r <= 0.075) return 1.0 + 30.55 * pow(0.292 - data.criticalCompressibility.value(), 1.72);
            if (mu_r > 0.075) return 1.0 + (30.55 * pow(0.292 - data.criticalCompressibility.value(), 1.72) * abs(0.96 + 0.1 * (tr - 0.7)));
            throw PCProps::PCPropsException("Invalid dipole moment value.");
        }();

        if (tr <= 1.0 && pressure <= data.equationOfState.saturationPressure(temperature)) {

            double A = 3.262 + 14.98 * pow(pr, 5.508);
            double B = 1.39 + 5.746 * pr;
            double Z2 = 0.6 + 0.76 * pow(pr, A) + (6.99 * pow(pr, B) - 0.6) * (1.0 - tr);
            double ksi = 0.176 * pow(data.criticalTemperature.value(), 1.0/6) * pow(data.molecularWeight.value(), -0.5) * pow(data.criticalPressure.value() / 1E5, -2.0/3);
            double eta_ig = data.saturatedVaporViscosityCorrelation(temperature);
            double Fp = (1.0 + (Fp_ig - 1.0) * pow(Z2/(ksi*eta_ig/1.0E-7), -3)) / Fp_ig;
            return Z2 * Fp * 1E-7 / ksi;
        }

        if (tr <= 40.0 && pr <= 100.0) {

            double A = 0.001245 / tr * exp(5.1726 * pow(tr, -0.3286));
            double B = A * (1.6553 * tr - 1.2723);
            double C = 0.4489 / tr * exp(3.0578 * pow(tr, -37.7332));
            double D = 1.7368 / tr * exp(2.2310 * pow(tr, -7.6351));
            double E = 1.3088;
            double F = 0.9425 * exp(-0.1853 * pow(tr, 0.4489));
            double Z2 = 1.0 + (A * pow(pr, E)) / (B * pow(pr, F) + pow(1.0 + C * pow(pr, D), -1));
            double eta_ig = data.saturatedVaporViscosityCorrelation(temperature);
            double Fp = (1.0 + (Fp_ig - 1.0) * pow(Z2, -3)) / Fp_ig;
            return eta_ig * Z2 * Fp;
        }

        throw PCProps::PCPropsException("Lucas Viscosity Estimation Error: Invalid temperature/pressure range.");

    }

}    // namespace

namespace PCProps
{
    // ===== Constructor, default
    PCComponent::PCComponent() = default;

    // ===== Constructor, taking a PCComponentData object as an argument
    PCComponent::PCComponent(const PCComponentData& data) : m_data(data)
    {
        if (!m_data.equationOfState) throw PCPropsException("Error: Invalid EOS object!");
        if (!m_data.idealGasCpCorrelation) throw PCPropsException("Error: Invalid Ideal Gas Cp object!");

        m_data.equationOfState.setProperties(m_data.criticalTemperature.value(), m_data.criticalPressure.value(), m_data.acentricFactor.value());

        m_data.equationOfState.setIdealGasCpFunction([&](double temperature) { return m_data.idealGasCpCorrelation.evaluateCp(temperature); });
        m_data.equationOfState.setIdealGasCpIntegralFunction([&](double temperature) { return m_data.idealGasCpCorrelation.integralOfCp(temperature); });
        m_data.equationOfState.setIdealGasCpOverTIntegralFunction([&](double temperature) { return m_data.idealGasCpCorrelation.integralOfCpOverT(temperature); });
    }

    // ===== Copy constructor
    PCComponent::PCComponent(const PCComponent& other) = default;

    // ===== Move constructor
    PCComponent::PCComponent(PCComponent&& other) noexcept = default;

    // ===== Destructor
    PCComponent::~PCComponent() = default;

    // ===== Copy assignment operator
    PCComponent& PCComponent::operator=(const PCComponent& other) = default;

    // ===== Move assignment operator
    PCComponent& PCComponent::operator=(PCComponent&& other) noexcept = default;

    // ===== Accessor to the PCComponentData member
    PCComponentData& PCComponent::data()
    {
        return m_data;
    }

    PCPhases PCComponent::flash(Utilities::Pressure pressure, Utilities::Temperature temperature) const
    {
        using std::get;
        auto results = PCPhase(m_data.equationOfState.flashPT(pressure.get(), temperature.get())[0]);
        switch (determinePhaseType(
            pressure.get(),
            m_data.criticalPressure.value(),
            m_data.saturatedLiquidVolumeCorrelation(temperature.get()),
            temperature.get(),
            m_data.criticalTemperature.value()))
        {
            case PhaseType::Liquid:
                results.setMolarVolume(m_data.saturatedLiquidVolumeCorrelation(temperature.get()));
                break;

            case PhaseType::Vapor:
            case PhaseType::Dense:
                break;

            default:
                throw PCPropsException("Something went wrong. Invalid phase properties");
        }

        results.setSurfaceTension(0.0);
        results.setThermalConductivity(0.0);
        results.setMolarWeight(m_data.molecularWeight.value());

        PCPhases result = { results.data() };
        computeViscosity(result);

        return { result };
    }

    PCPhases PCComponent::flash(Utilities::Pressure pressure, Utilities::VaporFraction vaporFraction) const
    {
        using std::get;
        PCPhases results;

        for (auto& phase : m_data.equationOfState.flashPx(pressure.get(), vaporFraction.get())) {
            auto result = PCPhase(phase);
            result.setSurfaceTension(0.0);
            result.setThermalConductivity(0.0);
            result.setMolarWeight(m_data.molecularWeight.value());
            results.emplace_back(result.data());
        }

        computeViscosity(results);
        return results;
    }

    PCPhases PCComponent::flash(Utilities::Temperature temperature, Utilities::VaporFraction vaporFraction) const
    {
        using std::get;
        PCPhases results;

        for (auto& phase : m_data.equationOfState.flashTx(temperature.get(), vaporFraction.get())) {
            auto result = PCPhase(phase);
            result.setSurfaceTension(0.0);
            result.setThermalConductivity(0.0);
            result.setMolarWeight(m_data.molecularWeight.value());
            results.emplace_back(result.data());
        }

        computeViscosity(results);
        return results;
    }

    PCPhases PCComponent::flash(Utilities::Pressure pressure, Utilities::Enthalpy enthalpy) const
    {
        using std::get;
        PCPhases results;

        for (auto& phase : m_data.equationOfState.flashPH(pressure.get(), enthalpy.get())) {
            auto result = PCPhase(phase);
            result.setSurfaceTension(0.0);
            result.setThermalConductivity(0.0);
            result.setMolarWeight(m_data.molecularWeight.value());
            results.emplace_back(result.data());
        }

        computeViscosity(results);
        return results;
    }

    PCPhases PCComponent::flash(Utilities::Pressure pressure, Utilities::Entropy entropy) const
    {
        using std::get;
        PCPhases results;

        for (auto& phase : m_data.equationOfState.flashPS(pressure.get(), entropy.get())) {
            auto result = PCPhase(phase);
            result.setSurfaceTension(0.0);
            result.setThermalConductivity(0.0);
            result.setMolarWeight(m_data.molecularWeight.value());
            results.emplace_back(result.data());
        }

        computeViscosity(results);
        return results;
    }

    PCPhases PCComponent::flash(Utilities::Temperature temperature, Utilities::Volume volume) const
    {
        return PCProps::PCPhases();
    }

    double PCComponent::saturationPressure(double temperature) const
    {
        return m_data.equationOfState.saturationPressure(temperature);
    }

    double PCComponent::saturationTemperature(double pressure) const
    {
        return m_data.equationOfState.saturationTemperature(pressure);
    }

    // ===== Get the component name
    const std::string& PCComponent::name() const
    {
        return m_data.name;
    }

    // ===== Get the component formula
    const std::string& PCComponent::formula() const
    {
        return m_data.formula;
    }

    // ===== Get the component CAS registration number
    const std::string& PCComponent::casrn() const
    {
        return m_data.casrn;
    }

    // ===== Get the component SMILES string
    const std::string& PCComponent::smiles() const
    {
        return m_data.smiles;
    }

    void PCComponent::computeViscosity(PCPhases& phases) const {

        if (phases.size() >= 2) {
            auto [min, max] = std::minmax_element(phases.begin(), phases.end(), [&](const auto& a, const auto&b) { return a[PCCompressibility] < b[PCCompressibility]; });
            (*min)[PCViscosity] = computeCompressedLiquidViscosity(m_data, (*min)[PCTemperature], (*min)[PCPressure]);
            (*max)[PCViscosity] = computeCompressedVaporViscosity(m_data, (*max)[PCTemperature], (*max)[PCPressure]);
        }

        else {
            auto psat = saturationPressure(phases[0][PCTemperature]);

            if (phases[0][PCPressure] <= psat)
                phases[0][PCViscosity] = computeCompressedVaporViscosity(m_data, phases[0][PCTemperature], phases[0][PCPressure]);
            else
                phases[0][PCViscosity] = computeCompressedLiquidViscosity(m_data, phases[0][PCTemperature], phases[0][PCPressure]);
        }
    }

} // namespace PCProps