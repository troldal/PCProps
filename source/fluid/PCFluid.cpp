//
// Created by Kenneth Balslev on 18/01/2021.
//

#include "PCFluid.hpp"

#include <json/json.hpp>

#include <stdexcept>

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

//    double computeCompressedLiquidViscosity(const PCProps::PCComponentData& data, double temperature, double pressure) {
//
//        using std::pow;
//        using std::max;
//
//        // ===== Calculate Tr and delta Pr. If the pressure is lower than the vapor pressure, the pressure is assumed
//        // ===== equal to the vapor pressure. If the pressure is lower than the vapor pressure, the fluid would be in
//        // ===== the gas phase rather than the liquid phase. However, at saturation conditions, calculations may show
//        // ===== that the liquid pressure is less than the vapor pressure, for numeric reasons.
//        double tr = temperature / data.criticalTemperature.value();
//        double dpr = (max(0.0, pressure - data.equationOfState.saturationPressure(temperature))) / data.criticalPressure.value();
//
//        double A = 0.9991 - (4.674E-4 / (1.0523 * pow(tr, -0.03877) - 1.0513));
//        double D = (0.3257 / pow(1.0039 - pow(tr, 2.573), 0.2906)) - 0.2086;
//        double C = -0.07921 +
//            2.1616 * tr -
//            13.4040 * pow(tr, 2) +
//            44.1706 * pow(tr, 3) -
//            84.8291 * pow(tr, 4) +
//            96.1209 * pow(tr, 5) -
//            59.8127 * pow(tr, 6) +
//            15.6719 * pow(tr, 7);
//
//        return (1.0 + D * pow(dpr/2.118, A)) / (1.0 + C * data.acentricFactor.value() * dpr) * data.saturatedLiquidViscosityCorrelation(temperature);
//    }

//    double computeCompressedVaporViscosity(const PCProps::PCComponentData& data, double temperature, double pressure) {
//
//        using std::pow;
//        using std::abs;
//
//        double mu_r = 52.46 * pow(data.dipoleMoment.value(), 2) * (data.criticalPressure.value() / 1E5) * pow(data.criticalTemperature.value(), -2);
//        double tr = temperature/data.criticalTemperature.value();
//        double pr = pressure/data.criticalPressure.value();
//
//        double Fp_ig = [&]()
//        {
//            if (mu_r >= 0.0 && mu_r <= 0.022) return 1.0;
//            if (mu_r > 0.22 && mu_r <= 0.075) return 1.0 + 30.55 * pow(0.292 - data.criticalCompressibility.value(), 1.72);
//            if (mu_r > 0.075) return 1.0 + (30.55 * pow(0.292 - data.criticalCompressibility.value(), 1.72) * abs(0.96 + 0.1 * (tr - 0.7)));
//            throw PCProps::PCPropsException("Invalid dipole moment value.");
//        }();
//
//        if (tr <= 1.0 && pressure <= data.equationOfState.saturationPressure(temperature)) {
//
//            double A = 3.262 + 14.98 * pow(pr, 5.508);
//            double B = 1.39 + 5.746 * pr;
//            double Z2 = 0.6 + 0.76 * pow(pr, A) + (6.99 * pow(pr, B) - 0.6) * (1.0 - tr);
//            double ksi = 0.176 * pow(data.criticalTemperature.value(), 1.0/6) * pow(data.molecularWeight.value(), -0.5) * pow(data.criticalPressure.value() / 1E5, -2.0/3);
//            double eta_ig = data.saturatedVaporViscosityCorrelation(temperature);
//            double Fp = (1.0 + (Fp_ig - 1.0) * pow(Z2/(ksi*eta_ig/1.0E-7), -3)) / Fp_ig;
//            return Z2 * Fp * 1E-7 / ksi;
//        }
//
//        if (tr <= 40.0 && pr <= 100.0) {
//
//            double A = 0.001245 / tr * exp(5.1726 * pow(tr, -0.3286));
//            double B = A * (1.6553 * tr - 1.2723);
//            double C = 0.4489 / tr * exp(3.0578 * pow(tr, -37.7332));
//            double D = 1.7368 / tr * exp(2.2310 * pow(tr, -7.6351));
//            double E = 1.3088;
//            double F = 0.9425 * exp(-0.1853 * pow(tr, 0.4489));
//            double Z2 = 1.0 + (A * pow(pr, E)) / (B * pow(pr, F) + pow(1.0 + C * pow(pr, D), -1));
//            double eta_ig = data.saturatedVaporViscosityCorrelation(temperature);
//            double Fp = (1.0 + (Fp_ig - 1.0) * pow(Z2, -3)) / Fp_ig;
//            return eta_ig * Z2 * Fp;
//        }
//
//        throw PCProps::PCPropsException("Lucas Viscosity Estimation Error: Invalid temperature/pressure range.");
//
//    }

}    // namespace


namespace PCProps
{
    PCFluid::PCFluid() = default;

    PCFluid::PCFluid(const PCPureComponent& pc, const PCEquationOfState& eos) : m_pureComponent{pc}, m_equationOfState{eos} {

        nlohmann::json obj;
        obj["Tc"] = m_pureComponent.criticalTemperature();
        obj["Pc"] = m_pureComponent.criticalPressure();
        obj["Omega"] = m_pureComponent.acentricFactor();

        m_equationOfState.setProperties(obj.dump());
        m_equationOfState.setIdealGasCpFunction([&](double t) {return m_pureComponent.idealGasCp(t);} );
    }

    PCFluid::PCFluid(const PCFluid& other) = default;

    PCFluid::PCFluid(PCFluid&& other) noexcept = default;

    PCFluid::~PCFluid() = default;

    PCFluid& PCFluid::operator=(const PCFluid& other) = default;

    PCFluid& PCFluid::operator=(PCFluid&& other) noexcept = default;

    PCPhases PCFluid::flash(Pressure pressure, Temperature temperature) const
    {
        using std::get;
        auto results = PCPhase(m_equationOfState.flashPT(pressure.get(), temperature.get())[0]);
        switch (determinePhaseType(
            pressure.get(),
            m_pureComponent.criticalPressure(),
            m_pureComponent.satLiquidVolume(temperature.get()),
            temperature.get(),
            m_pureComponent.criticalTemperature()))
        {
            case PhaseType::Liquid:
                results.setMolarVolume(m_pureComponent.satLiquidVolume(temperature.get()));
                break;

            case PhaseType::Vapor:
            case PhaseType::Dense:
                break;

            default:
                throw std::runtime_error("Something went wrong. Invalid phase properties");
        }

        results.setSurfaceTension(0.0);
        results.setThermalConductivity(0.0);
        results.setMolarWeight(m_pureComponent.molarWeight());

        PCPhases result = { results.data() };
//        computeViscosity(result);

        return { result };
    }

    PCPhases PCFluid::flash(Pressure pressure, VaporFraction vaporFraction) const
    {
        using std::get;
        PCPhases results;

        for (auto& phase : m_equationOfState.flashPx(pressure.get(), vaporFraction.get())) {
            auto result = PCPhase(phase);
            result.setSurfaceTension(0.0);
            result.setThermalConductivity(0.0);
            result.setMolarWeight(m_pureComponent.molarWeight());
            results.emplace_back(result.data());
        }

//        computeViscosity(results);
        return results;
    }

    PCPhases PCFluid::flash(Temperature temperature, VaporFraction vaporFraction) const
    {
        using std::get;
        PCPhases results;

        for (auto& phase : m_equationOfState.flashTx(temperature.get(), vaporFraction.get())) {
            auto result = PCPhase(phase);
            result.setSurfaceTension(0.0);
            result.setThermalConductivity(0.0);
            result.setMolarWeight(m_pureComponent.molarWeight());
            results.emplace_back(result.data());
        }

//        computeViscosity(results);
        return results;
    }

    PCPhases PCFluid::flash(Pressure pressure, Enthalpy enthalpy) const
    {
        using std::get;
        PCPhases results;

        for (auto& phase : m_equationOfState.flashPH(pressure.get(), enthalpy.get())) {
            auto result = PCPhase(phase);
            result.setSurfaceTension(0.0);
            result.setThermalConductivity(0.0);
            result.setMolarWeight(m_pureComponent.molarWeight());
            results.emplace_back(result.data());
        }

//        computeViscosity(results);
        return results;
    }

    PCPhases PCFluid::flash(Pressure pressure, Entropy entropy) const
    {
        using std::get;
        PCPhases results;

        for (auto& phase : m_equationOfState.flashPS(pressure.get(), entropy.get())) {
            auto result = PCPhase(phase);
            result.setSurfaceTension(0.0);
            result.setThermalConductivity(0.0);
            result.setMolarWeight(m_pureComponent.molarWeight());
            results.emplace_back(result.data());
        }

//        computeViscosity(results);
        return results;
    }

    PCPhases PCFluid::flash(Temperature temperature, MolarVolume volume) const
    {
        return PCProps::PCPhases();
    }

}