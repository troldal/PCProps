//
// Created by Kenneth Balslev on 18/01/2021.
//

#include "Fluid.hpp"

#include <PropertyLib.hpp>

#include <json/json.hpp>

#include <stdexcept>

namespace PCProps
{
    class Fluid::impl
    {
        IPureComponent   m_pureComponent {};
        IEquationOfState m_equationOfState {};

        std::function<double(double, double)> m_compressedLiquidVolume;
        std::function<double(double, double)> m_compressedLiquidViscosity;
        std::function<double(double, double)> m_compressedVaporViscosity;

        enum PhaseType { Liquid, Vapor, Dense, Undefined };

        PhaseType determinePhaseType(const PCProps::PCPhase& phase) const
        {
            PhaseType type = PhaseType::Undefined;

            auto tc = m_pureComponent.criticalTemperature();
            auto pc = m_pureComponent.criticalPressure();
            auto t = phase.temperature();
            auto p = phase.pressure();
            auto x = phase.molarFraction();
            auto z = phase.compressibility();
            auto psat = m_equationOfState.saturationPressure(t);

            if (t > tc && p > pc)
                return PhaseType::Dense;
            if (x < 1.0 && z > 0.5)
                return PhaseType::Vapor;
            if (x < 1.0 && z < 0.5)
                return PhaseType::Liquid;
            if ((t > tc && p <= pc) || (t <= tc && p <= psat))
                return PhaseType::Vapor;
            if ((x < 1.0 && z < 0.5) || (t <= tc && p > psat))
                return PhaseType::Liquid;

            return PhaseType::Undefined;
        }

        void computePhaseProperties(PCProps::PCPhase& phase) const
        {
            auto type = determinePhaseType(phase);
            auto t = phase.temperature();
            auto p = phase.pressure();

            phase.setMolarWeight(m_pureComponent.molarWeight());

            switch (type)
            {
                case PhaseType::Liquid:
                    phase.setMolarVolume(m_compressedLiquidVolume(t, p));
                    phase.setSurfaceTension(0.0);
                    phase.setThermalConductivity(0.0);
                    phase.setViscosity(m_compressedLiquidViscosity(t, p));

                    break;

                case PhaseType::Vapor:
                    phase.setSurfaceTension(0.0);
                    phase.setThermalConductivity(0.0);
                    phase.setViscosity(m_compressedVaporViscosity(t, p));
                    break;

                case PhaseType::Dense:
                    break;

                default:
                    throw std::runtime_error("Something went wrong. Invalid phase properties");
            }
        }

    public:
        impl(const IPureComponent& pc, const IEquationOfState& eos) : m_pureComponent { pc }, m_equationOfState { eos }
        {
            nlohmann::json obj;
            obj["Tc"]    = m_pureComponent.criticalTemperature();
            obj["Pc"]    = m_pureComponent.criticalPressure();
            obj["Omega"] = m_pureComponent.acentricFactor();

            m_equationOfState.setProperties(obj.dump());
            m_equationOfState.setIdealGasCpFunction([&](double t) { return m_pureComponent.idealGasCp(t); });

            using PCProps::LiquidVolume::Thomson;
            m_compressedLiquidVolume = Thomson(
                m_pureComponent.criticalTemperature(),
                m_pureComponent.criticalPressure(),
                m_pureComponent.acentricFactor(),
                [&](double t) { return m_pureComponent.satLiquidVolume(t); },
                [&](double t) { return m_equationOfState.saturationPressure(t); });

            using namespace PCProps::CompressedLiquidViscosity;
            m_compressedLiquidViscosity = CompressedLiquidViscosity::Lucas(
                m_pureComponent.criticalTemperature(),
                m_pureComponent.criticalPressure(),
                m_pureComponent.acentricFactor(),
                [&](double t) { return m_pureComponent.satLiquidViscosity(t); },
                [&](double t) { return m_equationOfState.saturationPressure(t); });

            using namespace PCProps::CompressedVaporViscosity;
            m_compressedVaporViscosity = CompressedVaporViscosity::Lucas(
                m_pureComponent.criticalTemperature(),
                m_pureComponent.criticalPressure(),
                m_pureComponent.criticalCompressibility(),
                m_pureComponent.molarWeight(),
                m_pureComponent.dipoleMoment(),
                [&](double t) { return m_pureComponent.satVaporViscosity(t); },
                [&](double t) { return m_equationOfState.saturationPressure(t); });
        }

        PCPhases flash(Pressure pressure, Temperature temperature) const
        {
            PCPhases results;
            for (auto& phase : m_equationOfState.flashPT(pressure.get(), temperature.get())) {
                auto result = PCPhase(phase);
                computePhaseProperties(result);
                results.emplace_back(result.data());
            }

            return results;
        }

        PCPhases flash(Pressure pressure, VaporFraction vaporFraction) const
        {
            PCPhases results;
            for (auto& phase : m_equationOfState.flashPx(pressure.get(), vaporFraction.get())) {
                auto result = PCPhase(phase);
                computePhaseProperties(result);
                results.emplace_back(result.data());
            }

            return results;
        }

        PCPhases flash(Temperature temperature, VaporFraction vaporFraction) const
        {
            PCPhases results;
            for (auto& phase : m_equationOfState.flashTx(temperature.get(), vaporFraction.get())) {
                auto result = PCPhase(phase);
                computePhaseProperties(result);
                results.emplace_back(result.data());
            }

            return results;
        }

        PCPhases flash(Pressure pressure, Enthalpy enthalpy) const
        {
            PCPhases results;
            for (auto& phase : m_equationOfState.flashPH(pressure.get(), enthalpy.get())) {
                auto result = PCPhase(phase);
                computePhaseProperties(result);
                results.emplace_back(result.data());
            }

            return results;
        }

        PCPhases flash(Pressure pressure, Entropy entropy) const
        {
            PCPhases results;
            for (auto& phase : m_equationOfState.flashPS(pressure.get(), entropy.get())) {
                auto result = PCPhase(phase);
                computePhaseProperties(result);
                results.emplace_back(result.data());
            }

            return results;
        }

        PCPhases flash(Temperature temperature, MolarVolume volume) const {

            return PCProps::PCPhases();
        }
    };

    // =====================================================================
    // PUBLIC INTERFACE
    // =====================================================================

    Fluid::Fluid() = default;

    Fluid::Fluid(const IPureComponent& pc, const IEquationOfState& eos) : m_impl(std::make_unique<impl>(pc, eos)) {}

    Fluid::Fluid(const Fluid& other) : m_impl(std::make_unique<impl>(*other.m_impl)) {};

    Fluid::Fluid(Fluid&& other) noexcept = default;

    Fluid::~Fluid() = default;

    Fluid& Fluid::operator=(const Fluid& other)
    {
        Fluid copy = other;
        *this      = std::move(copy);
        return *this;
    };

    Fluid& Fluid::operator=(Fluid&& other) noexcept = default;

    PCPhases Fluid::flash(Pressure pressure, Temperature temperature) const
    {
        return m_impl->flash(pressure, temperature);
    }

    PCPhases Fluid::flash(Pressure pressure, VaporFraction vaporFraction) const
    {
        return m_impl->flash(pressure, vaporFraction);
    }

    PCPhases Fluid::flash(Temperature temperature, VaporFraction vaporFraction) const
    {
        return m_impl->flash(temperature, vaporFraction);
    }

    PCPhases Fluid::flash(Pressure pressure, Enthalpy enthalpy) const
    {
        return m_impl->flash(pressure, enthalpy);
    }

    PCPhases Fluid::flash(Pressure pressure, Entropy entropy) const
    {
        return m_impl->flash(pressure, entropy);
    }

    PCPhases Fluid::flash(Temperature temperature, MolarVolume volume) const
    {
        return m_impl->flash(temperature, volume);
    }

}    // namespace PCProps