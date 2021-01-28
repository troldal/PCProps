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

        std::function<double(double, double, double, double)> m_compressedLiquidVolume;
        std::function<double(double, double, double, double)> m_compressedLiquidViscosity;
        std::function<double(double, double, double, double)> m_compressedVaporViscosity;

        mutable PCPhases m_phaseData {};

        enum PhaseType { Liquid, Vapor, Dense, Undefined };

        PhaseType determinePhaseType(const PCProps::PCPhase& phase) const
        {
            PhaseType type = PhaseType::Undefined;

            auto tc = m_pureComponent.criticalTemperature();
            auto pc = m_pureComponent.criticalPressure();
            auto t = phase[PCTemperature];
            auto p = phase[PCPressure];
            auto x = phase[PCMolarFlow];
            auto z = phase[PCCompressibility];
            auto psat = phase[PCVaporPressure]; //   m_equationOfState.saturationPressure(t);

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
            auto t = phase[PCTemperature];
            auto p = phase[PCPressure];

            phase[PCMolarWeight] = m_pureComponent.molarWeight();

            switch (type)
            {
                case PhaseType::Liquid:
                    phase[PCMolarWeight] = m_compressedLiquidVolume(t, p, m_equationOfState.saturationPressure(t), m_pureComponent.satLiquidVolume(t));
                    phase[PCSurfaceTension] = 0.0;
                    phase[PCThermalConductivity] = 0.0;
                    phase[PCViscosity] = m_compressedLiquidViscosity(t, p, m_equationOfState.saturationPressure(t), m_pureComponent.satLiquidViscosity(t));;

                    break;

                case PhaseType::Vapor:
                    phase[PCSurfaceTension] = 0.0;
                    phase[PCThermalConductivity] = 0.0;
                    phase[PCViscosity] = m_compressedVaporViscosity(t, p, m_equationOfState.saturationPressure(t), m_pureComponent.satVaporViscosity(t));
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

            using PCProps::HeatCapacity::AlyLee;
            m_equationOfState.setIdealGasCpFunction([pc=m_pureComponent](double t) { return pc.idealGasCp(t); } );

            using PCProps::LiquidVolume::Thomson;
            m_compressedLiquidVolume = Thomson(
                m_pureComponent.criticalTemperature(),
                m_pureComponent.criticalPressure(),
                m_pureComponent.acentricFactor());

            using namespace PCProps::CompressedLiquidViscosity;
            m_compressedLiquidViscosity = CompressedLiquidViscosity::Lucas(
                m_pureComponent.criticalTemperature(),
                m_pureComponent.criticalPressure(),
                m_pureComponent.acentricFactor());

            using namespace PCProps::CompressedVaporViscosity;
            m_compressedVaporViscosity = CompressedVaporViscosity::Lucas(
                m_pureComponent.criticalTemperature(),
                m_pureComponent.criticalPressure(),
                m_pureComponent.criticalCompressibility(),
                m_pureComponent.molarWeight(),
                m_pureComponent.dipoleMoment());
        }

        const PCPhases& flashPT(double pressure, double temperature) const
        {
            PCPhases results;
            auto temp = m_equationOfState.flashPT(pressure, temperature);
            for (auto& phase : temp) {
                auto result = PCPhase(phase);
                computePhaseProperties(result);
                results.emplace_back(result);
            }

            m_phaseData = std::move(results);
            return m_phaseData;
        }

        const PCPhases& flashPx(double pressure, double vaporFraction) const
        {
            PCPhases results;
            for (auto& phase : m_equationOfState.flashPx(pressure, vaporFraction)) {
                auto result = PCPhase(phase);
                computePhaseProperties(result);
                results.emplace_back(result);
            }

            m_phaseData = std::move(results);
            return m_phaseData;
        }

        const PCPhases& flashTx(double temperature, double vaporFraction) const
        {
            PCPhases results;
            for (auto& phase : m_equationOfState.flashTx(temperature, vaporFraction)) {
                auto result = PCPhase(phase);
                computePhaseProperties(result);
                results.emplace_back(result);
            }

            m_phaseData = std::move(results);
            return m_phaseData;
        }

        const PCPhases& flashPH(double pressure, double enthalpy) const
        {
            PCPhases results;
            for (auto& phase : m_equationOfState.flashPH(pressure, enthalpy)) {
                auto result = PCPhase(phase);
                computePhaseProperties(result);
                results.emplace_back(result);
            }

            m_phaseData = std::move(results);
            return m_phaseData;
        }

        const PCPhases& flashPS(double pressure, double entropy) const
        {
            PCPhases results;
            for (auto& phase : m_equationOfState.flashPS(pressure, entropy)) {
                auto result = PCPhase(phase);
                computePhaseProperties(result);
                results.emplace_back(result);
            }

            m_phaseData = std::move(results);
            return m_phaseData;
        }

        const PCPhases& flashTV(double temperature, double volume) const {

            return m_phaseData;
        }

        const PCPhases& getProperties() const {
            return m_phaseData;
        }
    };

    // =====================================================================
    // PUBLIC INTERFACE
    // =====================================================================

    Fluid::Fluid() = default;

    Fluid::Fluid(const IPureComponent& pc, const IEquationOfState& eos) : m_impl(std::make_unique<impl>(pc, eos)) {
        int i = 0;
    }

    Fluid::Fluid(const Fluid& other) : m_impl(std::make_unique<impl>(*other.m_impl)) {
        int i = 0;
    }

    Fluid::Fluid(Fluid&& other) noexcept = default;

    Fluid::~Fluid() = default;

    Fluid& Fluid::operator=(const Fluid& other)
    {
        Fluid copy = other;
        *this      = std::move(copy);
        return *this;
    };

    Fluid& Fluid::operator=(Fluid&& other) noexcept = default;

    const PCPhases& Fluid::flashPT(double pressure, double temperature) const
    {
        return m_impl->flashPT(pressure, temperature);
    }

    const PCPhases& Fluid::flashPx(double pressure, double vaporFraction) const
    {
        return m_impl->flashPx(pressure, vaporFraction);
    }

    const PCPhases& Fluid::flashTx(double temperature, double vaporFraction) const
    {
        return m_impl->flashTx(temperature, vaporFraction);
    }

    const PCPhases& Fluid::flashPH(double pressure, double enthalpy) const
    {
        return m_impl->flashPH(pressure, enthalpy);
    }

    const PCPhases& Fluid::flashPS(double pressure, double entropy) const
    {
        return m_impl->flashPS(pressure, entropy);
    }

    const PCPhases& Fluid::flashTV(double temperature, double volume) const
    {
        return m_impl->flashTV(temperature, volume);
    }

    const PCPhases& Fluid::properties() const
    {
        return m_impl->getProperties();
    }

}    // namespace PCProps