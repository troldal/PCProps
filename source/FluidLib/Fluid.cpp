//
// Created by Kenneth Balslev on 18/01/2021.
//

#include "Fluid.hpp"

#include <PropertyLib.hpp>

#include <json/json.hpp>
#include <common/PhaseProperties.hpp>

#include <stdexcept>
#include <tuple>

using JSONString = std::string;

namespace PCProps
{
    class Fluid::impl
    {

        // ===== Objects representing the pure component and equation of state.
        IPureComponent   m_pureComponent {};
        IEquationOfState m_equationOfState {};

        // ===== Function objects for calculating compressed fluid properties
        std::function<double(double, double, double, double)> m_compressedVaporViscosity;

        mutable std::vector<PhaseProperties> m_phaseProps {};

        enum PhaseType { Liquid, Vapor, Dense, Undefined };
        
        PhaseType determinePhaseType(const PhaseProperties& phase) const
        {
            auto tc = m_pureComponent.criticalTemperature();
            auto pc = m_pureComponent.criticalPressure();
            auto t = phase.Temperature;
            auto p = phase.Pressure;
            auto x = phase.MolarFlow;
            auto z = phase.Compressibility;
            auto psat = phase.VaporPressure; //   m_equationOfState.saturationPressure(t);

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

        void computePhaseProperties() const
        {
            for (auto& phase : m_phaseProps) {
                auto type = determinePhaseType(phase);
                auto t    = phase.Temperature;
                auto p    = phase.Pressure;

                phase.MolarWeight = m_pureComponent.molarWeight();

                switch (type) {
                    case PhaseType::Liquid:
                        phase.MolarVolume = m_pureComponent.compressedLiquidVolume({t, p, m_equationOfState.saturationPressure(t), m_pureComponent.satLiquidVolume(t)});
                        phase.SurfaceTension      = 0.0;
                        phase.ThermalConductivity = 0.0;
                        phase.Viscosity = m_pureComponent.compressedLiquidViscosity({t, p, m_equationOfState.saturationPressure(t), m_pureComponent.satLiquidViscosity(t)});

                        break;

                    case PhaseType::Vapor:
                        phase.SurfaceTension      = 0.0;
                        phase.ThermalConductivity = 0.0;
                        phase.Viscosity = m_pureComponent.compressedVaporViscosity({t, p, m_equationOfState.saturationPressure(t), m_pureComponent.satVaporViscosity(t)});
                        break;

                    case PhaseType::Dense:
                        break;

                    default:
                        throw std::runtime_error("Something went wrong. Invalid phase properties");
                }
            }
        }

    public:
        impl(const IPureComponent& pc, const IEquationOfState& eos) : m_pureComponent { pc }, m_equationOfState { eos }
        {
            m_equationOfState.init(std::make_tuple(
                m_pureComponent.criticalTemperature(),
                m_pureComponent.criticalPressure(),
                m_pureComponent.acentricFactor(),
                [&](double temp){return m_pureComponent.idealGasCp(temp);}));
        }

        const std::vector<PhaseProperties>& flashPT(double pressure, double temperature) const
        {
            m_phaseProps.clear();
            auto temp = nlohmann::json::parse(m_equationOfState.flashPT(pressure, temperature));
            for (auto& phase : temp) m_phaseProps.emplace_back(phase);

            computePhaseProperties();
            return m_phaseProps;
        }

        const std::vector<PhaseProperties>& flashPx(double pressure, double vaporFraction) const
        {
            m_phaseProps.clear();
            auto temp = nlohmann::json::parse(m_equationOfState.flashPx(pressure, vaporFraction));
            for (auto& phase : temp) m_phaseProps.emplace_back(phase);

            computePhaseProperties();
            return m_phaseProps;
        }

        const std::vector<PhaseProperties>& flashTx(double temperature, double vaporFraction) const
        {
            m_phaseProps.clear();
            auto temp = nlohmann::json::parse(m_equationOfState.flashTx(temperature, vaporFraction));
            for (auto& phase : temp) m_phaseProps.emplace_back(phase);

            computePhaseProperties();
            return m_phaseProps;
        }

        const std::vector<PhaseProperties>& flashPH(double pressure, double enthalpy) const
        {
            m_phaseProps.clear();
            auto temp = nlohmann::json::parse(m_equationOfState.flashPH(pressure, enthalpy));
            for (auto& phase : temp) m_phaseProps.emplace_back(phase);

            computePhaseProperties();
            return m_phaseProps;
        }

        const std::vector<PhaseProperties>& flashPS(double pressure, double entropy) const
        {
            m_phaseProps.clear();
            auto temp = nlohmann::json::parse(m_equationOfState.flashPS(pressure, entropy));
            for (auto& phase : temp) m_phaseProps.emplace_back(phase);

            computePhaseProperties();
            return m_phaseProps;
        }

        const std::vector<PhaseProperties>& flashTV(double temperature, double volume) const {

            m_phaseProps.clear();
            auto temp = nlohmann::json::parse(m_equationOfState.flashTV(temperature, volume));
            for (auto& phase : temp) m_phaseProps.emplace_back(phase);

            computePhaseProperties();
            return m_phaseProps;
        }

        const std::vector<PhaseProperties>& getProperties() const {
            return m_phaseProps;
        }
    };

    // =====================================================================
    // PUBLIC INTERFACE
    // =====================================================================

    Fluid::Fluid() = default;

    Fluid::Fluid(const IPureComponent& pc, const IEquationOfState& eos) : m_impl(std::make_unique<impl>(pc, eos)) {
    }

    Fluid::Fluid(const Fluid& other) : m_impl(std::make_unique<impl>(*other.m_impl)) {
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

    JSONString Fluid::flashPT(double pressure, double temperature) const
    {
        std::vector<nlohmann::json> result;
        for (const auto& phase : m_impl->flashPT(pressure, temperature)) result.emplace_back(phase.asJSON());
        return nlohmann::json(result).dump();
    }

    JSONString Fluid::flashPx(double pressure, double vaporFraction) const
    {
        std::vector<nlohmann::json> result;
        for (const auto& phase : m_impl->flashPx(pressure, vaporFraction)) result.emplace_back(phase.asJSON());
        return nlohmann::json(result).dump();
    }

    JSONString Fluid::flashTx(double temperature, double vaporFraction) const
    {
        std::vector<nlohmann::json> result;
        for (const auto& phase : m_impl->flashTx(temperature, vaporFraction)) result.emplace_back(phase.asJSON());
        return nlohmann::json(result).dump();
    }

    JSONString Fluid::flashPH(double pressure, double enthalpy) const
    {
        std::vector<nlohmann::json> result;
        for (const auto& phase : m_impl->flashPH(pressure, enthalpy)) result.emplace_back(phase.asJSON());
        return nlohmann::json(result).dump();
    }

    JSONString Fluid::flashPS(double pressure, double entropy) const
    {
        std::vector<nlohmann::json> result;
        for (const auto& phase : m_impl->flashPS(pressure, entropy)) result.emplace_back(phase.asJSON());
        return nlohmann::json(result).dump();
    }

    JSONString Fluid::flashTV(double temperature, double volume) const
    {
        std::vector<nlohmann::json> result;
        for (const auto& phase : m_impl->flashTV(temperature, volume)) result.emplace_back(phase.asJSON());
        return nlohmann::json(result).dump();
    }

    JSONString Fluid::properties() const
    {
        std::vector<nlohmann::json> result;
        for (const auto& phase : m_impl->getProperties()) result.emplace_back(phase.asJSON());
        return nlohmann::json(result).dump();
    }

}    // namespace PCProps