//
// Created by Kenneth Balslev on 18/01/2021.
//

#include "Fluid.hpp"

#include <PhaseProperties.hpp>
#include <json/json.hpp>

#include <stdexcept>
#include <tuple>

using JSONString = std::string;

namespace PCProps
{
    class Fluid::impl
    {
        // ===== Objects representing the pure component and equation of state.
        IPureComponent                       m_pureComponent {};
        IEquationOfState                     m_equationOfState {};
        mutable std::vector<PhaseProperties> m_phaseProps {};

        // ===== Enum class used in the determinePhaseType function.
        enum PhaseType { Liquid, Vapor, Dense, Undefined };

        /**
         * @brief Determine the phase type (liquid, vapor, etc), based on the phase properties.
         * @param phase The phase for which to determine the type.
         * @return A PhaseType enum representing the phase type.
         */
        PhaseType determinePhaseType(const PhaseProperties& phase) const
        {
            auto tc   = m_pureComponent.property("CriticalTemperature");
            auto pc   = m_pureComponent.property("CriticalPressure");
            auto t    = phase.Temperature;
            auto p    = phase.Pressure;
            auto x    = phase.MolarFlow;
            auto z    = phase.Compressibility;
            auto psat = phase.VaporPressure;

            if (t > tc && p > pc) return PhaseType::Dense;
            if (x < 1.0 && z > 0.5) return PhaseType::Vapor;
            if (x < 1.0 && z < 0.5) return PhaseType::Liquid;
            if ((t > tc && p <= pc) || (t <= tc && p <= psat)) return PhaseType::Vapor;
            if ((x < 1.0 && z < 0.5) || (t <= tc && p > psat)) return PhaseType::Liquid;

            return PhaseType::Undefined;
        }

        /**
         * @brief Compute the liquid properties for a phase.
         * @param liquid The (liquid) phase for which to compute properties.
         * @return The input phase updated with liquid properties.
         */
        PhaseProperties computeLiquidProperties(PhaseProperties liquid) const
        {
            liquid.MolarVolume = m_pureComponent.property(
                "CompressedLiquidVolume",
                { liquid.Temperature, liquid.Pressure, liquid.VaporPressure, m_pureComponent.property("SaturatedLiquidVolume", liquid.Temperature) });
            liquid.SurfaceTension      = 0.0;
            liquid.ThermalConductivity = 0.0;
            liquid.Viscosity = m_pureComponent.property(
                "CompressedLiquidViscosity",
                { liquid.Temperature, liquid.Pressure, liquid.VaporPressure, m_pureComponent.property("SaturatedLiquidViscosity", liquid.Temperature) });

            return liquid;
        }

        /**
         * @brief Compute the vapor properties for a phase.
         * @param vapor The (vapor) phase for which to compute properties.
         * @return The input phase updated with vapor properties.
         */
        PhaseProperties computeVaporProperties(PhaseProperties vapor) const
        {
            vapor.SurfaceTension      = 0.0;
            vapor.ThermalConductivity = 0.0;
            vapor.Viscosity = m_pureComponent.property(
                         "CompressedVaporViscosity",
                         { vapor.Temperature, vapor.Pressure, vapor.VaporPressure, m_pureComponent.property("SaturatedVaporViscosity", vapor.Temperature) });

            return vapor;
        }

        /**
         * @brief Compute properties for all phases in the fluid.
         */
        void computePhaseProperties() const
        {
            for (auto& phase : m_phaseProps) {
                phase.MolarWeight = m_pureComponent.property("MolarWeight");

                switch (determinePhaseType(phase)) {
                    case PhaseType::Liquid:
                        phase = computeLiquidProperties(phase);
                        break;

                    case PhaseType::Vapor:
                        phase = computeVaporProperties(phase);
                        break;

                    case PhaseType::Dense:
                        break;

                    default:
                        throw std::runtime_error("Something went wrong. Invalid phase properties");
                }
            }
        }

    public:
        /**
         *
         * @param pc
         * @param eos
         */
        impl(const IPureComponent& pc, const IEquationOfState& eos) : m_pureComponent { pc }, m_equationOfState { eos }
        {
            m_equationOfState.init(m_pureComponent);
        }

        /**
         *
         * @param pressure
         * @param temperature
         * @return
         */
        const std::vector<PhaseProperties>& flashPT(double pressure, double temperature) const
        {
            m_phaseProps.clear();
            auto temp = nlohmann::json::parse(m_equationOfState.flashPT(pressure, temperature));
            for (auto& phase : temp) m_phaseProps.emplace_back(phase.dump());

            computePhaseProperties();
            return m_phaseProps;
        }

        /**
         *
         * @param pressure
         * @param vaporFraction
         * @return
         */
        const std::vector<PhaseProperties>& flashPx(double pressure, double vaporFraction) const
        {
            m_phaseProps.clear();
            auto temp = nlohmann::json::parse(m_equationOfState.flashPx(pressure, vaporFraction));
            for (auto& phase : temp) m_phaseProps.emplace_back(phase.dump());

            computePhaseProperties();
            return m_phaseProps;
        }

        /**
         *
         * @param temperature
         * @param vaporFraction
         * @return
         */
        const std::vector<PhaseProperties>& flashTx(double temperature, double vaporFraction) const
        {
            m_phaseProps.clear();
            auto temp = nlohmann::json::parse(m_equationOfState.flashTx(temperature, vaporFraction));
            for (auto& phase : temp) m_phaseProps.emplace_back(phase.dump());

            computePhaseProperties();
            return m_phaseProps;
        }

        /**
         *
         * @param pressure
         * @param enthalpy
         * @return
         */
        const std::vector<PhaseProperties>& flashPH(double pressure, double enthalpy) const
        {
            m_phaseProps.clear();
            auto temp = nlohmann::json::parse(m_equationOfState.flashPH(pressure, enthalpy));
            for (auto& phase : temp) m_phaseProps.emplace_back(phase.dump());

            computePhaseProperties();
            return m_phaseProps;
        }

        /**
         *
         * @param pressure
         * @param entropy
         * @return
         */
        const std::vector<PhaseProperties>& flashPS(double pressure, double entropy) const
        {
            m_phaseProps.clear();
            auto temp = nlohmann::json::parse(m_equationOfState.flashPS(pressure, entropy));
            for (auto& phase : temp) m_phaseProps.emplace_back(phase.dump());

            computePhaseProperties();
            return m_phaseProps;
        }

        /**
         *
         * @param temperature
         * @param volume
         * @return
         */
        const std::vector<PhaseProperties>& flashTV(double temperature, double volume) const
        {
            m_phaseProps.clear();
            auto temp = nlohmann::json::parse(m_equationOfState.flashTV(temperature, volume));
            for (auto& phase : temp) m_phaseProps.emplace_back(phase.dump());

            computePhaseProperties();
            return m_phaseProps;
        }

        /**
         *
         * @return
         */
        const std::vector<PhaseProperties>& getProperties() const
        {
            return m_phaseProps;
        }
    };

    // =====================================================================
    // PUBLIC INTERFACE
    // =====================================================================

    Fluid::Fluid() = default;

    Fluid::Fluid(const IPureComponent& pc, const IEquationOfState& eos) : m_impl(std::make_unique<impl>(pc, eos)) {}

    Fluid::Fluid(const Fluid& other) : m_impl(std::make_unique<impl>(*other.m_impl)) {}

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
        for (const auto& phase : m_impl->flashPT(pressure, temperature)) result.emplace_back(nlohmann::json::parse(phase.asJSON()));
        return nlohmann::json(result).dump();
    }

    JSONString Fluid::flashPx(double pressure, double vaporFraction) const
    {
        std::vector<nlohmann::json> result;
        for (const auto& phase : m_impl->flashPx(pressure, vaporFraction)) result.emplace_back(nlohmann::json::parse(phase.asJSON()));
        return nlohmann::json(result).dump();
    }

    JSONString Fluid::flashTx(double temperature, double vaporFraction) const
    {
        std::vector<nlohmann::json> result;
        for (const auto& phase : m_impl->flashTx(temperature, vaporFraction)) result.emplace_back(nlohmann::json::parse(phase.asJSON()));
        return nlohmann::json(result).dump();
    }

    JSONString Fluid::flashPH(double pressure, double enthalpy) const
    {
        std::vector<nlohmann::json> result;
        for (const auto& phase : m_impl->flashPH(pressure, enthalpy)) result.emplace_back(nlohmann::json::parse(phase.asJSON()));
        return nlohmann::json(result).dump();
    }

    JSONString Fluid::flashPS(double pressure, double entropy) const
    {
        std::vector<nlohmann::json> result;
        for (const auto& phase : m_impl->flashPS(pressure, entropy)) result.emplace_back(nlohmann::json::parse(phase.asJSON()));
        return nlohmann::json(result).dump();
    }

    JSONString Fluid::flashTV(double temperature, double volume) const
    {
        std::vector<nlohmann::json> result;
        for (const auto& phase : m_impl->flashTV(temperature, volume)) result.emplace_back(nlohmann::json::parse(phase.asJSON()));
        return nlohmann::json(result).dump();
    }

    JSONString Fluid::properties() const
    {
        std::vector<nlohmann::json> result;
        for (const auto& phase : m_impl->getProperties()) result.emplace_back(nlohmann::json::parse(phase.asJSON()));
        return nlohmann::json(result).dump();
    }

}    // namespace PCProps