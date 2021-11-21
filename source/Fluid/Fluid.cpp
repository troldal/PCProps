//
// Created by Kenneth Balslev on 18/01/2021.
//

#include "Fluid.hpp"

#include <PhaseProperties.hpp>
#include <FluidProperties.hpp>

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
        mutable FluidProperties m_phaseProps {};

        // ===== Enum class used in the determinePhaseType function.
        enum PhaseType { Liquid, Vapor, Dense, Undefined };

        /**
         * @brief Determine the phase type (liquid, vapor, etc), based on the phase properties.
         * @param phase The phase for which to determine the type.
         * @return A PhaseType enum representing the phase type.
         * @todo How will this work for multi component mixtures?
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
            auto phases = m_phaseProps.phases();
            for (auto& phase : phases) {
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

            m_phaseProps = phases;
        }

    public:
        /**
         *
         * @param pureComponent
         * @param eos
         */
        impl(const IPureComponent& pureComponent, const IEquationOfState& eos) : m_pureComponent { pureComponent }, m_equationOfState { eos }
        {
            m_equationOfState.init(m_pureComponent);
        }

        /**
         *
         * @param pressure
         * @param temperature
         * @return
         */
        const FluidProperties& flashPT(double pressure, double temperature) const
        {
            m_phaseProps = FluidProperties(m_equationOfState.flashPT(pressure, temperature));
            computePhaseProperties();
            return m_phaseProps;
        }

        /**
         *
         * @param pressure
         * @param vaporFraction
         * @return
         */
        const FluidProperties& flashPx(double pressure, double vaporFraction) const
        {
            m_phaseProps = FluidProperties(m_equationOfState.flashPx(pressure, vaporFraction));
            computePhaseProperties();
            return m_phaseProps;
        }

        /**
         *
         * @param temperature
         * @param vaporFraction
         * @return
         */
        const FluidProperties& flashTx(double temperature, double vaporFraction) const
        {
            m_phaseProps = FluidProperties(m_equationOfState.flashTx(temperature, vaporFraction));
            computePhaseProperties();
            return m_phaseProps;
        }

        /**
         *
         * @param pressure
         * @param enthalpy
         * @return
         */
        const FluidProperties& flashPH(double pressure, double enthalpy) const
        {
            m_phaseProps = FluidProperties(m_equationOfState.flashPH(pressure, enthalpy));
            computePhaseProperties();
            return m_phaseProps;
        }

        /**
         *
         * @param pressure
         * @param entropy
         * @return
         */
        const FluidProperties& flashPS(double pressure, double entropy) const
        {
            m_phaseProps = FluidProperties(m_equationOfState.flashPS(pressure, entropy));
            computePhaseProperties();
            return m_phaseProps;
        }

        /**
         *
         * @param temperature
         * @param volume
         * @return
         */
        const FluidProperties& flashTV(double temperature, double volume) const
        {
            m_phaseProps = FluidProperties(m_equationOfState.flashTV(temperature, volume));
            computePhaseProperties();
            return m_phaseProps;
        }

        /**
         *
         * @return
         */
        const FluidProperties& getProperties() const
        {
            return m_phaseProps;
        }
    };

    // =====================================================================
    // PUBLIC INTERFACE
    // =====================================================================

    /**
     * Default constructor
     */
    Fluid::Fluid() = default;

    /**
     * Constructor
     */
    Fluid::Fluid(const IPureComponent& pureComponent, const IEquationOfState& eos) : m_impl(std::make_unique<impl>(pureComponent, eos)) {}

    /**
     * Copy constructor
     */
    Fluid::Fluid(const Fluid& other) : m_impl(std::make_unique<impl>(*other.m_impl)) {}

    /**
     * Move constructor
     */
    Fluid::Fluid(Fluid&& other) noexcept = default;

    /**
     * Destructor
     */
    Fluid::~Fluid() = default;

    /**
     * Copy assignment operator
     */
    Fluid& Fluid::operator=(const Fluid& other)
    {
        Fluid copy = other;
        *this      = std::move(copy);
        return *this;
    };

    /**
     * Move assignment operator
     */
    Fluid& Fluid::operator=(Fluid&& other) noexcept = default;

    /**
     *
     */
    JSONString Fluid::flashPT(double pressure, double temperature) const
    {
        return m_impl->flashPT(pressure, temperature).asJSON();
    }

    /**
     *
     */
    JSONString Fluid::flashPx(double pressure, double vaporFraction) const
    {
        return m_impl->flashPx(pressure, vaporFraction).asJSON();
    }

    /**
     *
     */
    JSONString Fluid::flashTx(double temperature, double vaporFraction) const
    {
        return m_impl->flashTx(temperature, vaporFraction).asJSON();
    }

    /**
     *
     */
    JSONString Fluid::flashPH(double pressure, double enthalpy) const
    {
        return m_impl->flashPH(pressure, enthalpy).asJSON();
    }

    /**
     *
     */
    JSONString Fluid::flashPS(double pressure, double entropy) const
    {
        return m_impl->flashPS(pressure, entropy).asJSON();
    }

    /**
     *
     */
    JSONString Fluid::flashTV(double temperature, double volume) const
    {
        return m_impl->flashTV(temperature, volume).asJSON();
    }

    /**
     *
     */
    JSONString Fluid::properties() const
    {
        return m_impl->getProperties().asJSON();
    }

}    // namespace PCProps