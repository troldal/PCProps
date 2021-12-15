//
// Created by Kenneth Balslev on 18/01/2021.
//

#include "Fluid.hpp"

#include <FluidProperties.hpp>
#include <PhaseProperties.hpp>

#include <stdexcept>
#include <tuple>

using JSONString = std::string;

namespace PCProps
{
    class Fluid::impl
    {
        // ===== Objects representing the pure component and equation of state.
        IPureComponent          m_pureComponent {};
        IEquationOfState        m_equationOfState {};
        mutable FluidProperties m_phaseProps {};

        /**
         * @brief Compute the liquid properties for a phase.
         * @param liquid The (liquid) phase for which to compute properties.
         * @return The input phase updated with liquid properties.
         */
        PhaseProperties computeLiquidProperties(PhaseProperties liquid) const
        {
//            liquid.MolarVolume = m_pureComponent.correlation(
//                "CompressedLiquidVolume",
//                { liquid.Temperature, liquid.Pressure, liquid.VaporPressure, m_pureComponent.correlation("SaturatedLiquidVolume", liquid.Temperature) });
//            liquid.SurfaceTension      = 0.0;
//            liquid.ThermalConductivity = 0.0;
//            liquid.Viscosity           = m_pureComponent.correlation(
//                "CompressedLiquidViscosity",
//                { liquid.Temperature, liquid.Pressure, liquid.VaporPressure, m_pureComponent.correlation("SaturatedLiquidViscosity", liquid.Temperature) });
//            //liquid.Cp = m_pureComponent.correlation("LiquidCp", liquid.Temperature);

            return liquid;
        }

        /**
         * @brief Compute the vapor properties for a phase.
         * @param vapor The (vapor) phase for which to compute properties.
         * @return The input phase updated with vapor properties.
         */
        PhaseProperties computeVaporProperties(PhaseProperties vapor) const
        {
//            vapor.SurfaceTension      = 0.0;
//            vapor.ThermalConductivity = 0.0;
//            vapor.Viscosity           = m_pureComponent.correlation(
//                "CompressedVaporViscosity",
//                { vapor.Temperature, vapor.Pressure, vapor.VaporPressure, m_pureComponent.correlation("SaturatedVaporViscosity", vapor.Temperature) });

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

                switch (phase.Type) {
                    case PhaseType::Liquid:
                        phase = computeLiquidProperties(phase);
                        break;

                    case PhaseType::Vapor:
                        phase = computeVaporProperties(phase);
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
            m_phaseProps = FluidProperties(m_equationOfState.flash("PT", pressure, temperature));
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
            m_phaseProps = FluidProperties(m_equationOfState.flash("Px",pressure, vaporFraction));
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
            m_phaseProps = FluidProperties(m_equationOfState.flash("Tx", temperature, vaporFraction));
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
            m_phaseProps = FluidProperties(m_equationOfState.flash("PH", pressure, enthalpy));
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
            m_phaseProps = FluidProperties(m_equationOfState.flash("PS", pressure, entropy));
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
            m_phaseProps = FluidProperties(m_equationOfState.flash("TV", temperature, volume));
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

    JSONString Fluid::flash(const std::string& specification, double var1, double var2) const
    {
        if (specification == "PT")
            return m_impl->flashPT(var1, var2).asJSON();

        if (specification == "Px")
            return m_impl->flashPx(var1, var2).asJSON();

        if (specification == "Tx")
            return m_impl->flashTx(var1, var2).asJSON();

        if (specification == "PH")
            return m_impl->flashPH(var1, var2).asJSON();

        if (specification == "PS")
            return m_impl->flashPS(var1, var2).asJSON();

        if (specification == "TV")
            return m_impl->flashTV(var1, var2).asJSON();
    }

    /**
     *
     */
    JSONString Fluid::properties() const
    {
        return m_impl->getProperties().asJSON();
    }

}    // namespace PCProps