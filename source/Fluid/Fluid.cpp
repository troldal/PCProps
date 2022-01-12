//
// Created by Kenneth Balslev on 18/01/2021.
//

#include "Fluid.hpp"

#include <FluidProperties.hpp>
#include <PhaseProperties.hpp>
#include <numerics.hpp>
#include <common/Globals.hpp>

#include <stdexcept>
#include <tuple>
#include <cmath>

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

        /**
         * @brief Compute the ideal gas enthalpy at the given T, relative the standard state.
         * @param temperature The temperature [K].
         * @return The ideal gas enthalpy [J/mol]
         */
        inline double idealGasEnthalpy(double temperature) const
        {
            using numeric::integrate;
            using PCProps::Globals::STANDARD_T;
            auto result = integrate([&](double t) { return m_pureComponent.correlation("IdealGasCp", t); }, STANDARD_T, temperature);
            if (std::isnan(result)) throw std::runtime_error("Numeric error: Ideal gas enthalpy could not be computed with T = " + std::to_string(temperature));
            return result;
        }

        /**
         * @brief Compute the ideal gas entropy at the given T and P, relative the standard state.
         * @param temperature The temperature [K].
         * @param pressure The pressure [Pa].
         * @return The ideal gas entropy [J/mol-K]
         */
        inline double idealGasEntropy(double temperature, double pressure) const
        {
            using numeric::integrate;
            using PCProps::Globals::STANDARD_T;
            using PCProps::Globals::STANDARD_P;
            using PCProps::Globals::R_CONST;
            auto result = integrate([&](double temp) { return m_pureComponent.correlation("IdealGasCp", temp) / temp; }, PCProps::Globals::STANDARD_T, temperature) - R_CONST * log(pressure / STANDARD_P);

            if (std::isnan(result))
                throw std::runtime_error("Numeric error: Ideal gas entropy could not be computed with T = " + std::to_string(temperature) + " and P = " + std::to_string(pressure));
            return result;
        }


        void computeProperties() const {
            using Globals::R_CONST;
            auto phases = m_phaseProps.phases();
            for (auto& phase : phases) {
                phase.MolarWeight                 = m_pureComponent.property("MolarWeight");
                phase.Enthalpy                    = idealGasEnthalpy(phase.Temperature) + phase.EnthalpyDeparture;
                phase.Entropy                     = idealGasEntropy(phase.Temperature, phase.Pressure) + phase.EntropyDeparture;
                phase.GibbsEnergy                 = idealGasEnthalpy(phase.Temperature) - phase.Temperature * idealGasEntropy(phase.Temperature, phase.Pressure) + phase.GibbsEnergyDeparture; //data.Enthalpy - temperature * data.Entropy;
                phase.InternalEnergy              = idealGasEnthalpy(phase.Temperature) - phase.Pressure * phase.MolarVolume + phase.InternalEnergyDeparture; //data.Enthalpy - pressure * data.MolarVolume;
                phase.HelmholzEnergy              = phase.InternalEnergy - phase.Temperature * phase.Entropy;
                phase.CriticalPressure            = m_pureComponent.property("CriticalPressure"); //m_criticalPressure;
                phase.CriticalTemperature         = m_pureComponent.property("CriticalTemperature"); //m_criticalTemperature;
                phase.NormalFreezingPoint         = m_pureComponent.property("NormalFreezingPoint"); //m_normalFreezingPoint;
                phase.NormalBoilingPoint          = m_pureComponent.property("NormalBoilingPoint"); //m_normalBoilingPoint;
                phase.Cp                          = m_pureComponent.correlation("IdealGasCp", phase.Temperature) + phase.CpDeparture; //computeCp(temperature, pressure, data.Compressibility);
                phase.Cv                          = m_pureComponent.correlation("IdealGasCp", phase.Temperature) - R_CONST + phase.CvDeparture; //computeCvDeparture(temperature, pressure, data.Compressibility) + data.Cp - R_CONST;
                phase.ThermalExpansionCoefficient = (1.0 / phase.MolarVolume) * phase.DVDT;
                phase.JouleThomsonCoefficient     = -1.0 / (phase.Cp) * (phase.Temperature * phase.DPDT / phase.DPDV + phase.MolarVolume);
                phase.IsothermalCompressibility   = (-1.0 / phase.MolarVolume) * phase.DVDP;
                phase.SpeedOfSound                = (phase.MolarVolume / (16*(-1.0 / phase.MolarVolume) * (phase.Cv / phase.Cp) / phase.DPDV));    // TODO: This calculation does not seem to yield correct results!

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
            m_phaseProps = FluidProperties(m_equationOfState.computeProperties(pressure, temperature));
            computeProperties();
            m_phaseProps = m_phaseProps.stablePhase();
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
            // ===== If the temperature < Tc
            if (temperature < m_pureComponent.property("CriticalTemperature")) {
                // ===== First, calculate the saturation pressure at the specified pressure.
                auto pressure = m_equationOfState.saturationPressure(temperature);
                m_phaseProps = FluidProperties(m_equationOfState.computeProperties(pressure, temperature));
                computeProperties();

                // ===== If the specified vapor fraction is 1.0 (or higher), the fluid is a saturated vapor.
                if (vaporFraction >= 1.0) {
                    m_phaseProps = m_phaseProps.lightPhase();
                    auto phase = m_phaseProps.phases().back();
                    phase.MolarFlow = 1.0;
                    m_phaseProps = { phase };
                    return m_phaseProps;
                }

                // ===== If the specified vapor fraction is 0.0 (or lower), the fluid is a saturated liquid.
                if (vaporFraction <= 0.0) {
                    m_phaseProps = m_phaseProps.heavyPhase();
                    auto phase = m_phaseProps.phases().front();
                    phase.MolarFlow = 1.0;
                    m_phaseProps = { phase };
                    return m_phaseProps;
                }

                // ===== If the vapor fraction is between 0.0 and 1.0, the fluid is two-phase.

                if (m_phaseProps.size() == 1) {
                    auto phase = m_phaseProps.phases().front();
                    phase.MolarFlow = 1.0;
                    m_phaseProps = { phase };
                    return m_phaseProps;
                }

                else {
                    auto heavy = m_phaseProps.phases().front();
                    heavy.MolarFlow = 1.0 - vaporFraction;

                    auto light = m_phaseProps.phases().back();
                    light.MolarFlow = vaporFraction;

                    m_phaseProps = { heavy, light };

                    return m_phaseProps;
                }
            }

            // ===== If the temperature >= Tc, calculate the hypothetical saturation conditions in the supercritical region.
            auto pressure = m_equationOfState.saturationPressure(temperature);
            return flashPT(pressure, temperature);


            //////
//            m_phaseProps = FluidProperties(m_equationOfState.flash("Tx", temperature, vaporFraction));
//            computePhaseProperties();
//            return m_phaseProps;
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