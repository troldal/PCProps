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

        /**
         * @brief
         * @param phases
         * @return
         */
        FluidProperties& computeProperties(FluidProperties& phases) const {
            using Globals::R_CONST;
            for (auto& phase : phases) {
                phase.MolarWeight                 = m_pureComponent.property("MolarWeight");
                phase.Enthalpy                    = idealGasEnthalpy(phase.Temperature) + phase.EnthalpyDeparture;
                phase.Entropy                     = idealGasEntropy(phase.Temperature, phase.Pressure) + phase.EntropyDeparture;
                phase.GibbsEnergy                 = idealGasEnthalpy(phase.Temperature) - phase.Temperature * idealGasEntropy(phase.Temperature, phase.Pressure) + phase.GibbsEnergyDeparture;
                phase.InternalEnergy              = idealGasEnthalpy(phase.Temperature) - phase.Pressure * phase.MolarVolume + phase.InternalEnergyDeparture;
                phase.HelmholzEnergy              = phase.GibbsEnergy - phase.Pressure * phase.MolarVolume;
                phase.CriticalPressure            = m_pureComponent.property("CriticalPressure");
                phase.CriticalTemperature         = m_pureComponent.property("CriticalTemperature");
                phase.NormalFreezingPoint         = m_pureComponent.property("NormalFreezingPoint");
                phase.NormalBoilingPoint          = m_pureComponent.property("NormalBoilingPoint");
                phase.Cp                          = m_pureComponent.correlation("IdealGasCp", phase.Temperature) + phase.CpDeparture;
                phase.Cv                          = m_pureComponent.correlation("IdealGasCp", phase.Temperature) - R_CONST + phase.CvDeparture;
                phase.ThermalExpansionCoefficient = (1.0 / phase.MolarVolume) * phase.DVDT;
                phase.JouleThomsonCoefficient     = -1.0 / (phase.Cp) * (phase.Temperature * phase.DPDT / phase.DPDV + phase.MolarVolume);
                phase.IsothermalCompressibility   = (-1.0 / phase.MolarVolume) * phase.DVDP;
                phase.SpeedOfSound                = (phase.MolarVolume / (phase.MolarWeight * (-1.0 / phase.MolarVolume) * (phase.Cv / phase.Cp) / phase.DPDV));    // TODO: This calculation does not seem to yield correct results!
            }
            return phases;
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
        const FluidProperties flashPT(double pressure, double temperature) const
        {
            auto phaseProps = FluidProperties(m_equationOfState.computePropertiesPT(pressure, temperature));
            return computeProperties(phaseProps).stablePhase();
        }

        /**
         *
         * @param temperature
         * @param vaporFraction
         * @return
         */
        const FluidProperties flashTx(double temperature, double vaporFraction) const
        {
            // ===== If the temperature < Tc
            if (temperature < m_pureComponent.property("CriticalTemperature")) {
                // ===== First, calculate the saturation pressure at the specified pressure.
                auto pressure = m_equationOfState.saturationPressure(temperature);
                auto phaseProps = FluidProperties(m_equationOfState.computePropertiesPT(pressure, temperature));
                computeProperties(phaseProps);

                // ===== If the specified vapor fraction is 1.0 (or higher), the fluid is a saturated vapor.
                if (vaporFraction >= 1.0) {
                    phaseProps = phaseProps.lightPhase();
                    phaseProps.back().MolarFlow = 1.0;
                    return phaseProps;
                }

                // ===== If the specified vapor fraction is 0.0 (or lower), the fluid is a saturated liquid.
                if (vaporFraction <= 0.0) {
                    phaseProps = phaseProps.heavyPhase();
                    phaseProps.front().MolarFlow = 1.0;
                    return phaseProps;
                }

                // ===== If the vapor fraction is between 0.0 and 1.0, the fluid is two-phase.

                if (phaseProps.size() == 1) {
                    phaseProps.front().MolarFlow = 1.0;
                    return phaseProps;
                }

                else {
                    phaseProps.front().MolarFlow = 1.0 - vaporFraction;
                    phaseProps.back().MolarFlow = vaporFraction;
                    return phaseProps;
                }
            }

            // ===== If the temperature >= Tc, calculate the hypothetical saturation conditions in the supercritical region.
            auto pressure = m_equationOfState.saturationPressure(temperature);
            return flashPT(pressure, temperature);
        }

        /**
         *
         * @param pressure
         * @param vaporFraction
         * @return
         */
        const FluidProperties flashPx(double pressure, double vaporFraction) const
        {
            auto temperature = m_equationOfState.saturationTemperature(pressure);
            return flashTx(temperature, vaporFraction);
        }

        /**
         *
         * @param pressure
         * @param enthalpy
         * @return
         */
        const FluidProperties flashPH(double pressure, double enthalpy) const
        {
            using std::get;

            // ===== Define objective function
            auto enthalpyObjFunction = [&](double t) {
                auto phase = FluidProperties(m_equationOfState.computePropertiesPT(pressure, t)).stablePhase().front();
                return idealGasEnthalpy(phase.Temperature) + phase.EnthalpyDeparture - enthalpy;
            };

            // ===== If the fluid is supercritical, compute like so... (single phase)
            if (pressure >= m_pureComponent.property("CriticalPressure")) {
                // ===== Use the (hypothetical) saturation temperature as first guess, and compute until the enthalpy is found.
                return flashPT(pressure, numeric::newton(enthalpyObjFunction, m_equationOfState.saturationTemperature(pressure)));
            }

            // ===== Otherwise, the fluid is sub-critical...
            // ===== First, calculate the saturation properties at the specified pressure.
            auto temperature  = m_equationOfState.saturationTemperature(pressure);
            auto phaseProps = FluidProperties(m_equationOfState.computePropertiesPT(pressure, temperature));

            auto h_v          = idealGasEnthalpy(phaseProps.back().Temperature) + phaseProps.back().EnthalpyDeparture;
            auto h_l          = idealGasEnthalpy(phaseProps.front().Temperature) + phaseProps.front().EnthalpyDeparture;

            // ===== If the specified enthalpy is lower than the saturated liquid enthalpy, the fluid is a compressed liquid.
            if (enthalpy < h_l) {
                return flashPT(pressure, numeric::newton(enthalpyObjFunction, temperature * 0.8));
            }

            // ===== If the specified enthalpy is higher than the saturated vapor entropy, the fluid is superheated vapor.
            if (enthalpy > h_v) {
                return flashPT(pressure, numeric::newton(enthalpyObjFunction, temperature * 1.2));
            }

            // ===== If the fluid is not a compressed liquid nor a superheated vapor, the fluid is two-phase.A
            auto vaporFraction = (h_l - enthalpy) / (h_l - h_v);
            return flashPx(pressure, vaporFraction);
        }

        /**
         *
         * @param pressure
         * @param entropy
         * @return
         */
        const FluidProperties flashPS(double pressure, double entropy) const
        {
            using std::get;

            // ===== Define objective function
            auto entropyObjFunction = [&](double t) {
                auto phase = FluidProperties(m_equationOfState.computePropertiesPT(pressure, t)).stablePhase().front();
                return idealGasEntropy(phase.Temperature, phase.Pressure) + phase.EntropyDeparture - entropy;
            };

            if (pressure >= m_pureComponent.property("CriticalPressure")) {
                return flashPT(pressure, numeric::newton(entropyObjFunction, m_equationOfState.saturationTemperature(pressure)));
            }

            // ===== First, calculate the saturation properties at the specified pressure.
            auto temperature  = m_equationOfState.saturationTemperature(pressure);
            auto phaseProps = FluidProperties(m_equationOfState.computePropertiesPT(pressure, temperature));

            auto s_v          = idealGasEntropy(phaseProps.back().Temperature, phaseProps.back().Pressure) + phaseProps.back().EntropyDeparture;
            auto s_l          = idealGasEntropy(phaseProps.front().Temperature, phaseProps.front().Pressure) + phaseProps.front().EntropyDeparture;

            // ===== If the specified entropy is lower than the saturated liquid entropy, the fluid is a compressed liquid.
            if (entropy < s_l) {
                return flashPT(pressure, numeric::newton(entropyObjFunction, temperature * 0.8));
            }

            // ===== If the specified entropy is higher than the saturated vapor entropy, the fluid is superheated vapor.
            if (entropy > s_v) {
                return flashPT(pressure, numeric::newton(entropyObjFunction, temperature * 1.2));
            }

            // ===== If the fluid is not a compressed liquid nor a superheated vapor, the fluid is two-phase.
            auto vaporFraction = (s_l - entropy) / (s_l - s_v);
            return flashPx(pressure, vaporFraction);
        }

        /**
         *
         * @param temperature
         * @param volume
         * @return
         */
        const FluidProperties flashTV(double temperature, double molarVolume) const
        {
            using std::pow;

            // ===== Fluid is supercritical
            if (temperature >= m_pureComponent.property("CriticalTemperature")) {
                auto phase = FluidProperties(m_equationOfState.computePropertiesTV(temperature, molarVolume)).stablePhase().front();
                return flashPT(phase.Pressure, temperature);
            }

            auto phaseProps = flashTx(temperature, 0.5);

            // ===== Fluid is a compressed liquid or super-heated vapor
            if (molarVolume < phaseProps.front().MolarVolume || molarVolume > phaseProps.back().MolarVolume) {
                auto phase = FluidProperties(m_equationOfState.computePropertiesTV(temperature, molarVolume)).stablePhase().front();
                return flashPT(phase.Pressure, temperature);
            }

            // ===== Fluid is multiphase
            auto vaporFraction = (molarVolume - phaseProps.front().MolarVolume) / (phaseProps.back().MolarVolume - phaseProps.front().MolarVolume);
            return flashTx(temperature, vaporFraction);
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
     * @brief
     * @param specification
     * @param var1
     * @param var2
     * @return
     */
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