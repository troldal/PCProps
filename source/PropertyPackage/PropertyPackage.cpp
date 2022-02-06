//
// Created by Kenneth Balslev on 18/01/2021.
//

#include "PropertyPackage.hpp"

#include <FluidProperties.hpp>
#include <numerics.hpp>
#include <common/Globals.hpp>

#include <stdexcept>
#include <tuple>
#include <cmath>

using JSONString = std::string;

namespace PCProps
{
    class PropertyPackage::impl
    {
        // ===== Objects representing the pure component and equation of state.
        IPureComponent          m_pureComponent {};
        IEquationOfState        m_equationOfState {};

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
         * @param liquid
         * @return
         */
        PhaseProperties& computeLiquidProperties(PhaseProperties& liquid) const
        {
            // ===== Compute the compressed molar volume.
            liquid.MolarVolume = m_pureComponent.correlation("CompressedLiquidVolume", { liquid.Temperature, liquid.Pressure, liquid.SaturationPressure, liquid.SaturationVolume });

            // ===== Compute the surface tension.
            liquid.SurfaceTension = 0.0;

            // ===== Compute the thermal conductivity (low pressure)
            liquid.ThermalConductivity = m_pureComponent.correlation("LiquidThermalConductivity", liquid.Temperature);

            // ===== Compute the compressed liquid viscosity
            liquid.Viscosity = m_pureComponent.correlation(
                "CompressedLiquidViscosity",
                { liquid.Temperature, liquid.Pressure, liquid.SaturationPressure, m_pureComponent.correlation("SaturatedLiquidViscosity", liquid.Temperature) });

            return liquid;
        }

        /**
         * @brief Compute the vapor properties for a phase.
         * @param vapor The (vapor) phase for which to calcResults properties.
         * @return The input phase updated with vapor properties.
         */
        PhaseProperties& computeVaporProperties(PhaseProperties& vapor) const
        {
            // ===== Compute the surface tension (not applicable for vapors; set to zero)
            vapor.SurfaceTension      = 0.0;

            // ===== Compute the thermal conductivity (low pressure)
            vapor.ThermalConductivity = m_pureComponent.correlation("VaporThermalConductivity", vapor.Temperature);

            // ===== Compute the compressed liquid viscosity
            vapor.Viscosity = m_pureComponent.correlation(
                "CompressedVaporViscosity",
                { vapor.Temperature, vapor.Pressure, vapor.SaturationPressure, m_pureComponent.correlation("SaturatedVaporViscosity", vapor.Temperature) });

            return vapor;
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

                if (phase.Type == PhaseType::Liquid) computeLiquidProperties(phase);
                if (phase.Type == PhaseType::Vapor) computeVaporProperties(phase);

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
            auto result =  computeProperties(phaseProps).stablePhase();
            result.front().MolarFlow = 1.0;
            result.front().MolarFraction = 1.0;
            return result;
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
                    phaseProps.back().MolarFraction = 1.0;
                    return phaseProps;
                }

                // ===== If the specified vapor fraction is 0.0 (or lower), the fluid is a saturated liquid.
                if (vaporFraction <= 0.0) {
                    phaseProps = phaseProps.heavyPhase();
                    phaseProps.front().MolarFlow = 1.0;
                    phaseProps.front().MolarFraction = 1.0;
                    return phaseProps;
                }

                // ===== If the vapor fraction is between 0.0 and 1.0, the fluid is two-phase.

                if (phaseProps.size() == 1) {
                    phaseProps.front().MolarFlow = 1.0;
                    phaseProps.front().MolarFraction = 1.0;
                    return phaseProps;
                }

                else {
                    phaseProps.front().MolarFlow = 1.0 - vaporFraction;
                    phaseProps.front().MolarFraction = 1.0 - vaporFraction;
                    phaseProps.back().MolarFlow = vaporFraction;
                    phaseProps.back().MolarFraction = vaporFraction;
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

            // ===== If the fluid is supercritical, calcResults like so... (single phase)
            if (pressure >= m_pureComponent.property("CriticalPressure")) {
                // ===== Use the (hypothetical) saturation temperature as first guess, and calcResults until the enthalpy is found.
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
                return flashPT(pressure, numeric::newton(enthalpyObjFunction, temperature * (1.0 - 1E-6)));
            }

            // ===== If the specified enthalpy is higher than the saturated vapor entropy, the fluid is superheated vapor.
            if (enthalpy > h_v) {
                return flashPT(pressure, numeric::newton(enthalpyObjFunction, temperature * (1.0 + 1E-6)));
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

            // ===== PropertyPackage is supercritical
            if (temperature >= m_pureComponent.property("CriticalTemperature")) {
                auto phase = FluidProperties(m_equationOfState.computePropertiesTV(temperature, molarVolume)).stablePhase().front();
                return flashPT(phase.Pressure, temperature);
            }

            auto phaseProps = flashTx(temperature, 0.5);

            // ===== PropertyPackage is a compressed liquid or super-heated vapor
            if (molarVolume < phaseProps.front().MolarVolume || molarVolume > phaseProps.back().MolarVolume) {
                auto phase = FluidProperties(m_equationOfState.computePropertiesTV(temperature, molarVolume)).stablePhase().front();
                return flashPT(phase.Pressure, temperature);
            }

            // ===== PropertyPackage is multiphase
            auto vaporFraction = (molarVolume - phaseProps.front().MolarVolume) / (phaseProps.back().MolarVolume - phaseProps.front().MolarVolume);
            return flashTx(temperature, vaporFraction);
        }
    };

    // =====================================================================
    // PUBLIC INTERFACE
    // =====================================================================

    /**
     * Default constructor
     */
    PropertyPackage::PropertyPackage() = default;

    /**
     * Constructor
     */
    PropertyPackage::PropertyPackage(const IPureComponent& pureComponent, const IEquationOfState& eos) : m_impl(std::make_unique<impl>(pureComponent, eos)) {}

    /**
     * Copy constructor
     */
    PropertyPackage::PropertyPackage(const PropertyPackage& other) : m_impl(std::make_unique<impl>(*other.m_impl)) {}

    /**
     * Move constructor
     */
    PropertyPackage::PropertyPackage(PropertyPackage&& other) noexcept = default;

    /**
     * Destructor
     */
    PropertyPackage::~PropertyPackage() = default;

    /**
     * Copy assignment operator
     */
    PropertyPackage& PropertyPackage::operator=(const PropertyPackage& other)
    {
        PropertyPackage copy = other;
        *this      = std::move(copy);
        return *this;
    };

    /**
     * Move assignment operator
     */
    PropertyPackage& PropertyPackage::operator=(PropertyPackage&& other) noexcept = default;

    /**
     * @brief
     * @param specification
     * @param var1
     * @param var2
     * @return
     */
    JSONString PropertyPackage::flash(const std::string& specification, double var1, double var2) const
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

}    // namespace PCProps