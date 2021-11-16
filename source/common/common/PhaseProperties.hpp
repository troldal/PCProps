//
// Created by Kenneth Balslev on 14/11/2021.
//

#ifndef PCPROPS_PHASEPROPERTIES_HPP
#define PCPROPS_PHASEPROPERTIES_HPP

#include <optional>
#include "PropertyData.hpp"

namespace PCProps
{
    struct PCPhaseProperties
    {
        double Pressure {0.0};
        double Temperature {0.0};
        double MolarVolume {0.0};
        double MolarWeight {0.0};
        double MolarFlow {0.0};
        double Compressibility {0.0};
        double FugacityCoefficient {0.0};
        double Viscosity {0.0};
        double SurfaceTension {0.0};
        double ThermalConductivity {0.0};
        double Cp {0.0};
        double Cv {0.0};
        double IsothermalCompressibility {0.0};
        double ThermalExpansionCoefficient {0.0};
        double JouleThomsonCoefficient {0.0};
        double VaporPressure {0.0};
        double Enthalpy {0.0};
        double Entropy {0.0};
        double InternalEnergy {0.0};
        double GibbsEnergy {0.0};
        double HelmholzEnergy {0.0};
    };

    class PhaseProperties {

        PCPhaseProperties m_data = {};

    public:

        // ===== Constructors & Assignment Operators ===== //

        /**
         * @brief Default constructor.
         */
        PhaseProperties() = default;

        /**
         * @brief Constructor, taking a PCComponentData object as an argument.
         * @param data A PCComponentData object with the component data.
         */
        explicit PhaseProperties(const PCPhaseProperties& data) : m_data(data) {}

        /**
         * @brief Copy constructor.
         */
        PhaseProperties(const PhaseProperties& other) = default;

        /**
         * @brief Move constructor.
         */
        PhaseProperties(PhaseProperties&& other) noexcept = default;

        /**
         * @brief Destructor.
         */
        ~PhaseProperties() = default;

        /**
         * @brief Copy assignment operator.
         */
        PhaseProperties& operator=(const PhaseProperties& other) = default;

        /**
         * @brief Move assignment operator.
         */
        PhaseProperties& operator=(PhaseProperties&& other) noexcept = default;

        /**
         * @brief
         * @return
         */
        double getPressure() const {
            return m_data.Pressure;
        }

        /**
         * @brief
         * @return
         */
        double getTemperature() const {
            return m_data.Temperature;
        }

        /**
         * @brief
         * @return
         */
        double getMolarVolume() const {
            return m_data.MolarVolume;
        }

        /**
         * @brief
         * @return
         */
        double getMolarWeight() const {
            return m_data.MolarWeight;
        }

        /**
         * @brief
         * @return
         */
        double getMolarFlow() const {
            return m_data.MolarFlow;
        }

        /**
         * @brief
         * @return
         */
        double getCompressibility() const {
            return m_data.Compressibility;
        }

        /**
         * @brief
         * @return
         */
        double getFugacityCoefficient() const {
            return m_data.FugacityCoefficient;
        }

        /**
         * @brief
         * @return
         */
        double getViscosity() const {
            return m_data.Viscosity;
        }

        /**
         * @brief
         * @return
         */
        double getSurfaceTension() const {
            return m_data.SurfaceTension;
        }

        /**
         * @brief
         * @return
         */
        double getThermalConductivity() const {
            return m_data.ThermalConductivity;
        }

        /**
         * @brief
         * @return
         */
        double getCp() const {
            return m_data.Cp;
        }

        /**
         * @brief
         * @return
         */
        double getCv() const {
            return m_data.Cv;
        }

        /**
         * @brief
         * @return
         */
        double getIsothermalCompressibility() const {
            return m_data.IsothermalCompressibility;
        }

        /**
         * @brief
         * @return
         */
        double getThermalExpansionCoefficient() const {
            return m_data.ThermalExpansionCoefficient;
        }

        /**
         * @brief
         * @return
         */
        double getJouleThomsonCoefficient() const {
            return m_data.JouleThomsonCoefficient;
        }

        /**
         * @brief
         * @return
         */
        double getVaporPressure() const {
            return m_data.VaporPressure;
        }

        /**
         * @brief
         * @return
         */
        double getEnthalpy() const {
            return m_data.Enthalpy;
        }

        /**
         * @brief
         * @return
         */
        double getEntropy() const {
            return m_data.Entropy;
        }

        /**
         * @brief
         * @return
         */
        double getInternalEnergy() const {
            return m_data.InternalEnergy;
        }

        /**
         * @brief
         * @return
         */
        double getGibbsEnergy() const {
            return m_data.GibbsEnergy;
        }

        /**
         * @brief
         * @return
         */
        double getHelmholzEnergy() const {
            return m_data.HelmholzEnergy;
        }

        PCPhase getPhaseData() const {

            PCPhase data;
            data[PCPressure] = getPressure();
            data[PCTemperature] = getTemperature();
            data[PCMolarVolume] = getMolarVolume();
            data[PCMolarWeight] = getMolarWeight();
            data[PCMolarFlow] = getMolarFlow();
            data[PCCompressibility] = getCompressibility();
            data[PCFugacityCoefficient] = getFugacityCoefficient();
            data[PCViscosity] = getViscosity();
            data[PCSurfaceTension] = getSurfaceTension();
            data[PCThermalConductivity] = getThermalConductivity();
            data[PCHeatCapacityCp] = getCp();
            data[PCHeatCapacityCv] = getCv();
            data[PCIsothermalCompressibility] = getIsothermalCompressibility();
            data[PCThermalExpansionCoefficient] = getThermalExpansionCoefficient();
            data[PCJouleThomsonCoefficient] = getJouleThomsonCoefficient();
            data[PCVaporPressure] = getVaporPressure();
            data[PCEnthalpy] = getEnthalpy();
            data[PCEntropy] = getEntropy();
            data[PCInternalEnergy] = getInternalEnergy();
            data[PCGibbsEnergy] = getGibbsEnergy();
            data[PCHelmholzEnergy] = getHelmholzEnergy();

            return data;
        }

    };

}
#endif    // PCPROPS_PHASEPROPERTIES_HPP
