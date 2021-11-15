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
        std::optional<double> Pressure;
        std::optional<double> Temperature;
        std::optional<double> MolarVolume;
        std::optional<double> MolarWeight;
        std::optional<double> MolarFlow;
        std::optional<double> Compressibility;
        std::optional<double> FugacityCoefficient;
        std::optional<double> Viscosity;
        std::optional<double> SurfaceTension;
        std::optional<double> ThermalConductivity;
        std::optional<double> Cp;
        std::optional<double> Cv;
        std::optional<double> IsothermalCompressibility;
        std::optional<double> ThermalExpansionCoefficient;
        std::optional<double> JouleThomsonCoefficient;
        std::optional<double> VaporPressure;
        std::optional<double> Enthalpy;
        std::optional<double> Entropy;
        std::optional<double> InternalEnergy;
        std::optional<double> GibbsEnergy;
        std::optional<double> HelmholzEnergy;
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
            return m_data.Pressure.value();
        }

        /**
         * @brief
         * @return
         */
        double getTemperature() const {
            return m_data.Temperature.value();
        }

        /**
         * @brief
         * @return
         */
        double getMolarVolume() const {
            return m_data.MolarVolume.value();
        }

        /**
         * @brief
         * @return
         */
        double getMolarWeight() const {
            return m_data.MolarWeight.value();
        }

        /**
         * @brief
         * @return
         */
        double getMolarFlow() const {
            return m_data.MolarFlow.value();
        }

        /**
         * @brief
         * @return
         */
        double getCompressibility() const {
            return m_data.Compressibility.value();
        }

        /**
         * @brief
         * @return
         */
        double getFugacityCoefficient() const {
            return m_data.FugacityCoefficient.value();
        }

        /**
         * @brief
         * @return
         */
        double getViscosity() const {
            return m_data.Viscosity.value();
        }

        /**
         * @brief
         * @return
         */
        double getSurfaceTension() const {
            return m_data.SurfaceTension.value();
        }

        /**
         * @brief
         * @return
         */
        double getThermalConductivity() const {
            return m_data.ThermalConductivity.value();
        }

        /**
         * @brief
         * @return
         */
        double getCp() const {
            return m_data.Cp.value();
        }

        /**
         * @brief
         * @return
         */
        double getCv() const {
            return m_data.Cv.value();
        }

        /**
         * @brief
         * @return
         */
        double getIsothermalCompressibility() const {
            return m_data.IsothermalCompressibility.value();
        }

        /**
         * @brief
         * @return
         */
        double getThermalExpansionCoefficient() const {
            return m_data.ThermalExpansionCoefficient.value();
        }

        /**
         * @brief
         * @return
         */
        double getJouleThomsonCoefficient() const {
            return m_data.JouleThomsonCoefficient.value();
        }

        /**
         * @brief
         * @return
         */
        double getVaporPressure() const {
            return m_data.VaporPressure.value();
        }

        /**
         * @brief
         * @return
         */
        double getEnthalpy() const {
            return m_data.Enthalpy.value();
        }

        /**
         * @brief
         * @return
         */
        double getEntropy() const {
            return m_data.Entropy.value();
        }

        /**
         * @brief
         * @return
         */
        double getInternalEnergy() const {
            return m_data.InternalEnergy.value();
        }

        /**
         * @brief
         * @return
         */
        double getGibbsEnergy() const {
            return m_data.GibbsEnergy.value();
        }

        /**
         * @brief
         * @return
         */
        double getHelmholzEnergy() const {
            return m_data.HelmholzEnergy.value();
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
