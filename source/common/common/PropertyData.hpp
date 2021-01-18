/*

8888888b.   .d8888b.  8888888b.
888   Y88b d88P  Y88b 888   Y88b
888    888 888    888 888    888
888   d88P 888        888   d88P 888d888 .d88b.  88888b.  .d8888b
8888888P"  888        8888888P"  888P"  d88""88b 888 "88b 88K
888        888    888 888        888    888  888 888  888 "Y8888b.
888        Y88b  d88P 888        888    Y88..88P 888 d88P      X88
888         "Y8888P"  888        888     "Y88P"  88888P"   88888P'
                                                 888
                                                 888
                                                 888

Copyright (c) 2020 Kenneth Troldal Balslev

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#ifndef PCPROPS_PROPERTYDATA_HPP
#define PCPROPS_PROPERTYDATA_HPP

#include <array>
#include <iomanip>
#include <iostream>
#include <vector>

namespace PCProps
{
    using PCPhaseData = std::array<double, 20>;
    using PCPhases    = std::vector<PCPhaseData>;

    enum PCPhaseDataElement {
        PCPressure,
        PCTemperature,
        PCMolarVolume,
        PCMolarWeight,
        PCMolarFraction,
        PCCompressibility,
        PCFugacityCoefficient,
        PCViscosity,
        PCSurfaceTension,
        PCThermalConductivity,
        PCHeatCapacityCp,
        PCHeatCapacityCv,
        PCIsothermalCompressibility,
        PCThermalExpansionCoefficient,
        PCJouleThomsonCoefficient,
        PCEnthalpy,
        PCEntropy,
        PCInternalEnergy,
        PCGibbsEnergy,
        PCHelmholzEnergy
    };

    /**
     * @brief
     */
    class PCPhase
    {
        PCPhaseData m_data;

    public:
        PCPhase() : m_data() {}

        explicit PCPhase(const PCPhaseData& data) : m_data(data) {}

        inline double pressure() const
        {
            return m_data[PCPressure];
        }

        inline void setPressure(double pressure)
        {
            m_data[PCPressure] = pressure;
        }

        inline double temperature() const
        {
            return m_data[PCTemperature];
        }

        inline void setTemperature(double temperature)
        {
            m_data[PCTemperature] = temperature;
        }

        inline double molarVolume() const
        {
            return m_data[PCMolarVolume];
        }

        inline void setMolarVolume(double molarVolume)
        {
            m_data[PCMolarVolume] = molarVolume;
        }

        inline double molarWeight() const
        {
            return m_data[PCMolarWeight];
        }

        inline void setMolarWeight(double molarWeight)
        {
            m_data[PCMolarWeight] = molarWeight;
        }

        inline double molarFraction() const
        {
            return m_data[PCMolarFraction];
        }

        inline void setMolarFraction(double molarFraction)
        {
            m_data[PCMolarFraction] = molarFraction;
        }

        inline double compressibility() const
        {
            return m_data[PCCompressibility];
        }

        inline void setCompressibility(double compressibility)
        {
            m_data[PCCompressibility] = compressibility;
        }

        inline double fugacityCoefficient() const
        {
            return m_data[PCFugacityCoefficient];
        }

        inline void setFugacityCoefficient(double fugacityCoefficient)
        {
            m_data[PCFugacityCoefficient] = fugacityCoefficient;
        }

        inline double viscosity() const
        {
            return m_data[PCViscosity];
        }

        inline void setViscosity(double viscosity)
        {
            m_data[PCViscosity] = viscosity;
        }

        inline double surfaceTension() const
        {
            return m_data[PCSurfaceTension];
        }

        inline void setSurfaceTension(double surfaceTension)
        {
            m_data[PCSurfaceTension] = surfaceTension;
        }

        inline double thermalConductivity() const
        {
            return m_data[PCThermalConductivity];
        }

        inline void setThermalConductivity(double thermalConductivity)
        {
            m_data[PCThermalConductivity] = thermalConductivity;
        }

        inline double heatCapacityCp() const
        {
            return m_data[PCHeatCapacityCp];
        }

        inline void setHeatCapacityCp(double heatCapacity)
        {
            m_data[PCHeatCapacityCp] = heatCapacity;
        }

        inline double heatCapacityCv() const
        {
            return m_data[PCHeatCapacityCv];
        }

        inline void setHeatCapacityCv(double heatCapacity)
        {
            m_data[PCHeatCapacityCv] = heatCapacity;
        }

        inline double isothermalCompressibility() const {
            return m_data[PCIsothermalCompressibility];
        }

        inline void setIsothermalCompressibility(double isothermalCompressibility) {
            m_data[PCIsothermalCompressibility] = isothermalCompressibility;
        }

        inline double thermalExpansionCoefficient() const {
            return m_data[PCThermalExpansionCoefficient];
        }

        inline void setThermalExpansionCoefficient(double thermalExpansionCoefficient) {
            m_data[PCThermalExpansionCoefficient] = thermalExpansionCoefficient;
        }

        inline double jouleThomsonCoefficient() const {
            return m_data[PCJouleThomsonCoefficient];
        }

        inline void setJouleThomsonCoefficient(double jouleThomsonCoefficient) {
            m_data[PCJouleThomsonCoefficient] = jouleThomsonCoefficient;
        }

        inline double enthalpy() const
        {
            return m_data[PCEnthalpy];
        }

        inline void setEnthalpy(double enthalpy)
        {
            m_data[PCEnthalpy] = enthalpy;
        }

        inline double entropy() const
        {
            return m_data[PCEntropy];
        }

        inline void setEntropy(double entropy)
        {
            m_data[PCEntropy] = entropy;
        }

        inline double internalEnergy() const
        {
            return m_data[PCInternalEnergy];
        }

        inline void setInternalEnergy(double internalEnergy)
        {
            m_data[PCInternalEnergy] = internalEnergy;
        }

        inline double gibbsEnergy() const
        {
            return m_data[PCGibbsEnergy];
        }

        inline void setGibbsEnergy(double gibbsEnergy)
        {
            m_data[PCGibbsEnergy] = gibbsEnergy;
        }

        inline double helmholzEnergy() const
        {
            return m_data[PCHelmholzEnergy];
        }

        inline void setHelmholzEnergy(double helmholzEnergy)
        {
            m_data[PCHelmholzEnergy] = helmholzEnergy;
        }

        inline const PCPhaseData& data() const
        {
            return m_data;
        }
    };

    /**
     * @brief
     * @param stream
     * @param properties
     * @return
     */
    inline std::ostream& operator<<(std::ostream& stream, const PCProps::PCPhase& properties)
    {
        return stream << std::setprecision(8) << std::fixed
                      << "Molar Fraction                : " << std::right << std::setw(20) << properties.molarFraction() << std::endl
                      << "Molar Volume                  : " << std::right << std::setw(20) << properties.molarVolume() << " m3/mol" << std::endl
                      << "Surface Tension               : " << std::right << std::setw(20) << properties.surfaceTension() << " N/m" << std::endl
                      << "Thermal Conductivity          : " << std::right << std::setw(20) << properties.thermalConductivity() << " W/m-K" << std::endl
                      << "Viscosity                     : " << std::right << std::setw(20) << properties.viscosity() << " Pa-s" << std::endl
                      << "Cp                            : " << std::right << std::setw(20) << properties.heatCapacityCp() << " J/mol-K" << std::endl
                      << "Cv                            : " << std::right << std::setw(20) << properties.heatCapacityCv() << " J/mol-K" << std::endl
                      << "Isothermal Compressibility    : " << std::right << std::setw(20) << properties.isothermalCompressibility() << " 1/Pa" << std::endl
                      << "Thermal Expansion Coefficient : " << std::right << std::setw(20) << properties.thermalExpansionCoefficient() << " 1/K" << std::endl
                      << "Joule-Thomson Coefficient     : " << std::right << std::setw(20) << properties.jouleThomsonCoefficient() << " K/Pa" << std::endl
                      << "Molecular Weight              : " << std::right << std::setw(20) << properties.molarWeight() << " g/mol" << std::endl
                      << "Temperature                   : " << std::right << std::setw(20) << properties.temperature() << " K" << std::endl
                      << "Pressure                      : " << std::right << std::setw(20) << properties.pressure() << " Pa" << std::endl
                      << "Compressibility               : " << std::right << std::setw(20) << properties.compressibility() << std::endl
                      << "Fugacity Coefficient          : " << std::right << std::setw(20) << properties.fugacityCoefficient() << std::endl
                      << "Enthalpy                      : " << std::right << std::setw(20) << properties.enthalpy() << " J/mol" << std::endl
                      << "Entropy                       : " << std::right << std::setw(20) << properties.entropy() << " J/mol-K" << std::endl
                      << "Internal Energy               : " << std::right << std::setw(20) << properties.internalEnergy() << " J/mol" << std::endl
                      << "Gibbs Energy                  : " << std::right << std::setw(20) << properties.gibbsEnergy() << " J/mol" << std::endl
                      << "Helmholz Energy               : " << std::right << std::setw(20) << properties.helmholzEnergy() << " J/mol" << std::endl;
    }

    inline std::ostream& operator<<(std::ostream& stream, const PCProps::PCPhaseData& properties)
    {
        return stream << PCPhase(properties);
    }

}    // namespace PCProps

#endif    // PCPROPS_PROPERTYDATA_HPP
