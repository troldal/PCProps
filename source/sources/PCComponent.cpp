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

#include "PCComponent.hpp"

namespace PCProps {

    // ===== Constructor, default
    PCComponent::PCComponent() = default;

    // ===== Constructor, taking a PCComponentData object as an argument
    PCComponent::PCComponent(const PCComponentData& data) : m_data(data) {}

    // ===== Copy constructor
    PCComponent::PCComponent(const PCComponent& other) = default;

    // ===== Move constructor
    PCComponent::PCComponent(PCComponent&& other) noexcept = default;

    // ===== Destructor
    PCComponent::~PCComponent() = default;

    // ===== Copy assignment operator
    PCComponent& PCComponent::operator=(const PCComponent& other) = default;

    // ===== Move assignment operator
    PCComponent& PCComponent::operator=(PCComponent&& other)  noexcept = default;

    // ===== Accessor to the PCComponentData member
    PCComponentData PCComponent::data() const
    {
        return m_data;
    }

    // ===== Get the component name
    const std::string& PCComponent::name() const
    {
        return m_data.name;
    }

    // ===== Get the component formula
    const std::string& PCComponent::formula() const
    {
        return m_data.formula;
    }

    // ===== Get the component CAS registration number
    const std::string& PCComponent::casrn() const
    {
        return m_data.casrn;
    }

    // ===== Get the component SMILES string
    const std::string& PCComponent::smiles() const
    {
        return m_data.smiles;
    }

    // ===== Check if the molecular weight has been set
    bool PCComponent::hasMolecularWeight() const
    {
        return m_data.molecularWeight.has_value();
    }

    // ===== Check if the boiling point temperature has been set
    bool PCComponent::hasBoilingTemperature() const
    {
        return m_data.boilingTemperature.has_value();
    }

    // ===== Check if the freezing temperature has been set
    bool PCComponent::hasFreezingTemperature() const
    {
        return m_data.freezingTemperature.has_value();
    }

    // ===== Check if the critical temperature has been set
    bool PCComponent::hasCriticalTemperature() const
    {
        return m_data.criticalTemperature.has_value();
    }

    // ===== Check if the critical pressure has been set
    bool PCComponent::hasCriticalPressure() const
    {
        return m_data.criticalPressure.has_value();
    }

    // ===== Check if the critical volume has been set
    bool PCComponent::hasCriticalVolume() const
    {
        return m_data.criticalVolume.has_value();
    }

    // ===== Check if the critical density has been set
    bool PCComponent::hasCriticalDensity() const
    {
        return m_data.criticalDensity.has_value();
    }

    // ===== Check if the critical compressibility has been set
    bool PCComponent::hasCriticalCompressibility() const
    {
        return m_data.criticalCompressibility.has_value();
    }

    // ===== Check if the acentric factor has been set
    bool PCComponent::hasAcentricFactor() const
    {
        return m_data.acentricFactor.has_value();
    }

    // ===== Get molecular weight
    double PCComponent::molecularWeight() const
    {
        if (!m_data.molecularWeight) throw PCPropsException("Error: Molecular weight value not set for component.");
        return m_data.molecularWeight.value();
    }

    // ===== Get normal boiling temperature
    double PCComponent::boilingTemperature() const
    {
        if (!m_data.boilingTemperature) throw PCPropsException("Error: Boiling temperature value not set for component.");
        return m_data.boilingTemperature.value();
    }

    // ===== Get freezing temperature
    double PCComponent::freezingTemperature() const
    {
        if (!m_data.freezingTemperature) throw PCPropsException("Error: Freezing temperature value not set for component.");
        return m_data.freezingTemperature.value();
    }

    // ===== Get critical temperature
    double PCComponent::criticalTemperature() const
    {
        if (!m_data.criticalTemperature) throw PCPropsException("Error: Critical temperature value not set for component.");
        return m_data.criticalTemperature.value();
    }

    // ===== Get critical pressure
    double PCComponent::criticalPressure() const
    {
        if (!m_data.criticalPressure) throw PCPropsException("Error: Critical pressure value not set for component.");
        return m_data.criticalPressure.value();
    }

    // ===== Get critical volume
    double PCComponent::criticalVolume() const
    {
        if (!m_data.criticalVolume) throw PCPropsException("Error: Critical volume value not set for component.");
        return m_data.criticalVolume.value();
    }

    // ===== Get critical density
    double PCComponent::criticalDensity() const
    {
        if (!m_data.criticalDensity) throw PCPropsException("Error: Critical density value not set for component.");
        return m_data.criticalDensity.value();
    }

    // ===== Get critical compressibility factor
    double PCComponent::criticalCompressibility() const
    {
        if (!m_data.criticalCompressibility) throw PCPropsException("Error: Critical compressibility factor not set for component.");
        return m_data.criticalCompressibility.value();
    }

    // ===== Get acentric factor
    double PCComponent::acentricFactor() const
    {
        if (!m_data.acentricFactor) throw PCPropsException("Error: Acentric factor not set for component.");
        return m_data.acentricFactor.value();
    }

    // ===== Check if the vapor pressure function has been set
    bool PCComponent::hasVaporPressureFunction() const
    {
        return static_cast<bool>(m_data.vaporPressureFunction);
    }

    // ===== Check if the liquid density function has been set
    bool PCComponent::hasLiquidDensityFunction() const
    {
        return static_cast<bool>(m_data.liquidDensityFunction);
    }

    // ===== Check if the surface tension function has been set
    bool PCComponent::hasSurfaceTensionFunction() const
    {
        return static_cast<bool>(m_data.surfaceTensionFunction);
    }

    // ===== Check if the heat of vaporization function has been set
    bool PCComponent::hasHeatOfVaporizationFunction() const
    {
        return static_cast<bool>(m_data.heatOfVaporizationFunction);
    }

    // ===== Check if the vapor thermal conductivity function has been set
    bool PCComponent::hasVaporThermalConductivityFunction() const
    {
        return static_cast<bool>(m_data.vaporThermalConductivityFunction);
    }

    // ===== Check if the liquid thermal conductivity function has been set
    bool PCComponent::hasLiquidThermalConductivityFunction() const
    {
        return static_cast<bool>(m_data.liquidThermalConductivityFunction);
    }

    // ===== Check if the vapor viscosity function has been set
    bool PCComponent::hasVaporViscosityFunction() const
    {
        return static_cast<bool>(m_data.vaporViscosityFunction);
    }

    // ===== Check if the liquid viscosity function has been set
    bool PCComponent::hasLiquidViscosityFunction() const
    {
        return static_cast<bool>(m_data.liquidViscosityFunction);
    }

    // ===== Check if the ideal gas Cp function has been set
    bool PCComponent::hasIdealGasCpFunction() const
    {
        return static_cast<bool>(m_data.idealGasCpFunction);
    }

    // ===== Check if the liquid Cp function has been set
    bool PCComponent::hasLiquidCpFunction() const
    {
        return static_cast<bool>(m_data.liquidCpFunction);
    }

    // ===== Compute the vapor pressure at the specified temperature
    double PCComponent::vaporPressure(double temperature) const
    {
        if (!m_data.vaporPressureFunction) throw PCPropsException("Error: Vapor pressure function not set for component.");
        return m_data.vaporPressureFunction(temperature);
    }

    // ===== Compute the liquid density at the specified temperature
    double PCComponent::liquidDensity(double temperature) const
    {
        if (!m_data.liquidDensityFunction) throw PCPropsException("Error: Liquid density function not set for component.");
        return m_data.liquidDensityFunction(temperature);
    }

    // ===== Compute the surface tension at the specified temperature
    double PCComponent::surfaceTension(double temperature) const
    {
        if (!m_data.surfaceTensionFunction) throw PCPropsException("Error: Surface tension function not set for component.");
        return m_data.surfaceTensionFunction(temperature);
    }

    // ===== Compute the latent heat at the specified temperature
    double PCComponent::heatOfVaporization(double temperature) const
    {
        if (!m_data.heatOfVaporizationFunction) throw PCPropsException("Error: Latent heat function not set for component.");
        return m_data.heatOfVaporizationFunction(temperature);
    }

    // ===== Compute the vapor thermal conductivity at the specified temperature
    double PCComponent::vaporThermalConductivity(double temperature) const
    {
        if (!m_data.vaporThermalConductivityFunction) throw PCPropsException("Error: Vapor thermal conductivity function not set for component.");
        return m_data.vaporThermalConductivityFunction(temperature);
    }

    // ===== Compute the liquid thermal conductivity at the specified temperature
    double PCComponent::liquidThermalConductivity(double temperature) const
    {
        if (!m_data.liquidThermalConductivityFunction) throw PCPropsException("Error: Liquid thermal conductivity function not set for component.");
        return m_data.liquidThermalConductivityFunction(temperature);
    }

    // ===== Compute the vapor viscosity at the specified temperature
    double PCComponent::vaporViscosity(double temperature) const
    {
        if (!m_data.vaporViscosityFunction) throw PCPropsException("Error: Vapor viscosity function not set for component.");
        return m_data.vaporViscosityFunction(temperature);
    }

    // ===== Compute the liquid viscosity at the specified temperature
    double PCComponent::liquidViscosity(double temperature) const
    {
        if (!m_data.liquidViscosityFunction) throw PCPropsException("Error: Liquid viscosity function not set for component.");
        return m_data.liquidViscosityFunction(temperature);
    }

    // ===== Compute the ideal gas Cp at the specified temperature
    double PCComponent::idealGasCp(double temperature) const
    {
        if (!m_data.idealGasCpFunction) throw PCPropsException("Error: Ideal gas Cp function not set for component.");
        return m_data.idealGasCpFunction(temperature);
    }

    // ===== Compute the liquid Cp at the specified temperature
    double PCComponent::liquidCp(double temperature) const
    {
        if (!m_data.liquidCpFunction) throw PCPropsException("Error: Liquid Cp function not set for component.");
        return m_data.liquidCpFunction(temperature);
    }

} // namespace PCProps