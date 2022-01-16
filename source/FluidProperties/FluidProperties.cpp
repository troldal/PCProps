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


#include "FluidProperties.hpp"

#include <json/json.hpp>

#include <sstream>

namespace PCProps {

    /**
     * @details
     */
    FluidProperties::FluidProperties() = default;

    /**
     * @details
     */
    FluidProperties::FluidProperties(const std::string& JSONData) {

        using nlohmann::json;

        auto parse_json = [&](const json& data) {

            // ===== Embedded lambda function for converting phase type string to enum type.
            auto AsType = [&](const std::string& type) {
                if (type == "VAPOR") return PhaseType::Vapor;
                if (type == "LIQUID") return PhaseType::Liquid;
                return PhaseType::Undefined;
            };

            PhaseProperties phase;

            for (auto item = data.begin(); item != data.end(); ++item) {

                auto key = item.key();

                if (key == "Type") phase.Type = AsType(item.value());
                else if (key == "Name") phase.Name = item.value();
                else if (key == "CAS") phase.CAS = item.value();

                else if (key == "NormalFreezingPoint") phase.NormalFreezingPoint = item.value();
                else if (key == "NormalBoilingPoint") phase.NormalBoilingPoint = item.value();
                else if (key == "CriticalTemperature") phase.CriticalTemperature = item.value();
                else if (key == "CriticalPressure") phase.CriticalPressure = item.value();

                else if (key == "Pressure") phase.Pressure = item.value();
                else if (key == "Temperature") phase.Temperature = item.value();
                else if (key == "MolarVolume") phase.MolarVolume = item.value();
                else if (key == "MolarWeight") phase.MolarWeight = item.value();
                else if (key == "MolarFlow") phase.MolarFlow = item.value();
                else if (key == "Compressibility") phase.Compressibility = item.value();
                else if (key == "FugacityCoefficient") phase.FugacityCoefficient = item.value();
                else if (key == "Viscosity") phase.Viscosity = item.value();
                else if (key == "SurfaceTension") phase.SurfaceTension = item.value();
                else if (key == "ThermalConductivity") phase.ThermalConductivity = item.value();
                else if (key == "CpDeparture") phase.CpDeparture = item.value();
                else if (key == "CvDeparture") phase.CvDeparture = item.value();
                else if (key == "EnthalpyDeparture") phase.EnthalpyDeparture = item.value();
                else if (key == "EntropyDeparture") phase.EntropyDeparture = item.value();
                else if (key == "InternalEnergyDeparture") phase.InternalEnergyDeparture = item.value();
                else if (key == "GibbsEnergyDeparture") phase.GibbsEnergyDeparture = item.value();
                else if (key == "HelmholzEnergyDeparture") phase.HelmholzEnergyDeparture = item.value();

                else if (key == "DPDV") phase.DPDV = item.value();
                else if (key == "DPDT") phase.DPDT = item.value();
                else if (key == "DVDP") phase.DVDP = item.value();
                else if (key == "DVDT") phase.DVDT = item.value();
                else if (key == "DTDV") phase.DTDV = item.value();
                else if (key == "DTDV") phase.DTDV = item.value();
                else if (key == "DTDP") phase.DTDP = item.value();

                else if (key == "Cp") phase.Cp = item.value();
                else if (key == "Cv") phase.Cv = item.value();
                else if (key == "IsothermalCompressibility") phase.IsothermalCompressibility = item.value();
                else if (key == "ThermalExpansionCoefficient") phase.ThermalExpansionCoefficient = item.value();
                else if (key == "JouleThomsonCoefficient") phase.JouleThomsonCoefficient = item.value();
                else if (key == "SpeedOfSound") phase.SpeedOfSound = item.value();
                else if (key == "SaturationPressure") phase.SaturationPressure = item.value();
                else if (key == "SaturationVolume") phase.SaturationVolume = item.value();
                else if (key == "Enthalpy") phase.Enthalpy = item.value();
                else if (key == "Entropy") phase.Entropy = item.value();
                else if (key == "InternalEnergy") phase.InternalEnergy = item.value();
                else if (key == "GibbsEnergy") phase.GibbsEnergy = item.value();
                else if (key == "HelmholzEnergy") phase.HelmholzEnergy = item.value();
            }
            return phase;
        };

        auto fluid = nlohmann::json::parse(JSONData);
        for (const auto& phase : fluid) m_phases.emplace_back(parse_json(phase));
    }

    /**
     * @details
     */
    FluidProperties::FluidProperties(const std::vector<PhaseProperties>& fluidProps) : m_phases(fluidProps) {}

    /**
     * @details
     */
    FluidProperties::FluidProperties(const FluidProperties& other) = default;

    /**
     * @details
     */
    FluidProperties::FluidProperties(FluidProperties&& other) noexcept = default;

    /**
     * @details
     */
    FluidProperties::~FluidProperties() = default;

    /**
     * @details
     */
    FluidProperties& FluidProperties::operator=(const FluidProperties& other) = default;

    /**
     * @details
     */
    FluidProperties& FluidProperties::operator=(FluidProperties&& other) noexcept = default;

    /**
     * @details
     */
    FluidProperties& FluidProperties::operator=(const std::vector<PhaseProperties>& fluidProps)
    {
        m_phases = fluidProps;
        return *this;
    }

    /**
     * @details
     */
    const PhaseProperties& FluidProperties::operator[](int index) const
    {
        return m_phases[index];
    }

    /**
     * @details
     */
    const std::vector<PhaseProperties>& FluidProperties::phases() const
    {
        return m_phases;
    }

    /**
     * @details
     */
    size_t FluidProperties::size() const
    {
        return m_phases.size();
    }

    /**
     * @details
     */
    void FluidProperties::erase(int index) {
        auto iter = m_phases.begin();
        std::advance(iter, index);
        m_phases.erase(iter);
    }

    /**
     * @details
     */
    FluidProperties FluidProperties::stablePhase() const
    {
        return FluidProperties(std::vector {*std::min_element(
            m_phases.begin(),
            m_phases.end(),
            [](const PhaseProperties& a, const PhaseProperties& b){return a.FugacityCoefficient < b.FugacityCoefficient;})});
    }

    /**
     * @details
     */
    FluidProperties FluidProperties::heavyPhase() const
    {
        return FluidProperties(std::vector {*std::min_element(
            m_phases.begin(),
            m_phases.end(),
            [](const PhaseProperties& a, const PhaseProperties& b){return a.MolarVolume < b.MolarVolume;})});
    }

    /**
     * @details
     */
    FluidProperties FluidProperties::lightPhase() const
    {
        return FluidProperties(std::vector {*std::max_element(
            m_phases.begin(),
            m_phases.end(),
            [](const PhaseProperties& a, const PhaseProperties& b){return a.MolarVolume < b.MolarVolume;})});
    }

    /**
     * @details
     */
    JSONString FluidProperties::asJSON() const
    {

        auto make_json = [&](const PhaseProperties& props) {
                    nlohmann::json data;

                    auto TypeAsString = [&](const PhaseType type) {
                        if (type == PhaseType::Vapor) return "VAPOR";
                        if (type == PhaseType::Liquid) return "LIQUID";
                        return "UNDEFINED";
                    };

                    data["Type"] = TypeAsString(props.Type);
                    if (!props.Name.empty()) data["Name"] = props.Name;
                    if (!props.CAS.empty()) data["CAS"] = props.CAS;
                    if (props.NormalFreezingPoint != 0.0) data["NormalFreezingPoint"] = props.NormalFreezingPoint;
                    if (props.NormalBoilingPoint != 0.0) data["NormalBoilingPoint"] = props.NormalBoilingPoint;
                    if (props.CriticalTemperature != 0.0) data["CriticalTemperature"] = props.CriticalTemperature;
                    if (props.CriticalPressure != 0.0) data["CriticalPressure"] = props.CriticalPressure;
                    if (props.Pressure != 0.0) data["Pressure"] = props.Pressure;
                    if (props.Temperature != 0.0) data["Temperature"] = props.Temperature;
                    if (props.MolarVolume != 0.0) data["MolarVolume"] = props.MolarVolume;
                    if (props.MolarWeight != 0.0) data["MolarWeight"] = props.MolarWeight;
                    if (props.MolarFlow != 0.0) data["MolarFlow"] = props.MolarFlow;
                    if (props.Compressibility != 0.0) data["Compressibility"] = props.Compressibility;
                    if (props.FugacityCoefficient != 0.0) data["FugacityCoefficient"] = props.FugacityCoefficient;
                    if (props.Viscosity != 0.0) data["Viscosity"] = props.Viscosity;
                    if (props.SurfaceTension != 0.0) data["SurfaceTension"] = props.SurfaceTension;
                    if (props.ThermalConductivity != 0.0) data["ThermalConductivity"] = props.ThermalConductivity;
                    if (props.CpDeparture != 0.0) data["CpDeparture"] = props.CpDeparture;
                    if (props.CvDeparture != 0.0) data["CvDeparture"] = props.CvDeparture;
                    if (props.EnthalpyDeparture != 0.0) data["EnthalpyDeparture"] = props.EnthalpyDeparture;
                    if (props.EntropyDeparture != 0.0) data["EntropyDeparture"] = props.EntropyDeparture;
                    if (props.InternalEnergyDeparture != 0.0) data["InternalEnergyDeparture"] = props.InternalEnergyDeparture;
                    if (props.GibbsEnergyDeparture != 0.0) data["GibbsEnergyDeparture"] = props.GibbsEnergyDeparture;
                    if (props.HelmholzEnergyDeparture != 0.0) data["HelmholzEnergyDeparture"] = props.HelmholzEnergyDeparture;
                    if (props.DPDV != 0.0) data["DPDV"] = props.DPDV;
                    if (props.DPDT != 0.0) data["DPDT"] = props.DPDT;
                    if (props.DVDP != 0.0) data["DVDP"] = props.DVDP;
                    if (props.DVDT != 0.0) data["DVDT"] = props.DVDT;
                    if (props.DTDV != 0.0) data["DTDV"] = props.DTDV;
                    if (props.DTDP != 0.0) data["DTDP"] = props.DTDP;
                    if (props.Cp != 0.0) data["Cp"] = props.Cp;
                    if (props.Cv != 0.0) data["Cv"] = props.Cv;
                    if (props.IsothermalCompressibility != 0.0) data["IsothermalCompressibility"] = props.IsothermalCompressibility;
                    if (props.ThermalExpansionCoefficient != 0.0) data["ThermalExpansionCoefficient"] = props.ThermalExpansionCoefficient;
                    if (props.JouleThomsonCoefficient != 0.0) data["JouleThomsonCoefficient"] = props.JouleThomsonCoefficient;
                    if (props.SpeedOfSound != 0.0) data["SpeedOfSound"] = props.SpeedOfSound;
                    if (props.SaturationPressure != 0.0) data["SaturationPressure"] = props.SaturationPressure;
                    if (props.SaturationVolume != 0.0) data["SaturationVolume"] = props.SaturationVolume;
                    if (props.Enthalpy != 0.0) data["Enthalpy"] = props.Enthalpy;
                    if (props.Entropy != 0.0) data["Entropy"] = props.Entropy;
                    if (props.InternalEnergy != 0.0) data["InternalEnergy"] = props.InternalEnergy;
                    if (props.GibbsEnergy != 0.0) data["GibbsEnergy"] = props.GibbsEnergy;
                    if (props.HelmholzEnergy != 0.0) data["HelmholzEnergy"] = props.HelmholzEnergy;

                    return data;
        };

        std::vector<nlohmann::json> data;
        for (const auto& phase : m_phases) data.emplace_back(make_json(phase));
        return nlohmann::json(data).dump();
    }

    /**
     * @details
     */
    std::vector<PhaseProperties>::iterator FluidProperties::begin()
    {
        return m_phases.begin();
    }

    /**
     * @details
     */
    std::vector<PhaseProperties>::iterator FluidProperties::end()
    {
        return m_phases.end();
    }

    /**
     * @details
     */
    PhaseProperties& FluidProperties::back()
    {
        return m_phases.back();
    }

    /**
     * @details
     */
    const PhaseProperties& FluidProperties::back() const
    {
        return m_phases.back();
    }

    /**
     * @details
     */
    PhaseProperties& FluidProperties::front()
    {
        return m_phases.front();
    }

    /**
     * @details
     */
    const PhaseProperties& FluidProperties::front() const
    {
        return m_phases.front();
    }

    /**
     * @details
     */
    void FluidProperties::print(std::ostream& stream) {
        stream << std::setprecision(6) << std::fixed;
        auto width = 25;

        auto TypeAsString = [&](const PhaseType type) {
            if (type == PhaseType::Vapor) return "VAPOR";
            if (type == PhaseType::Liquid) return "LIQUID";
            return "UNDEFINED";
        };

        stream << "Type                                : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << TypeAsString(phase.Type);
        stream << std::endl;

        stream << "Molar Flow                      [-] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.MolarFlow;
        stream << std::endl;

        stream << "Molar Volume                [kg/m3] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << 1 / (phase.MolarVolume / phase.MolarWeight * 1000);
        stream << std::endl;

        stream << "Surface Tension               [N/m] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.SurfaceTension;
        stream << std::endl;

        stream << "Thermal Conductivity        [W/m-K] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.ThermalConductivity;
        stream << std::endl;

        stream << "Viscosity                    [Pa-s] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.Viscosity;
        stream << std::endl;

        stream << "CpDeparture                [J/kg-K] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.CpDeparture / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "CvDeparture                [J/kg-K] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.CvDeparture / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "EnthalpyDeparture            [J/kg] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.EnthalpyDeparture / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "EntropyDeparture           [J/kg-K] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.EntropyDeparture / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "InternalEnergyDeparture      [J/kg] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.InternalEnergyDeparture / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "GibbsEnergyDeparture         [J/kg] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.GibbsEnergyDeparture / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "HelmholzEnergyDeparture      [J/kg] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.HelmholzEnergyDeparture / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "DPDV                        [Pa/m3] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.DPDV;
        stream << std::endl;

        stream << "DPDT                         [Pa/K] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.DPDT;
        stream << std::endl;

        stream << "DVDP                        [m3/Pa] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.DVDP;
        stream << std::endl;

        stream << "DVDT                         [m3/K] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.DVDT;
        stream << std::endl;

        stream << "DTDV                         [K/m3] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.DTDV;
        stream << std::endl;

        stream << "DTDP                         [K/Pa] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.DTDP;
        stream << std::endl;

        stream << "Cp                         [J/kg-K] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.Cp / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "Cv                         [J/kg-K] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.Cv / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "Isothermal Compressibility   [1/Pa] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.IsothermalCompressibility;
        stream << std::endl;

        stream << "Thermal Expansion Coefficient [1/K] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.ThermalExpansionCoefficient;
        stream << std::endl;

        stream << "Joule-Thomson Coefficient    [K/Pa] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.JouleThomsonCoefficient;
        stream << std::endl;

        stream << "Speed of Sound                [m/s] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.SpeedOfSound;
        stream << std::endl;

        stream << "Molecular Weight            [g/mol] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.MolarWeight;
        stream << std::endl;

        stream << "Temperature                     [K] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.Temperature;
        stream << std::endl;

        stream << "Pressure                       [Pa] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.Pressure;
        stream << std::endl;

        stream << "Compressibility                 [-] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.Compressibility;
        stream << std::endl;

        stream << "Fugacity Coefficient            [-] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.FugacityCoefficient;
        stream << std::endl;

        stream << "Saturation Pressure            [Pa] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.SaturationPressure;
        stream << std::endl;

        stream << "Saturation Volume           [kg/m3] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << 1 / (phase.SaturationVolume / phase.MolarWeight * 1000);
        stream << std::endl;

        stream << "Enthalpy                     [J/kg] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.Enthalpy / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "Entropy                    [J/kg-K] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.Entropy / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "Internal Energy              [J/kg] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.InternalEnergy / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "Gibbs Energy                 [J/kg] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.GibbsEnergy / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "Helmholz Energy              [J/kg] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.HelmholzEnergy / phase.MolarWeight * 1000;
        stream << std::endl;
    }
}