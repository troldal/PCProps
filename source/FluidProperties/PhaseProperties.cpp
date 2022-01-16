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

#include "PhaseProperties.hpp"

#include <json/json.hpp>

namespace PCProps
{

    /**
     * @details
     */
    PhaseProperties::PhaseProperties() = default;

    /**
     * @details
     */
    PhaseProperties::PhaseProperties(const std::string& JSONData)
    {
        auto data = nlohmann::json::parse(JSONData);

        // ===== Embedded lambda function for converting phase type string to enum type.
        auto AsType = [&](const std::string& type) {
            if (type == "VAPOR") return PhaseType::Vapor;
            if (type == "LIQUID") return PhaseType::Liquid;
            return PhaseType::Undefined;
        };

        for (auto item = data.begin(); item != data.end(); ++item) {

            auto key = item.key();

            if (key == "Type") Type = AsType(item.value());
            else if (key == "Name") Name = item.value();
            else if (key == "CAS") CAS = item.value();

            else if (key == "NormalFreezingPoint") NormalFreezingPoint = item.value();
            else if (key == "NormalBoilingPoint") NormalBoilingPoint = item.value();
            else if (key == "CriticalTemperature") CriticalTemperature = item.value();
            else if (key == "CriticalPressure") CriticalPressure = item.value();

            else if (key == "Pressure") Pressure = item.value();
            else if (key == "Temperature") Temperature = item.value();
            else if (key == "MolarVolume") MolarVolume = item.value();
            else if (key == "MolarWeight") MolarWeight = item.value();
            else if (key == "MolarFlow") MolarFlow = item.value();
            else if (key == "Compressibility") Compressibility = item.value();
            else if (key == "FugacityCoefficient") FugacityCoefficient = item.value();
            else if (key == "Viscosity") Viscosity = item.value();
            else if (key == "SurfaceTension") SurfaceTension = item.value();
            else if (key == "ThermalConductivity") ThermalConductivity = item.value();
            else if (key == "CpDeparture") CpDeparture = item.value();
            else if (key == "CvDeparture") CvDeparture = item.value();
            else if (key == "EnthalpyDeparture") EnthalpyDeparture = item.value();
            else if (key == "EntropyDeparture") EntropyDeparture = item.value();
            else if (key == "InternalEnergyDeparture") InternalEnergyDeparture = item.value();
            else if (key == "GibbsEnergyDeparture") GibbsEnergyDeparture = item.value();
            else if (key == "HelmholzEnergyDeparture") HelmholzEnergyDeparture = item.value();

            else if (key == "DPDV") DPDV = item.value();
            else if (key == "DPDT") DPDT = item.value();
            else if (key == "DVDP") DVDP = item.value();
            else if (key == "DVDT") DVDT = item.value();
            else if (key == "DTDV") DTDV = item.value();
            else if (key == "DTDV") DTDV = item.value();
            else if (key == "DTDP") DTDP = item.value();

            else if (key == "Cp") Cp = item.value();
            else if (key == "Cv") Cv = item.value();
            else if (key == "IsothermalCompressibility") IsothermalCompressibility = item.value();
            else if (key == "ThermalExpansionCoefficient") ThermalExpansionCoefficient = item.value();
            else if (key == "JouleThomsonCoefficient") JouleThomsonCoefficient = item.value();
            else if (key == "SpeedOfSound") SpeedOfSound = item.value();
            else if (key == "SaturationPressure") SaturationPressure = item.value();
            else if (key == "SaturationVolume") SaturationVolume = item.value();
            else if (key == "Enthalpy") Enthalpy = item.value();
            else if (key == "Entropy") Entropy = item.value();
            else if (key == "InternalEnergy") InternalEnergy = item.value();
            else if (key == "GibbsEnergy") GibbsEnergy = item.value();
            else if (key == "HelmholzEnergy") HelmholzEnergy = item.value();
        }
    }

    /**
     * @details
     */
    PhaseProperties::PhaseProperties(const PhaseProperties& other) = default;

    /**
     * @details
     */
    PhaseProperties::PhaseProperties(PhaseProperties&& other) noexcept = default;

    /**
     * @details
     */
    PhaseProperties::~PhaseProperties() = default;

    /**
     * @details
     */
    PhaseProperties& PhaseProperties::operator=(const PhaseProperties& other) = default;

    /**
     * @details
     */
    PhaseProperties& PhaseProperties::operator=(PhaseProperties&& other) noexcept = default;

    /**
     * @details
     */
    JSONString PhaseProperties::asJSON() const
    {
        nlohmann::json data;

        auto TypeAsString = [&](const PhaseType type) {
            if (type == PhaseType::Vapor) return "VAPOR";
            if (type == PhaseType::Liquid) return "LIQUID";
            return "UNDEFINED";
        };

        data["Type"] = TypeAsString(Type);
        if (!Name.empty()) data["Name"] = Name;
        if (!CAS.empty()) data["CAS"] = CAS;
        if (NormalFreezingPoint != 0.0) data["NormalFreezingPoint"] = NormalFreezingPoint;
        if (NormalBoilingPoint != 0.0) data["NormalBoilingPoint"] = NormalBoilingPoint;
        if (CriticalTemperature != 0.0) data["CriticalTemperature"] = CriticalTemperature;
        if (CriticalPressure != 0.0) data["CriticalPressure"] = CriticalPressure;
        if (Pressure != 0.0) data["Pressure"] = Pressure;
        if (Temperature != 0.0) data["Temperature"] = Temperature;
        if (MolarVolume != 0.0) data["MolarVolume"] = MolarVolume;
        if (MolarWeight != 0.0) data["MolarWeight"] = MolarWeight;
        if (MolarFlow != 0.0) data["MolarFlow"] = MolarFlow;
        if (Compressibility != 0.0) data["Compressibility"] = Compressibility;
        if (FugacityCoefficient != 0.0) data["FugacityCoefficient"] = FugacityCoefficient;
        if (Viscosity != 0.0) data["Viscosity"] = Viscosity;
        if (SurfaceTension != 0.0) data["SurfaceTension"] = SurfaceTension;
        if (ThermalConductivity != 0.0) data["ThermalConductivity"] = ThermalConductivity;
        if (CpDeparture != 0.0) data["CpDeparture"] = CpDeparture;
        if (CvDeparture != 0.0) data["CvDeparture"] = CvDeparture;
        if (EnthalpyDeparture != 0.0) data["EnthalpyDeparture"] = EnthalpyDeparture;
        if (EntropyDeparture != 0.0) data["EntropyDeparture"] = EntropyDeparture;
        if (InternalEnergyDeparture != 0.0) data["InternalEnergyDeparture"] = InternalEnergyDeparture;
        if (GibbsEnergyDeparture != 0.0) data["GibbsEnergyDeparture"] = GibbsEnergyDeparture;
        if (HelmholzEnergyDeparture != 0.0) data["HelmholzEnergyDeparture"] = HelmholzEnergyDeparture;
        if (DPDV != 0.0) data["DPDV"] = DPDV;
        if (DPDT != 0.0) data["DPDT"] = DPDT;
        if (DVDP != 0.0) data["DVDP"] = DVDP;
        if (DVDT != 0.0) data["DVDT"] = DVDT;
        if (DTDV != 0.0) data["DTDV"] = DTDV;
        if (DTDP != 0.0) data["DTDP"] = DTDP;
        if (Cp != 0.0) data["Cp"] = Cp;
        if (Cv != 0.0) data["Cv"] = Cv;
        if (IsothermalCompressibility != 0.0) data["IsothermalCompressibility"] = IsothermalCompressibility;
        if (ThermalExpansionCoefficient != 0.0) data["ThermalExpansionCoefficient"] = ThermalExpansionCoefficient;
        if (JouleThomsonCoefficient != 0.0) data["JouleThomsonCoefficient"] = JouleThomsonCoefficient;
        if (SpeedOfSound != 0.0) data["SpeedOfSound"] = SpeedOfSound;
        if (SaturationPressure != 0.0) data["SaturationPressure"] = SaturationPressure;
        if (SaturationVolume != 0.0) data["SaturationVolume"] = SaturationVolume;
        if (Enthalpy != 0.0) data["Enthalpy"] = Enthalpy;
        if (Entropy != 0.0) data["Entropy"] = Entropy;
        if (InternalEnergy != 0.0) data["InternalEnergy"] = InternalEnergy;
        if (GibbsEnergy != 0.0) data["GibbsEnergy"] = GibbsEnergy;
        if (HelmholzEnergy != 0.0) data["HelmholzEnergy"] = HelmholzEnergy;

        return data.dump();
    }

}    // namespace PCProps
