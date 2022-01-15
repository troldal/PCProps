//
// Created by Kenneth Balslev on 19/11/2021.
//

#include "PhaseProperties.hpp"

#include <json/json.hpp>

namespace PCProps
{

    /**
     * @brief Default constructor.
     */
    PhaseProperties::PhaseProperties() = default;

    /**
     * @brief Constructor.
     */
    PhaseProperties::PhaseProperties(const std::string& JSONData)
    {
        auto data = nlohmann::json::parse(JSONData);

        auto AsType = [&](const std::string& type) {
            if (type == "VAPOR") return PhaseType::Vapor;
            if (type == "LIQUID") return PhaseType::Liquid;
            return PhaseType::Undefined;
        };

        Type = data["Type"].is_null() ? PhaseType::Undefined : AsType(data["Type"].get<std::string>());
        Name = data["Name"].is_null() ? "" : data["Name"].get<std::string>();
        CAS  = data["CAS"].is_null() ? "" : data["CAS"].get<std::string>();

        NormalFreezingPoint = data["NormalFreezingPoint"].is_null() ? std::nan("") : data["NormalFreezingPoint"].get<double>();
        NormalBoilingPoint  = data["NormalBoilingPoint"].is_null() ? std::nan("") : data["NormalBoilingPoint"].get<double>();
        CriticalTemperature = data["CriticalTemperature"].is_null() ? std::nan("") : data["CriticalTemperature"].get<double>();
        CriticalPressure    = data["CriticalPressure"].is_null() ? std::nan("") : data["CriticalPressure"].get<double>();

        Pressure                    = data["Pressure"].is_null() ? std::nan("") : data["Pressure"].get<double>();
        Temperature                 = data["Temperature"].is_null() ? std::nan("") : data["Temperature"].get<double>();
        MolarVolume                 = data["MolarVolume"].is_null() ? std::nan("") : data["MolarVolume"].get<double>();
        MolarWeight                 = data["MolarWeight"].is_null() ? std::nan("") : data["MolarWeight"].get<double>();
        MolarFlow                   = data["MolarFlow"].is_null() ? std::nan("") : data["MolarFlow"].get<double>();
        Compressibility             = data["Compressibility"].is_null() ? std::nan("") : data["Compressibility"].get<double>();
        FugacityCoefficient         = data["FugacityCoefficient"].is_null() ? std::nan("") : data["FugacityCoefficient"].get<double>();
        Viscosity                   = data["Viscosity"].is_null() ? std::nan("") : data["Viscosity"].get<double>();
        SurfaceTension              = data["SurfaceTension"].is_null() ? std::nan("") : data["SurfaceTension"].get<double>();
        ThermalConductivity         = data["ThermalConductivity"].is_null() ? std::nan("") : data["ThermalConductivity"].get<double>();
        CpDeparture                 = data["CpDeparture"].is_null() ? std::nan("") : data["CpDeparture"].get<double>();
        CvDeparture                 = data["CvDeparture"].is_null() ? std::nan("") : data["CvDeparture"].get<double>();
        EnthalpyDeparture           = data["EnthalpyDeparture"].is_null() ? std::nan("") : data["EnthalpyDeparture"].get<double>();
        EntropyDeparture            = data["EntropyDeparture"].is_null() ? std::nan("") : data["EntropyDeparture"].get<double>();
        InternalEnergyDeparture     = data["InternalEnergyDeparture"].is_null() ? std::nan("") : data["InternalEnergyDeparture"].get<double>();
        GibbsEnergyDeparture        = data["GibbsEnergyDeparture"].is_null() ? std::nan("") : data["GibbsEnergyDeparture"].get<double>();
        HelmholzEnergyDeparture     = data["HelmholzEnergyDeparture"].is_null() ? std::nan("") : data["HelmholzEnergyDeparture"].get<double>();
        DPDV                        = data["DPDV"].is_null() ? std::nan("") : data["DPDV"].get<double>();
        DPDT                        = data["DPDT"].is_null() ? std::nan("") : data["DPDT"].get<double>();
        DVDP                        = data["DVDP"].is_null() ? std::nan("") : data["DVDP"].get<double>();
        DVDT                        = data["DVDT"].is_null() ? std::nan("") : data["DVDT"].get<double>();
        DTDV                        = data["DTDV"].is_null() ? std::nan("") : data["DTDV"].get<double>();
        DTDP                        = data["DTDP"].is_null() ? std::nan("") : data["DTDP"].get<double>();
        Cp                          = data["Cp"].is_null() ? std::nan("") : data["Cp"].get<double>();
        Cv                          = data["Cv"].is_null() ? std::nan("") : data["Cv"].get<double>();
        IsothermalCompressibility   = data["IsothermalCompressibility"].is_null() ? std::nan("") : data["IsothermalCompressibility"].get<double>();
        ThermalExpansionCoefficient = data["ThermalExpansionCoefficient"].is_null() ? std::nan("") : data["ThermalExpansionCoefficient"].get<double>();
        JouleThomsonCoefficient     = data["JouleThomsonCoefficient"].is_null() ? std::nan("") : data["JouleThomsonCoefficient"].get<double>();
        SpeedOfSound                = data["SpeedOfSound"].is_null() ? std::nan("") : data["SpeedOfSound"].get<double>();
        SaturationPressure          = data["SaturationPressure"].is_null() ? std::nan("") : data["SaturationPressure"].get<double>();
        SaturationVolume            = data["SaturationVolume"].is_null() ? std::nan("") : data["SaturationVolume"].get<double>();
        Enthalpy                    = data["Enthalpy"].is_null() ? std::nan("") : data["Enthalpy"].get<double>();
        Entropy                     = data["Entropy"].is_null() ? std::nan("") : data["Entropy"].get<double>();
        InternalEnergy              = data["InternalEnergy"].is_null() ? std::nan("") : data["InternalEnergy"].get<double>();
        GibbsEnergy                 = data["GibbsEnergy"].is_null() ? std::nan("") : data["GibbsEnergy"].get<double>();
        HelmholzEnergy              = data["HelmholzEnergy"].is_null() ? std::nan("") : data["HelmholzEnergy"].get<double>();
    }

    /**
     * @brief Copy constructor.
     */
    PhaseProperties::PhaseProperties(const PhaseProperties& other) = default;

    /**
     * @brief Move constructor.
     */
    PhaseProperties::PhaseProperties(PhaseProperties&& other) noexcept = default;

    /**
     * @brief Destructor.
     */
    PhaseProperties::~PhaseProperties() = default;

    /**
     * @brief Copy assignment operator.
     */
    PhaseProperties& PhaseProperties::operator=(const PhaseProperties& other) = default;

    /**
     * @brief Move assignment operator.
     */
    PhaseProperties& PhaseProperties::operator=(PhaseProperties&& other) noexcept = default;

    /**
     *
     */
    PhaseProperties::JSONString PhaseProperties::asJSON() const

    {
        nlohmann::json data;

        auto TypeAsString = [&](const PhaseType type) {
            if (type == PhaseType::Vapor) return "VAPOR";
            if (type == PhaseType::Liquid) return "LIQUID";
            return "UNDEFINED";
        };

        data["Type"]                        = TypeAsString(Type);
        data["Name"]                        = Name;
        data["CAS"]                         = CAS;
        data["NormalFreezingPoint"]         = NormalFreezingPoint;
        data["NormalBoilingPoint"]          = NormalBoilingPoint;
        data["CriticalTemperature"]         = CriticalTemperature;
        data["CriticalPressure"]            = CriticalPressure;
        data["Pressure"]                    = Pressure;
        data["Temperature"]                 = Temperature;
        data["MolarVolume"]                 = MolarVolume;
        data["MolarWeight"]                 = MolarWeight;
        data["MolarFlow"]                   = MolarFlow;
        data["Compressibility"]             = Compressibility;
        data["FugacityCoefficient"]         = FugacityCoefficient;
        data["Viscosity"]                   = Viscosity;
        data["SurfaceTension"]              = SurfaceTension;
        data["ThermalConductivity"]         = ThermalConductivity;
        data["CpDeparture"]                 = CpDeparture;
        data["CvDeparture"]                 = CvDeparture;
        data["EnthalpyDeparture"]           = EnthalpyDeparture;
        data["EntropyDeparture"]            = EntropyDeparture;
        data["InternalEnergyDeparture"]     = InternalEnergyDeparture;
        data["GibbsEnergyDeparture"]        = GibbsEnergyDeparture;
        data["HelmholzEnergyDeparture"]     = HelmholzEnergyDeparture;
        data["DPDV"]                        = DPDV;
        data["DPDT"]                        = DPDT;
        data["DVDP"]                        = DVDP;
        data["DVDT"]                        = DVDT;
        data["DTDV"]                        = DTDV;
        data["DTDP"]                        = DTDP;
        data["Cp"]                          = Cp;
        data["Cv"]                          = Cv;
        data["IsothermalCompressibility"]   = IsothermalCompressibility;
        data["ThermalExpansionCoefficient"] = ThermalExpansionCoefficient;
        data["JouleThomsonCoefficient"]     = JouleThomsonCoefficient;
        data["SpeedOfSound"]                = SpeedOfSound;
        data["SaturationPressure"]          = SaturationPressure;
        data["SaturationVolume"]            = SaturationVolume;
        data["Enthalpy"]                    = Enthalpy;
        data["Entropy"]                     = Entropy;
        data["InternalEnergy"]              = InternalEnergy;
        data["GibbsEnergy"]                 = GibbsEnergy;
        data["HelmholzEnergy"]              = HelmholzEnergy;

        return data.dump();
    }

}    // namespace PCProps
