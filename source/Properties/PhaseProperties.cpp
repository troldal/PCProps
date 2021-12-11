//
// Created by Kenneth Balslev on 19/11/2021.
//

#include "PhaseProperties.hpp"

#include <json/json.hpp>

namespace PCProps {

    /**
     * @brief Default constructor.
     */
    PhaseProperties::PhaseProperties() = default;

    /**
     * @brief Constructor.
     */
    PhaseProperties::PhaseProperties(const std::string& JSONData) {

        auto data = nlohmann::json::parse(JSONData);

        auto AsType = [&](const std::string& type) {
            if (type == "VAPOR") return PhaseType::Vapor;
            if (type == "LIQUID") return PhaseType::Liquid;
            return PhaseType::Undefined;
        };

        Type                        = data["Type"].is_null() ? PhaseType::Undefined : AsType(data["Type"].get<std::string>());
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
        Cp                          = data["Cp"].is_null() ? std::nan("") : data["Cp"].get<double>();
        Cv                          = data["Cv"].is_null() ? std::nan("") : data["Cv"].get<double>();
        IsothermalCompressibility   = data["IsothermalCompressibility"].is_null() ? std::nan("") : data["IsothermalCompressibility"].get<double>();
        ThermalExpansionCoefficient = data["ThermalExpansionCoefficient"].is_null() ? std::nan("") : data["ThermalExpansionCoefficient"].get<double>();
        JouleThomsonCoefficient     = data["JouleThomsonCoefficient"].is_null() ? std::nan("") : data["JouleThomsonCoefficient"].get<double>();
        SpeedOfSound                = data["SpeedOfSound"].is_null() ? std::nan("") : data["SpeedOfSound"].get<double>();
        VaporPressure               = data["VaporPressure"].is_null() ? std::nan("") : data["VaporPressure"].get<double>();
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
            data["Cp"]                          = Cp;
            data["Cv"]                          = Cv;
            data["IsothermalCompressibility"]   = IsothermalCompressibility;
            data["ThermalExpansionCoefficient"] = ThermalExpansionCoefficient;
            data["JouleThomsonCoefficient"]     = JouleThomsonCoefficient;
            data["SpeedOfSound"]                = SpeedOfSound;
            data["VaporPressure"]               = VaporPressure;
            data["Enthalpy"]                    = Enthalpy;
            data["Entropy"]                     = Entropy;
            data["InternalEnergy"]              = InternalEnergy;
            data["GibbsEnergy"]                 = GibbsEnergy;
            data["HelmholzEnergy"]              = HelmholzEnergy;

            return data.dump();
        }

}

