//
// Created by Kenneth Balslev on 14/11/2021.
//

#ifndef PCPROPS_PHASEPROPERTIES_HPP
#define PCPROPS_PHASEPROPERTIES_HPP

#include <iomanip>
#include <iostream>

#include <json/json.hpp>

namespace PCProps
{
    class PCPhaseProperties
    {
        using json = nlohmann::json;
    public:

        double Pressure { 0.0 };
        double Temperature { 0.0 };
        double MolarVolume { 0.0 };
        double MolarWeight { 0.0 };
        double MolarFlow { 0.0 };
        double Compressibility { 0.0 };
        double FugacityCoefficient { 0.0 };
        double Viscosity { 0.0 };
        double SurfaceTension { 0.0 };
        double ThermalConductivity { 0.0 };
        double Cp { 0.0 };
        double Cv { 0.0 };
        double IsothermalCompressibility { 0.0 };
        double ThermalExpansionCoefficient { 0.0 };
        double JouleThomsonCoefficient { 0.0 };
        double VaporPressure { 0.0 };
        double Enthalpy { 0.0 };
        double Entropy { 0.0 };
        double InternalEnergy { 0.0 };
        double GibbsEnergy { 0.0 };
        double HelmholzEnergy { 0.0 };

        PCPhaseProperties() = default;

        explicit PCPhaseProperties(const nlohmann::json& data) {
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
            VaporPressure               = data["VaporPressure"].is_null() ? std::nan("") : data["VaporPressure"].get<double>();
            Enthalpy                    = data["Enthalpy"].is_null() ? std::nan("") : data["Enthalpy"].get<double>();
            Entropy                     = data["Entropy"].is_null() ? std::nan("") : data["Entropy"].get<double>();
            InternalEnergy              = data["InternalEnergy"].is_null() ? std::nan("") : data["InternalEnergy"].get<double>();
            GibbsEnergy                 = data["GibbsEnergy"].is_null() ? std::nan("") : data["GibbsEnergy"].get<double>();
            HelmholzEnergy              = data["HelmholzEnergy"].is_null() ? std::nan("") : data["HelmholzEnergy"].get<double>();
        }

        explicit PCPhaseProperties(const std::string& data) : PCPhaseProperties(nlohmann::json::parse(data)) {}

        /**
         * @brief Copy constructor.
         */
        PCPhaseProperties(const PCPhaseProperties& other) = default;

        /**
         * @brief Move constructor.
         */
        PCPhaseProperties(PCPhaseProperties&& other) noexcept = default;

        /**
         * @brief Destructor.
         */
        ~PCPhaseProperties() = default;

        /**
         * @brief Copy assignment operator.
         */
        PCPhaseProperties& operator=(const PCPhaseProperties& other) = default;

        /**
         * @brief Move assignment operator.
         */
        PCPhaseProperties& operator=(PCPhaseProperties&& other) noexcept = default;

        json asJSON() const
        {
            json data;

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
            data["VaporPressure"]               = VaporPressure;
            data["Enthalpy"]                    = Enthalpy;
            data["Entropy"]                     = Entropy;
            data["InternalEnergy"]              = InternalEnergy;
            data["GibbsEnergy"]                 = GibbsEnergy;
            data["HelmholzEnergy"]              = HelmholzEnergy;

            return data;
        }
    };

        inline std::ostream& operator<<(std::ostream& stream, const PCProps::PCPhaseProperties& properties)
        {
            return stream << std::setprecision(8) << std::fixed << "Molar Flow                    : " << std::right << std::setw(20) << properties.MolarFlow << std::endl
                          << "Molar Volume                  : " << std::right << std::setw(20) << properties.MolarVolume << " m3/mol" << std::endl
                          << "Surface Tension               : " << std::right << std::setw(20) << properties.SurfaceTension << " N/m" << std::endl
                          << "Thermal Conductivity          : " << std::right << std::setw(20) << properties.ThermalConductivity << " W/m-K" << std::endl
                          << "Viscosity                     : " << std::right << std::setw(20) << properties.Viscosity << " Pa-s" << std::endl
                          << "Cp                            : " << std::right << std::setw(20) << properties.Cp << " J/mol-K" << std::endl
                          << "Cv                            : " << std::right << std::setw(20) << properties.Cv << " J/mol-K" << std::endl
                          << "Isothermal Compressibility    : " << std::right << std::setw(20) << properties.IsothermalCompressibility << " 1/Pa" << std::endl
                          << "Thermal Expansion Coefficient : " << std::right << std::setw(20) << properties.ThermalExpansionCoefficient << " 1/K" << std::endl
                          << "Joule-Thomson Coefficient     : " << std::right << std::setw(20) << properties.JouleThomsonCoefficient << " K/Pa" << std::endl
                          << "Molecular Weight              : " << std::right << std::setw(20) << properties.MolarWeight << " g/mol" << std::endl
                          << "Temperature                   : " << std::right << std::setw(20) << properties.Temperature << " K" << std::endl
                          << "Pressure                      : " << std::right << std::setw(20) << properties.Pressure << " Pa" << std::endl
                          << "Compressibility               : " << std::right << std::setw(20) << properties.Compressibility << " -" << std::endl
                          << "Fugacity Coefficient          : " << std::right << std::setw(20) << properties.FugacityCoefficient << " -" << std::endl
                          << "Vapor Pressure                : " << std::right << std::setw(20) << properties.VaporPressure << " Pa" << std::endl
                          << "Enthalpy                      : " << std::right << std::setw(20) << properties.Enthalpy << " J/mol" << std::endl
                          << "Entropy                       : " << std::right << std::setw(20) << properties.Entropy << " J/mol-K" << std::endl
                          << "Internal Energy               : " << std::right << std::setw(20) << properties.InternalEnergy << " J/mol" << std::endl
                          << "Gibbs Energy                  : " << std::right << std::setw(20) << properties.GibbsEnergy << " J/mol" << std::endl
                          << "Helmholz Energy               : " << std::right << std::setw(20) << properties.HelmholzEnergy << " J/mol" << std::endl;
        }

}    // namespace PCProps
#endif    // PCPROPS_PHASEPROPERTIES_HPP
