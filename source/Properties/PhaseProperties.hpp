//
// Created by Kenneth Balslev on 14/11/2021.
//

#ifndef PCPROPS_PHASEPROPERTIES_HPP
#define PCPROPS_PHASEPROPERTIES_HPP

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace PCProps
{
    class PhaseProperties
    {
        using JSONString = std::string;

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

        PhaseProperties();


        explicit PhaseProperties(const std::string& JSONData);

        /**
         * @brief Copy constructor.
         */
        PhaseProperties(const PhaseProperties& other);

        /**
         * @brief Move constructor.
         */
        PhaseProperties(PhaseProperties&& other) noexcept;

        /**
         * @brief Destructor.
         */
        ~PhaseProperties();

        /**
         * @brief Copy assignment operator.
         */
        PhaseProperties& operator=(const PhaseProperties& other);

        /**
         * @brief Move assignment operator.
         */
        PhaseProperties& operator=(PhaseProperties&& other) noexcept;

        /**
         *
         * @return
         */
        JSONString asJSON() const;
    };

    inline std::ostream& operator<<(std::ostream& stream, const PCProps::PhaseProperties& properties)
    {
        return stream << std::setprecision(8) << std::fixed
                      << "Molar Flow                    : " << std::right << std::setw(20) << properties.MolarFlow << std::endl
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
