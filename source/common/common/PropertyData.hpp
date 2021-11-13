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
    using PCPhase     = std::array<double, 21>;
    using PCPhases    = std::vector<PCPhase>;

    enum PCPhaseDataElement {
        PCPressure,
        PCTemperature,
        PCMolarVolume,
        PCMolarWeight,
        PCMolarFlow,
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
        PCVaporPressure,
        PCEnthalpy,
        PCEntropy,
        PCInternalEnergy,
        PCGibbsEnergy,
        PCHelmholzEnergy
    };

    inline std::ostream& operator<<(std::ostream& stream, const PCProps::PCPhase& properties)
    {

        return stream << std::setprecision(8) << std::fixed
                      << "Molar Flow                    : " << std::right << std::setw(20) << properties[PCMolarFlow] << std::endl
                      << "Molar Volume                  : " << std::right << std::setw(20) << properties[PCMolarVolume] << " m3/mol" << std::endl
                      << "Surface Tension               : " << std::right << std::setw(20) << properties[PCSurfaceTension] << " N/m" << std::endl
                      << "Thermal Conductivity          : " << std::right << std::setw(20) << properties[PCThermalConductivity] << " W/m-K" << std::endl
                      << "Viscosity                     : " << std::right << std::setw(20) << properties[PCViscosity] << " Pa-s" << std::endl
                      << "Cp                            : " << std::right << std::setw(20) << properties[PCHeatCapacityCp] << " J/mol-K" << std::endl
                      << "Cv                            : " << std::right << std::setw(20) << properties[PCHeatCapacityCv] << " J/mol-K" << std::endl
                      << "Isothermal Compressibility    : " << std::right << std::setw(20) << properties[PCIsothermalCompressibility] << " 1/Pa" << std::endl
                      << "Thermal Expansion Coefficient : " << std::right << std::setw(20) << properties[PCThermalExpansionCoefficient] << " 1/K" << std::endl
                      << "Joule-Thomson Coefficient     : " << std::right << std::setw(20) << properties[PCJouleThomsonCoefficient] << " K/Pa" << std::endl
                      << "Molecular Weight              : " << std::right << std::setw(20) << properties[PCMolarWeight] << " g/mol" << std::endl
                      << "Temperature                   : " << std::right << std::setw(20) << properties[PCTemperature] << " K" << std::endl
                      << "Pressure                      : " << std::right << std::setw(20) << properties[PCPressure] << " Pa" << std::endl
                      << "Compressibility               : " << std::right << std::setw(20) << properties[PCCompressibility] << " -" << std::endl
                      << "Fugacity Coefficient          : " << std::right << std::setw(20) << properties[PCFugacityCoefficient] << " -" << std::endl
                      << "Vapor Pressure                : " << std::right << std::setw(20) << properties[PCVaporPressure] << " Pa" << std::endl
                      << "Enthalpy                      : " << std::right << std::setw(20) << properties[PCEnthalpy] << " J/mol" << std::endl
                      << "Entropy                       : " << std::right << std::setw(20) << properties[PCEntropy] << " J/mol-K" << std::endl
                      << "Internal Energy               : " << std::right << std::setw(20) << properties[PCInternalEnergy] << " J/mol" << std::endl
                      << "Gibbs Energy                  : " << std::right << std::setw(20) << properties[PCGibbsEnergy] << " J/mol" << std::endl
                      << "Helmholz Energy               : " << std::right << std::setw(20) << properties[PCHelmholzEnergy] << " J/mol" << std::endl;

    }

}    // namespace PCProps

#endif    // PCPROPS_PROPERTYDATA_HPP
