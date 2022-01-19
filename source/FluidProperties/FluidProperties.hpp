/*

888    d8P  8888888b.
888   d8P   888   Y88b
888  d8P    888    888
888d88K     888   d88P 888d888 .d88b.  88888b.  .d8888b
8888888b    8888888P"  888P"  d88""88b 888 "88b 88K
888  Y88b   888        888    888  888 888  888 "Y8888b.
888   Y88b  888        888    Y88..88P 888 d88P      X88
888    Y88b 888        888     "Y88P"  88888P"   88888P'
                                       888
                                       888
                                       888

Copyright (c) 2022 Kenneth Troldal Balslev

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

#ifndef PCPROPS_FLUIDPROPERTIES_HPP
#define PCPROPS_FLUIDPROPERTIES_HPP

// ===== Standard Library headers ===== //
#include <string>
#include <iomanip>
#include <iostream>
#include <vector>

namespace PCProps
{

    // ===== Enum definitions ===== //
    enum class PhaseType { Vapor, Liquid, Undefined };

    /**
     * @brief The PhaseProperties class encapsulates the properties of a single phase of a fluid (liquid or vapor).
     * This class is essentially just a struct with fields corresponding to the individual property elements
     */
    class PhaseProperties
    {
    public:
        PhaseType Type { PhaseType::Undefined };

        std::string Name {};
        std::string CAS {};

        double NormalFreezingPoint { 0.0 };
        double NormalBoilingPoint { 0.0 };
        double CriticalTemperature { 0.0 };
        double CriticalPressure { 0.0 };

        double Pressure { 0.0 };
        double Temperature { 0.0 };
        double MolarVolume { 0.0 };
        double MolarWeight { 0.0 };
        double MolarFlow { 0.0 };
        double Compressibility { 0.0 };
        double FugacityCoefficient { 0.0 };
        double CpDeparture { 0.0 };
        double CvDeparture { 0.0 };
        double EnthalpyDeparture { 0.0 };
        double EntropyDeparture { 0.0 };
        double InternalEnergyDeparture { 0.0 };
        double GibbsEnergyDeparture { 0.0 };
        double HelmholzEnergyDeparture { 0.0 };
        double DPDV { 0.0 };
        double DPDT { 0.0 };
        double DVDP { 0.0 };
        double DVDT { 0.0 };
        double DTDV { 0.0 };
        double DTDP { 0.0 };
        double Viscosity { 0.0 };
        double SurfaceTension { 0.0 };
        double ThermalConductivity { 0.0 };
        double Cp { 0.0 };
        double Cv { 0.0 };
        double IsothermalCompressibility { 0.0 };
        double ThermalExpansionCoefficient { 0.0 };
        double JouleThomsonCoefficient { 0.0 };
        double SpeedOfSound { 0.0 };
        double SaturationPressure { 0.0 };
        double SaturationVolume { 0.0 };
        double Enthalpy { 0.0 };
        double Entropy { 0.0 };
        double InternalEnergy { 0.0 };
        double GibbsEnergy { 0.0 };
        double HelmholzEnergy { 0.0 };

        /**
         * @brief Default constructor.
         */
        PhaseProperties() = default;

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
    };


    /**
     * @brief
     */
    class FluidProperties
    {
        using JSONString = std::string;

    public:

        /**
         * @brief Default constructor.
         */
        FluidProperties();

        /**
         * @brief Constructor taking a JSON string as an argument.
         * @param JSONData PropertyPackage property data, encoded as a JSON string.
         */
        explicit FluidProperties(const JSONString& JSONData);

        /**
         * @brief Constructor taking a std::vector of PhaseProperties objects as an argument.
         * @param fluidProps a std::vector of PhaseProperties.
         */
        explicit FluidProperties(const std::vector<PhaseProperties>& fluidProps);

        /**
         * @brief Copy constructor.
         */
        FluidProperties(const FluidProperties& other);

        /**
         * @brief Move constructor.
         */
        FluidProperties(FluidProperties&& other) noexcept;

        /**
         * @brief Destructor.
         */
        ~FluidProperties();

        /**
         * @brief Copy assignment operator.
         */
        FluidProperties& operator=(const FluidProperties& other);

        /**
         * @brief Move assignment operator.
         */
        FluidProperties& operator=(FluidProperties&& other) noexcept;

        /**
         * @brief Assignment operator, taking a std::vector of PhaseProperties objects as an argument.
         * @param fluidProps a std::vector of PhaseProperties.
         */
        FluidProperties& operator=(const std::vector<PhaseProperties>& fluidProps);

        /**
         * @brief Array index operator.
         * @param index Index of PhaseProperties object to retrieve.
         * @return The PhaseProperties at the given index.
         */
        const PhaseProperties& operator[](int index) const;

        /**
         * @brief Get a std::vector with the PhaseProperties objects.
         * @return A reference to the std::vector of PhaseProperties objects.
         */
        const std::vector<PhaseProperties>& phases() const;

        /**
         * @brief Get the size/count of fluid phases.
         * @return The size/count of flud phases.
         */
        size_t size() const;

        /**
         * @brief
         * @param index
         */
        void erase(int index);

        /**
         * @brief
         * @return
         */
        FluidProperties stablePhase() const;

        /**
         * @brief
         * @return
         */
        FluidProperties heavyPhase() const;

        /**
         * @brief
         * @return
         */
        FluidProperties lightPhase() const;

        /**
         * @brief Get the fluid data as a JSON array.
         * @return A JSONString (std::string) with the fluid data.
         */
        JSONString asJSON() const;

        /**
         * @brief
         * @return
         */
        std::vector<PhaseProperties>::iterator begin();

        /**
         * @brief
         * @return
         */
        std::vector<PhaseProperties>::iterator end();

        /**
         * @brief
         * @return
         */
        PhaseProperties& back();

        /**
         * @brief
         * @return
         */
        const PhaseProperties& back() const;

        /**
         * @brief
         * @return
         */
        PhaseProperties& front();

        /**
         * @brief
         * @return
         */
        const PhaseProperties& front() const;

        /**
         * @brief Prints the fluid data to an ostream object.
         * @param stream An std::ostream object, e.g. std::cout.
         */
        void print(std::ostream& stream);

    private:
        class impl;
        std::unique_ptr<impl> m_impl;
    };

    /**
     * @brief
     * @param stream
     * @param properties
     * @return
     */
    inline std::ostream& operator<<(std::ostream& stream, const PhaseProperties& properties)
    {
        auto TypeAsString = [&](const PhaseType type) {
            if (type == PhaseType::Vapor) return "VAPOR";
            if (type == PhaseType::Liquid) return "LIQUID";
            return "UNDEFINED";
        };

        return stream << std::setprecision(8) << std::fixed << "Type                          : " << std::right << std::setw(20) << TypeAsString(properties.Type) << std::endl
                      << "Molar Flow                    : " << std::right << std::setw(20) << properties.MolarFlow << std::endl
                      << "Molar Volume                  : " << std::right << std::setw(20) << properties.MolarVolume << " m3/mol" << std::endl
                      << "Surface Tension               : " << std::right << std::setw(20) << properties.SurfaceTension << " N/m" << std::endl
                      << "Thermal Conductivity          : " << std::right << std::setw(20) << properties.ThermalConductivity << " W/m-K" << std::endl
                      << "Viscosity                     : " << std::right << std::setw(20) << properties.Viscosity << " Pa-s" << std::endl
                      << "CpDeparture                   : " << std::right << std::setw(20) << properties.CpDeparture << " J/mol-K" << std::endl
                      << "CvDeparture                   : " << std::right << std::setw(20) << properties.CvDeparture << " J/mol-K" << std::endl
                      << "EnthalpyDeparture             : " << std::right << std::setw(20) << properties.EnthalpyDeparture << " J/mol" << std::endl
                      << "EntropyDeparture              : " << std::right << std::setw(20) << properties.EntropyDeparture << " J/mol-K" << std::endl
                      << "InternalEnergyDeparture       : " << std::right << std::setw(20) << properties.InternalEnergyDeparture << " J/mol" << std::endl
                      << "GibbsEnergyDeparture          : " << std::right << std::setw(20) << properties.GibbsEnergyDeparture << " J/mol" << std::endl
                      << "HelmholzEnergyDeparture       : " << std::right << std::setw(20) << properties.HelmholzEnergyDeparture << " J/mol" << std::endl
                      << "DPDV                          : " << std::right << std::setw(20) << properties.DPDV << " Pa/m3" << std::endl
                      << "DPDT                          : " << std::right << std::setw(20) << properties.DPDT << " Pa/K" << std::endl
                      << "DVDP                          : " << std::right << std::setw(20) << properties.DVDP << " m3/Pa" << std::endl
                      << "DVDT                          : " << std::right << std::setw(20) << properties.DVDP << " m3/K" << std::endl
                      << "DTDV                          : " << std::right << std::setw(20) << properties.DTDV << " K/m3" << std::endl
                      << "DTDP                          : " << std::right << std::setw(20) << properties.DTDP << " K/Pa" << std::endl
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
                      << "Saturation Pressure           : " << std::right << std::setw(20) << properties.SaturationPressure << " Pa" << std::endl
                      << "Saturation Volume             : " << std::right << std::setw(20) << properties.SaturationVolume << " m3/mol" << std::endl
                      << "Enthalpy                      : " << std::right << std::setw(20) << properties.Enthalpy << " J/mol" << std::endl
                      << "Entropy                       : " << std::right << std::setw(20) << properties.Entropy << " J/mol-K" << std::endl
                      << "Internal Energy               : " << std::right << std::setw(20) << properties.InternalEnergy << " J/mol" << std::endl
                      << "Gibbs Energy                  : " << std::right << std::setw(20) << properties.GibbsEnergy << " J/mol" << std::endl
                      << "Helmholz Energy               : " << std::right << std::setw(20) << properties.HelmholzEnergy << " J/mol" << std::endl;
    }

    /**
     * @brief
     * @param stream
     * @param properties
     * @return
     */
    inline std::ostream& operator<<(std::ostream& stream, const FluidProperties& properties) {
        stream << std::setprecision(8) << std::fixed;

        auto TypeAsString = [&](const PhaseType type) {
            if (type == PhaseType::Vapor) return "VAPOR";
            if (type == PhaseType::Liquid) return "LIQUID";
            return "UNDEFINED";
        };

        stream << "Type                                : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << TypeAsString(phase.Type);
        stream << std::endl;

        stream << "Molar Flow                      [-] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.MolarFlow;
        stream << std::endl;

        stream << "Molar Volume               [m3/mol] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.MolarVolume;
        stream << std::endl;

        stream << "Surface Tension               [N/m] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.SurfaceTension;
        stream << std::endl;

        stream << "Thermal Conductivity        [W/m-K] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.ThermalConductivity;
        stream << std::endl;

        stream << "Viscosity                    [Pa-s] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.Viscosity;
        stream << std::endl;

        stream << "Cp                        [J/mol-K] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.Cp;
        stream << std::endl;

        stream << "Cv                        [J/mol-K] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.Cv;
        stream << std::endl;

        stream << "Isothermal Compressibility   [1/Pa] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.IsothermalCompressibility;
        stream << std::endl;

        stream << "Thermal Expansion Coefficient [1/K] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.ThermalExpansionCoefficient;
        stream << std::endl;

        stream << "Joule-Thomson Coefficient    [K/Pa] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.JouleThomsonCoefficient;
        stream << std::endl;

        stream << "Speed of Sound                [m/s] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.SpeedOfSound;
        stream << std::endl;

        stream << "Molecular Weight            [g/mol] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.ThermalExpansionCoefficient;
        stream << std::endl;

        stream << "Temperature                     [K] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.Temperature;
        stream << std::endl;

        stream << "Pressure                       [Pa] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.Pressure;
        stream << std::endl;

        stream << "Compressibility                 [-] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.Compressibility;
        stream << std::endl;

        stream << "Fugacity Coefficient            [-] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.FugacityCoefficient;
        stream << std::endl;

        stream << "Saturation Pressure            [Pa] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.SaturationPressure;
        stream << std::endl;

        stream << "Saturation Volume          [m3/mol] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.SaturationVolume;
        stream << std::endl;

        stream << "Enthalpy                    [J/mol] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.Enthalpy;
        stream << std::endl;

        stream << "Entropy                   [J/mol-K] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.Entropy;
        stream << std::endl;

        stream << "Internal Energy             [J/mol] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.InternalEnergy;
        stream << std::endl;

        stream << "Gibbs Energy                [J/mol] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.GibbsEnergy;
        stream << std::endl;

        stream << "Helmholz Energy             [J/mol] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.HelmholzEnergy;
        stream << std::endl;

        return stream;
    }
}


#endif    // PCPROPS_FLUIDPROPERTIES_HPP
