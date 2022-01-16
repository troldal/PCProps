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


#ifndef PCPROPS_FLUIDPROPERTIES_HPP
#define PCPROPS_FLUIDPROPERTIES_HPP

#include "PhaseProperties.hpp"

// ===== External headers
#include <vector>

namespace PCProps
{
    // ===== Alias declarations
    using JSONString = std::string;

    /**
     * @brief
     */
    class FluidProperties
    {
    private:

        std::vector<PhaseProperties> m_phases;

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
    };

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
