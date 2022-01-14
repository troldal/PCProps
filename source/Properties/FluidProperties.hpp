//
// Created by Kenneth Balslev on 20/11/2021.
//

#ifndef PCPROPS_FLUIDPROPERTIES_HPP
#define PCPROPS_FLUIDPROPERTIES_HPP

#include "PhaseProperties.hpp"

#include <vector>

namespace PCProps
{
    /**
     * @brief
     */
    class FluidProperties
    {
    private:
        using JSONString = std::string;
        std::vector<PhaseProperties> m_phases;

    public:

        /**
         * @brief Default constructor.
         */
        FluidProperties();

        /**
         * @brief Constructor taking a JSON string as an argument.
         * @param JSONData Fluid property data, encoded as a JSON string.
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

        FluidProperties stablePhase() const;

        FluidProperties heavyPhase() const;

        FluidProperties lightPhase() const;

        /**
         * @brief Get the fluid data as a JSON array.
         * @return A JSONString (std::string) with the fluid data.
         */
        JSONString asJSON() const;

        std::vector<PhaseProperties>::iterator begin() {
            return m_phases.begin();
        }

        std::vector<PhaseProperties>::iterator end() {
            return m_phases.end();
        }

        PhaseProperties& back() {
            return m_phases.back();
        }

        const PhaseProperties& back() const {
            return m_phases.back();
        }

        PhaseProperties& front() {
            return m_phases.front();
        }

        const PhaseProperties& front() const {
            return m_phases.front();
        }

        /**
         * @brief Prints the fluid data to an ostream object.
         * @param stream An std::ostream object, e.g. std::cout.
         */
        inline void print(std::ostream& stream) {
            stream << std::setprecision(5) << std::fixed;
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
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << 1/(phase.MolarVolume/phase.MolarWeight*1000);
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
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.CpDeparture/phase.MolarWeight*1000;
            stream << std::endl;

            stream << "CvDeparture                [J/kg-K] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.CvDeparture/phase.MolarWeight*1000;
            stream << std::endl;

            stream << "EnthalpyDeparture            [J/kg] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.EnthalpyDeparture/phase.MolarWeight*1000;
            stream << std::endl;

            stream << "EntropyDeparture           [J/kg-K] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.EntropyDeparture/phase.MolarWeight*1000;
            stream << std::endl;

            stream << "InternalEnergyDeparture      [J/kg] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.InternalEnergyDeparture/phase.MolarWeight*1000;
            stream << std::endl;

            stream << "GibbsEnergyDeparture         [J/kg] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.GibbsEnergyDeparture/phase.MolarWeight*1000;
            stream << std::endl;

            stream << "HelmholzEnergyDeparture      [J/kg] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.HelmholzEnergyDeparture/phase.MolarWeight*1000;
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
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.Cp/phase.MolarWeight*1000;
            stream << std::endl;

            stream << "Cv                         [J/kg-K] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.Cv/phase.MolarWeight*1000;
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

            stream << "Vapor Pressure                 [Pa] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.VaporPressure;
            stream << std::endl;

            stream << "Enthalpy                     [J/kg] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.Enthalpy/phase.MolarWeight*1000;
            stream << std::endl;

            stream << "Entropy                    [J/kg-K] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.Entropy/phase.MolarWeight*1000;
            stream << std::endl;

            stream << "Internal Energy              [J/kg] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.InternalEnergy/phase.MolarWeight*1000;
            stream << std::endl;

            stream << "Gibbs Energy                 [J/kg] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.GibbsEnergy/phase.MolarWeight*1000;
            stream << std::endl;

            stream << "Helmholz Energy              [J/kg] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.HelmholzEnergy/phase.MolarWeight*1000;
            stream << std::endl;
        }
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

        stream << "Vapor Pressure                 [Pa] : ";
        for (const PhaseProperties& phase : properties.phases()) stream << std::right << std::setw(20) << phase.VaporPressure;
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
