//
// Created by Kenneth Balslev on 20/11/2021.
//

#ifndef PCPROPS_FLUIDPROPERTIES_HPP
#define PCPROPS_FLUIDPROPERTIES_HPP

#include "PhaseProperties.hpp"

#include <vector>

namespace PCProps
{
    class FluidProperties
    {
        using JSONString = std::string;
        std::vector<PhaseProperties> m_phases;

    public:

        /**
         *
         */
        FluidProperties();

        /**
         *
         * @param JSONData
         */
        explicit FluidProperties(const std::string& JSONData);

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
         *
         * @param index
         * @return
         */
        const PhaseProperties& operator[](int index) const;

        /**
         *
         * @return
         */
        const std::vector<PhaseProperties>& phases() const;

        /**
         *
         * @return
         */
        JSONString asJSON() const;

        inline void print(std::ostream& stream) {
            stream << std::setprecision(8) << std::fixed;

            stream << "Molar Flow                      [-] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.MolarFlow;
            stream << std::endl;

            stream << "Molar Volume                [kg/m3] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << 1/(phase.MolarVolume/phase.MolarWeight*1000);
            stream << std::endl;

            stream << "Surface Tension               [N/m] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.SurfaceTension;
            stream << std::endl;

            stream << "Thermal Conductivity        [W/m-K] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.ThermalConductivity;
            stream << std::endl;

            stream << "Viscosity                    [Pa-s] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.Viscosity;
            stream << std::endl;

            stream << "Cp                         [J/kg-K] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.Cp/phase.MolarWeight*1000;
            stream << std::endl;

            stream << "Cv                         [J/kg-K] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.Cv/phase.MolarWeight*1000;
            stream << std::endl;

            stream << "Isothermal Compressibility   [1/Pa] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.IsothermalCompressibility;
            stream << std::endl;

            stream << "Thermal Expansion Coefficient [1/K] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.ThermalExpansionCoefficient;
            stream << std::endl;

            stream << "Joule-Thomson Coefficient    [K/Pa] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.JouleThomsonCoefficient;
            stream << std::endl;

            stream << "Molecular Weight            [g/mol] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.MolarWeight;
            stream << std::endl;

            stream << "Temperature                     [K] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.Temperature;
            stream << std::endl;

            stream << "Pressure                       [Pa] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.Pressure;
            stream << std::endl;

            stream << "Compressibility                 [-] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.Compressibility;
            stream << std::endl;

            stream << "Fugacity Coefficient            [-] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.FugacityCoefficient;
            stream << std::endl;

            stream << "Vapor Pressure                 [Pa] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.VaporPressure;
            stream << std::endl;

            stream << "Enthalpy                     [J/kg] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.Enthalpy/phase.MolarWeight*1000;
            stream << std::endl;

            stream << "Entropy                    [J/kg-K] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.Entropy/phase.MolarWeight*1000;
            stream << std::endl;

            stream << "Internal Energy              [J/kg] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.InternalEnergy/phase.MolarWeight*1000;
            stream << std::endl;

            stream << "Gibbs Energy                 [J/kg] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.GibbsEnergy/phase.MolarWeight*1000;
            stream << std::endl;

            stream << "Helmholz Energy              [J/kg] : ";
            for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(20) << phase.HelmholzEnergy/phase.MolarWeight*1000;
            stream << std::endl;
        }
    };

    inline std::ostream& operator<<(std::ostream& stream, const FluidProperties& properties) {
        stream << std::setprecision(8) << std::fixed;

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
