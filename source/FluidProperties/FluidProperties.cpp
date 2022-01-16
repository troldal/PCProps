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


#include "FluidProperties.hpp"

#include <json/json.hpp>

namespace PCProps {

    /**
     * @details
     */
    FluidProperties::FluidProperties() = default;

    /**
     * @details
     */
    FluidProperties::FluidProperties(const std::string& JSONData) {

        auto fluid = nlohmann::json::parse(JSONData);
        for (const auto& phase : fluid) m_phases.emplace_back(phase.dump());
    }

    /**
     * @details
     */
    FluidProperties::FluidProperties(const std::vector<PhaseProperties>& fluidProps) : m_phases(fluidProps) {}

    /**
     * @details
     */
    FluidProperties::FluidProperties(const FluidProperties& other) = default;

    /**
     * @details
     */
    FluidProperties::FluidProperties(FluidProperties&& other) noexcept = default;

    /**
     * @details
     */
    FluidProperties::~FluidProperties() = default;

    /**
     * @details
     */
    FluidProperties& FluidProperties::operator=(const FluidProperties& other) = default;

    /**
     * @details
     */
    FluidProperties& FluidProperties::operator=(FluidProperties&& other) noexcept = default;

    /**
     * @details
     */
    FluidProperties& FluidProperties::operator=(const std::vector<PhaseProperties>& fluidProps)
    {
        m_phases = fluidProps;
        return *this;
    }

    /**
     * @details
     */
    const PhaseProperties& FluidProperties::operator[](int index) const
    {
        return m_phases[index];
    }

    /**
     * @details
     */
    const std::vector<PhaseProperties>& FluidProperties::phases() const
    {
        return m_phases;
    }

    /**
     * @details
     */
    size_t FluidProperties::size() const
    {
        return m_phases.size();
    }

    /**
     * @details
     */
    void FluidProperties::erase(int index) {
        auto iter = m_phases.begin();
        std::advance(iter, index);
        m_phases.erase(iter);
    }

    /**
     * @details
     */
    FluidProperties FluidProperties::stablePhase() const
    {
        return FluidProperties(std::vector {*std::min_element(
            m_phases.begin(),
            m_phases.end(),
            [](const PhaseProperties& a, const PhaseProperties& b){return a.FugacityCoefficient < b.FugacityCoefficient;})});
    }

    /**
     * @details
     */
    FluidProperties FluidProperties::heavyPhase() const
    {
        return FluidProperties(std::vector {*std::min_element(
            m_phases.begin(),
            m_phases.end(),
            [](const PhaseProperties& a, const PhaseProperties& b){return a.MolarVolume < b.MolarVolume;})});
    }

    /**
     * @details
     */
    FluidProperties FluidProperties::lightPhase() const
    {
        return FluidProperties(std::vector {*std::max_element(
            m_phases.begin(),
            m_phases.end(),
            [](const PhaseProperties& a, const PhaseProperties& b){return a.MolarVolume < b.MolarVolume;})});
    }

    /**
     * @details
     */
    JSONString FluidProperties::asJSON() const
    {
        std::vector<nlohmann::json> data;
        for (const auto& phase : m_phases) data.emplace_back(nlohmann::json::parse(phase.asJSON()));
        return nlohmann::json(data).dump();
    }

    /**
     * @details
     */
    std::vector<PhaseProperties>::iterator FluidProperties::begin()
    {
        return m_phases.begin();
    }

    /**
     * @details
     */
    std::vector<PhaseProperties>::iterator FluidProperties::end()
    {
        return m_phases.end();
    }

    /**
     * @details
     */
    PhaseProperties& FluidProperties::back()
    {
        return m_phases.back();
    }

    /**
     * @details
     */
    const PhaseProperties& FluidProperties::back() const
    {
        return m_phases.back();
    }

    /**
     * @details
     */
    PhaseProperties& FluidProperties::front()
    {
        return m_phases.front();
    }

    /**
     * @details
     */
    const PhaseProperties& FluidProperties::front() const
    {
        return m_phases.front();
    }

    /**
     * @details
     */
    void FluidProperties::print(std::ostream& stream) {
        stream << std::setprecision(6) << std::fixed;
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
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << 1 / (phase.MolarVolume / phase.MolarWeight * 1000);
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
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.CpDeparture / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "CvDeparture                [J/kg-K] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.CvDeparture / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "EnthalpyDeparture            [J/kg] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.EnthalpyDeparture / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "EntropyDeparture           [J/kg-K] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.EntropyDeparture / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "InternalEnergyDeparture      [J/kg] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.InternalEnergyDeparture / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "GibbsEnergyDeparture         [J/kg] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.GibbsEnergyDeparture / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "HelmholzEnergyDeparture      [J/kg] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.HelmholzEnergyDeparture / phase.MolarWeight * 1000;
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
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.Cp / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "Cv                         [J/kg-K] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.Cv / phase.MolarWeight * 1000;
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

        stream << "Saturation Pressure            [Pa] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.SaturationPressure;
        stream << std::endl;

        stream << "Saturation Volume           [kg/m3] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << 1 / (phase.SaturationVolume / phase.MolarWeight * 1000);
        stream << std::endl;

        stream << "Enthalpy                     [J/kg] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.Enthalpy / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "Entropy                    [J/kg-K] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.Entropy / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "Internal Energy              [J/kg] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.InternalEnergy / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "Gibbs Energy                 [J/kg] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.GibbsEnergy / phase.MolarWeight * 1000;
        stream << std::endl;

        stream << "Helmholz Energy              [J/kg] : ";
        for (const PhaseProperties& phase : this->phases()) stream << std::right << std::setw(width) << phase.HelmholzEnergy / phase.MolarWeight * 1000;
        stream << std::endl;
    }
}