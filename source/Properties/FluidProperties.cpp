//
// Created by Kenneth Balslev on 20/11/2021.
//

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

    FluidProperties FluidProperties::stablePhase() const
    {
        return FluidProperties(std::vector {*std::min_element(
            m_phases.begin(),
            m_phases.end(),
            [](const PhaseProperties& a, const PhaseProperties& b){return a.FugacityCoefficient < b.FugacityCoefficient;})});
    }

    FluidProperties FluidProperties::heavyPhase() const
    {
        return FluidProperties(std::vector {*std::min_element(
            m_phases.begin(),
            m_phases.end(),
            [](const PhaseProperties& a, const PhaseProperties& b){return a.MolarVolume < b.MolarVolume;})});
    }

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
    FluidProperties::JSONString FluidProperties::asJSON() const
    {
        std::vector<nlohmann::json> data;
        for (const auto& phase : m_phases) data.emplace_back(nlohmann::json::parse(phase.asJSON()));
        return nlohmann::json(data).dump();
    }
}