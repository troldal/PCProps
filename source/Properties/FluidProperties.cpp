//
// Created by Kenneth Balslev on 20/11/2021.
//

#include "FluidProperties.hpp"

#include <json/json.hpp>

namespace PCProps {

    FluidProperties::FluidProperties() = default;

    FluidProperties::FluidProperties(const std::string& JSONData) {

        auto fluid = nlohmann::json::parse(JSONData);
        for (const auto& phase : fluid) m_phases.emplace_back(phase.dump());
    }

    FluidProperties::FluidProperties(const std::vector<PhaseProperties>& fluidProps) : m_phases(fluidProps) {}

    FluidProperties::FluidProperties(const FluidProperties& other) = default;

    FluidProperties::FluidProperties(FluidProperties&& other) noexcept = default;

    FluidProperties::~FluidProperties() = default;

    FluidProperties& FluidProperties::operator=(const FluidProperties& other) = default;

    FluidProperties& FluidProperties::operator=(FluidProperties&& other) noexcept = default;

    FluidProperties& FluidProperties::operator=(const std::vector<PhaseProperties>& fluidProps)
    {
        m_phases = fluidProps;
        return *this;
    }

    const PhaseProperties& FluidProperties::operator[](int index) const
    {
        return m_phases[index];
    }

    const std::vector<PhaseProperties>& FluidProperties::phases() const
    {
        return m_phases;
    }

    FluidProperties::JSONString FluidProperties::asJSON() const
    {
        std::vector<nlohmann::json> data;
        for (const auto& phase : m_phases) data.emplace_back(nlohmann::json::parse(phase.asJSON()));
        return nlohmann::json(data).dump();
    }

}