//
// Created by Kenneth Balslev on 25/11/2021.
//

#include "PureComponentFactory.hpp"

#include <json/json.hpp>
#include <PropertyLib.hpp>

using namespace nlohmann;

namespace PCProps {

    class PureComponentFactory::impl
    {
    public:

        impl() = default;

        explicit impl(const JSONString& pcdata) : m_pcdata(json::parse(pcdata)) {}

        ~impl() = default;

        impl(const impl& other) = default;

        impl(impl&& other) noexcept = default;

        impl& operator=(const impl& other) = default;

        impl& operator=(impl&& other) noexcept = default;

        PureComponent makeComponent(const std::string& CAS) {

            PureComponent result;
            auto component = *find_if(m_pcdata.begin(), m_pcdata.end(), [&](const json& rec){ return rec["CAS"].get<std::string>() == CAS;});

            result.addDataItem("Name", component["Name"].get<std::string>());
            result.addDataItem("CAS", component["CAS"].get<std::string>());
            result.addDataItem("MolarWeight", component["MolarWeight"].get<double>());
            result.addDataItem("MeltingTemperature", component["MeltingTemperature"].get<double>());
            result.addDataItem("BoilingTemperature", component["BoilingTemperature"].get<double>());
            result.addDataItem("CriticalTemperature", component["CriticalTemperature"].get<double>());
            result.addDataItem("CriticalPressure", component["CriticalPressure"].get<double>());
            result.addDataItem("CriticalVolume", component["CriticalVolume"].get<double>());
            result.addDataItem("CriticalDensity", component["CriticalDensity"].get<double>());
            result.addDataItem("CriticalCompressibility", component["CriticalCompressibility"].get<double>());
            result.addDataItem("AcentricFactor", component["AcentricFactor"].get<double>());
            result.addDataItem("DipoleMoment", component["DipoleMoment"].get<double>());

            result.addDataItem("CompressedLiquidVolume",  LiquidVolume::Thomson(result));
            result.addDataItem("CompressedLiquidViscosity", CompressedLiquidViscosity::Lucas(result));
            result.addDataItem("CompressedVaporViscosity", CompressedVaporViscosity::Lucas(result));

            return result;
        }

    private:
        json m_pcdata;


    };

    PureComponentFactory::PureComponentFactory() : m_impl(nullptr) {}

    PureComponentFactory::PureComponentFactory(const JSONString& pcdata) : m_impl(std::make_unique<impl>(pcdata)) {}

    PureComponentFactory::~PureComponentFactory() = default;

    PureComponentFactory::PureComponentFactory(const PureComponentFactory& other) : m_impl(other.m_impl ? std::make_unique<impl>(*other.m_impl) : nullptr) {}

    PureComponentFactory::PureComponentFactory(PureComponentFactory&& other) noexcept = default;

    PureComponentFactory& PureComponentFactory::operator=(const PureComponentFactory& other)
    {
        PureComponentFactory copy = other;
        *this             = std::move(copy);
        return *this;
    }

    PureComponentFactory& PureComponentFactory::operator=(PureComponentFactory&& other) noexcept = default;

    void PureComponentFactory::init(const PureComponentFactory::JSONString& pcdata) {
        m_impl = std::make_unique<impl>(pcdata);
    }

    PureComponent PureComponentFactory::makeComponent(const std::string& CAS) const
    {
        return m_impl->makeComponent(CAS);
    }

} // namespace PCProps
