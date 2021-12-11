//
// Created by Kenneth Balslev on 25/11/2021.
//

#include "PureComponentFactory.hpp"

#include <PropertyLib.hpp>
#include <json/json.hpp>

using namespace nlohmann;

namespace
{

    /**
     * @brief
     * @param obj
     * @return
     */
    std::function<double(double)> createIdealGasCpFunction(const json& obj)
    {
        if (obj["IdealGasCp"]["Equation"] == "DIPPR-100")
            return PCProps::HeatCapacity::Polynomial(PCProps::HeatCapacity::Polynomial::CreateFromDIPPR { obj["IdealGasCp"]["C1"].get<double>(),
                                                                                                          obj["IdealGasCp"]["C2"].get<double>(),
                                                                                                          obj["IdealGasCp"]["C3"].get<double>(),
                                                                                                          obj["IdealGasCp"]["C4"].get<double>(),
                                                                                                          obj["IdealGasCp"]["C5"].get<double>() });

        if (obj["IdealGasCp"]["Equation"] == "DIPPR-107")
            return PCProps::HeatCapacity::AlyLee(PCProps::HeatCapacity::AlyLee::CreateFromDIPPR { obj["IdealGasCp"]["C1"].get<double>(),
                                                                                                  obj["IdealGasCp"]["C2"].get<double>(),
                                                                                                  obj["IdealGasCp"]["C3"].get<double>(),
                                                                                                  obj["IdealGasCp"]["C4"].get<double>(),
                                                                                                  obj["IdealGasCp"]["C5"].get<double>() });
    }

    /**
     * @brief
     * @param obj
     * @return
     */
    std::function<double(double)> createLiquidCpFunction(const json& obj)
    {
        if (obj["LiquidCp"]["Equation"] == "DIPPR-100")
            return PCProps::HeatCapacity::Polynomial(PCProps::HeatCapacity::Polynomial::CreateFromDIPPR { obj["LiquidCp"]["C1"].get<double>(),
                                                                                                          obj["LiquidCp"]["C2"].get<double>(),
                                                                                                          obj["LiquidCp"]["C3"].get<double>(),
                                                                                                          obj["LiquidCp"]["C4"].get<double>(),
                                                                                                          obj["LiquidCp"]["C5"].get<double>() });

        if (obj["LiquidCp"]["Equation"] == "DIPPR-114")
            return PCProps::HeatCapacity::PPDSLiquid(PCProps::HeatCapacity::PPDSLiquid::CreateFromDIPPR { obj["LiquidCp"]["C1"].get<double>(),
                                                                                                          obj["LiquidCp"]["C2"].get<double>(),
                                                                                                          obj["LiquidCp"]["C3"].get<double>(),
                                                                                                          obj["LiquidCp"]["C4"].get<double>(),
                                                                                                          obj["CriticalTemperature"].get<double>() });

        if (obj["LiquidCp"]["Equation"] == "VDI")
            return PCProps::HeatCapacity::PPDSLiquid(PCProps::HeatCapacity::PPDSLiquid::CreateFromPPDS { obj["LiquidCp"]["C1"].get<double>(),
                                                                                                         obj["LiquidCp"]["C2"].get<double>(),
                                                                                                         obj["LiquidCp"]["C3"].get<double>(),
                                                                                                         obj["LiquidCp"]["C4"].get<double>(),
                                                                                                         obj["LiquidCp"]["C5"].get<double>(),
                                                                                                         0.0,
                                                                                                         obj["CriticalTemperature"].get<double>(),
                                                                                                         0.0 });
    }

    /**
     * @brief
     * @param obj
     * @return
     */
    std::function<double(double)> createVaporPressureFunction(const json& obj)
    {
        if (obj["VaporPressure"]["Equation"] == "DIPPR-101")
            return PCProps::VaporPressure::AntoineExtended(PCProps::VaporPressure::AntoineExtended::CreateFromDIPPR { obj["VaporPressure"]["C1"].get<double>(),
                                                                                                                      obj["VaporPressure"]["C2"].get<double>(),
                                                                                                                      obj["VaporPressure"]["C3"].get<double>(),
                                                                                                                      obj["VaporPressure"]["C4"].get<double>(),
                                                                                                                      obj["VaporPressure"]["C5"].get<double>() });
    }

    /**
     * @brief
     * @param obj
     * @return
     */
    std::function<double(double)> createVaporViscosityFunction(const json& obj)
    {
        if (obj["SaturatedVaporViscosity"]["Equation"] == "DIPPR-102")
            return PCProps::Viscosity::DIPPR102({ obj["SaturatedVaporViscosity"]["C1"].get<double>(),
                                                  obj["SaturatedVaporViscosity"]["C2"].get<double>(),
                                                  obj["SaturatedVaporViscosity"]["C3"].get<double>(),
                                                  obj["SaturatedVaporViscosity"]["C4"].get<double>() });
    }

    /**
     * @brief
     * @param obj
     * @return
     */
    std::function<double(double)> createLiquidViscosityFunction(const json& obj)
    {
        if (obj["SaturatedLiquidViscosity"]["Equation"] == "DIPPR-101")
            return PCProps::Viscosity::KirchhoffExtended(PCProps::Viscosity::KirchhoffExtended::CreateFromDIPPR { obj["SaturatedLiquidViscosity"]["C1"].get<double>(),
                                                                                                                  obj["SaturatedLiquidViscosity"]["C2"].get<double>(),
                                                                                                                  obj["SaturatedLiquidViscosity"]["C3"].get<double>(),
                                                                                                                  obj["SaturatedLiquidViscosity"]["C4"].get<double>(),
                                                                                                                  obj["SaturatedLiquidViscosity"]["C5"].get<double>() });
    }

    /**
     * @brief
     * @param obj
     * @return
     */
    std::function<double(double)> createLiquidVolumeFunction(const json& obj)
    {
        if (obj["SaturatedLiquidVolume"]["Equation"] == "DIPPR-105")
            return PCProps::LiquidVolume::Rackett(PCProps::LiquidVolume::Rackett::CreateFromDIPPR { obj["SaturatedLiquidVolume"]["C1"].get<double>(),
                                                                                                    obj["SaturatedLiquidVolume"]["C2"].get<double>(),
                                                                                                    obj["SaturatedLiquidVolume"]["C3"].get<double>(),
                                                                                                    obj["SaturatedLiquidVolume"]["C4"].get<double>() });
    }

}    // namespace

namespace PCProps
{

    /**
     * @brief
     */
    class PureComponentFactory::impl
    {
    public:

        /**
         * @brief
         */
        impl() = default;

        /**
         * @brief
         * @param pcdata
         */
        explicit impl(const JSONString& pcdata) : m_pcdata(json::parse(pcdata)) {}

        /**
         * @brief
         */
        ~impl() = default;

        /**
         * @brief
         * @param other
         */
        impl(const impl& other) = default;

        /**
         * @brief
         * @param other
         */
        impl(impl&& other) noexcept = default;

        /**
         * @brief
         * @param other
         * @return
         */
        impl& operator=(const impl& other) = default;

        /**
         * @brief
         * @param other
         * @return
         */
        impl& operator=(impl&& other) noexcept = default;

        /**
         * @brief
         * @param CAS
         * @return
         */
        PureComponent makeComponent(const std::string& CAS)
        {
            PureComponent result;
            auto          component = *find_if(m_pcdata.begin(), m_pcdata.end(), [&](const json& rec) { return rec["CAS"].get<std::string>() == CAS; });

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

            result.addDataItem("IdealGasCp", createIdealGasCpFunction(component));
            result.addDataItem("LiquidCp", createLiquidCpFunction(component));
            result.addDataItem("VaporPressure", createVaporPressureFunction(component));
            result.addDataItem("SaturatedVaporViscosity", createVaporViscosityFunction(component));
            result.addDataItem("SaturatedLiquidViscosity", createLiquidViscosityFunction(component));
            result.addDataItem("SaturatedLiquidVolume", createLiquidVolumeFunction(component));

            result.addDataItem("CompressedLiquidVolume", LiquidVolume::Thomson(result));
            result.addDataItem("CompressedLiquidViscosity", CompressedLiquidViscosity::Lucas(result));
            result.addDataItem("CompressedVaporViscosity", CompressedVaporViscosity::Lucas(result));

            return result;
        }

    private:
        json m_pcdata;
    };

    // =====
    PureComponentFactory::PureComponentFactory() : m_impl(nullptr) {}

    // =====
    PureComponentFactory::PureComponentFactory(const JSONString& pcdata) : m_impl(std::make_unique<impl>(pcdata)) {}

    // =====
    PureComponentFactory::~PureComponentFactory() = default;

    // =====
    PureComponentFactory::PureComponentFactory(const PureComponentFactory& other) : m_impl(other.m_impl ? std::make_unique<impl>(*other.m_impl) : nullptr) {}

    // =====
    PureComponentFactory::PureComponentFactory(PureComponentFactory&& other) noexcept = default;

    // =====
    PureComponentFactory& PureComponentFactory::operator=(const PureComponentFactory& other)
    {
        PureComponentFactory copy = other;
        *this                     = std::move(copy);
        return *this;
    }

    // =====
    PureComponentFactory& PureComponentFactory::operator=(PureComponentFactory&& other) noexcept = default;

    // =====
    void PureComponentFactory::init(const PureComponentFactory::JSONString& pcdata)
    {
        m_impl = std::make_unique<impl>(pcdata);
    }

    // =====
    PureComponent PureComponentFactory::makeComponent(const std::string& CAS) const
    {
        return m_impl->makeComponent(CAS);
    }

}    // namespace PCProps
