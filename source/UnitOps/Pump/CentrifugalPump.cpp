//
// Created by Kenneth Balslev on 22/01/2021.
//

#include "CentrifugalPump.hpp"

#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

#include <FluidProperties.hpp>
#include <Common/Globals.hpp>

namespace PCProps::UnitOps {

    using JSONString = std::string;

    class CentrifugalPump::impl {

        enum class SpecType { OutletPressure, DeltaPressure, HydraulicPower, EffectivePower };

        SpecType m_pumpSpecType;
        const Stream* m_inletStream;
        mutable Stream m_outletStream;
        double m_pumpEfficiency {};

        mutable FluidProperties m_fluidProps;
        mutable double m_outletPressure;
        mutable double m_deltaPressure;
        mutable double m_pumpHead;
        mutable double m_hydraulicPower;
        mutable double m_effectivePower;

    public:

        // Head
        // Hydraulic Power
        // Actual Power
        // Viscosity Correction
        // Temperature Rise

        // Required Input:
        // Pump Efficiency

        // Optional Input:
        // Pump Curves

        impl(const JSONString& specification) {

            using rapidjson::Document;
            Document pumpspec;
            pumpspec.Parse(specification.c_str());

            for (const auto& item : pumpspec.GetObject()) {
                std::string key = item.name.GetString();

                if (key == "PumpEfficiency") {
                    m_pumpEfficiency = item.value.GetDouble();
                }
                if (key == "OutletPressure") {
                    m_pumpSpecType = SpecType::OutletPressure;
                    m_outletPressure = item.value.GetDouble();
                }
                else if (key == "DeltaPressure") {
                    m_pumpSpecType = SpecType::DeltaPressure;
                    m_deltaPressure = item.value.GetDouble();
                }
                else if (key == "HydraulicPower") {
                    m_pumpSpecType = SpecType::HydraulicPower;
                    m_hydraulicPower = item.value.GetDouble();
                }
                else if (key == "EffectivePower") {
                    m_pumpSpecType = SpecType::EffectivePower;
                    m_effectivePower = item.value.GetDouble();
                }
            }
        }

        double computeDensity() const {

            double result {};
            for (auto iter = m_fluidProps.begin(); iter != m_fluidProps.end(); ++iter) {
                result += 1.0 / (iter->MolarVolume / iter->MolarWeight * 1000) * iter->MolarFraction;
            }

            return result;
        }

        double computeMassFlow() const {

            double result {};
            for (auto iter = m_fluidProps.begin(); iter != m_fluidProps.end(); ++iter) {
                result += iter->MolarFlow * iter->MolarWeight;
            }

            return result / 1000.0;
        }

        double computeCp() const {

            double result {};
            for (auto iter = m_fluidProps.begin(); iter != m_fluidProps.end(); ++iter) {
                result += (iter->Cp / iter->MolarWeight) * 1000 * iter->MolarFraction;
            }

            return result;
        }

        void calcResults() const {
            using namespace PCProps::Globals;

            switch (m_pumpSpecType) {
                case SpecType::OutletPressure:
                {
                    m_deltaPressure = m_outletPressure - m_fluidProps.front().Pressure;
                    m_pumpHead = m_deltaPressure / (computeDensity() * G_ACCL);
                    m_hydraulicPower = computeMassFlow() * m_pumpHead * G_ACCL;
                    m_effectivePower = m_hydraulicPower / m_pumpEfficiency;
                }
                    break;
                case SpecType::DeltaPressure:
                {
                    m_outletPressure = m_fluidProps.front().Pressure + m_deltaPressure;
                    m_pumpHead = m_deltaPressure / (computeDensity() * G_ACCL);
                    m_hydraulicPower = computeMassFlow() * m_pumpHead * G_ACCL;
                    m_effectivePower = m_hydraulicPower / m_pumpEfficiency;
                }
                    break;
                case SpecType::HydraulicPower:
                {
                    m_effectivePower = m_hydraulicPower / m_pumpEfficiency;
                    m_pumpHead = m_hydraulicPower / (computeMassFlow() * G_ACCL);
                    m_deltaPressure = m_pumpHead * computeDensity() * G_ACCL;
                    m_outletPressure = m_fluidProps.front().Pressure + m_deltaPressure;
                }
                    break;
                case SpecType::EffectivePower:
                {
                    m_hydraulicPower = m_effectivePower * m_pumpEfficiency;
                    m_pumpHead = m_hydraulicPower / (computeMassFlow() * G_ACCL);
                    m_deltaPressure = m_pumpHead * computeDensity() * G_ACCL;
                    m_outletPressure = m_fluidProps.front().Pressure + m_deltaPressure;
                }
                    break;
            }
        }

//        const Stream& operator()() const {
//
//            using namespace PCProps::Globals;
//
//            m_fluidProps = FluidProperties(m_inletStream->properties());
//            calcResults();
//
//            auto temperatureRise = G_ACCL * m_pumpHead * (1.0 / m_pumpEfficiency - 1.0) / computeCp();
//
//            m_outletStream = *m_inletStream;
//            m_outletStream.flash("PT", m_outletPressure, m_fluidProps.front().Temperature + temperatureRise);
//            return m_outletStream;
//        }


        void setInletStream(const Stream* stream) {
            m_inletStream = stream;
        }


        void setSpecification(const JSONString& specification) {

            m_outletPressure = 0.0;
            m_deltaPressure = 0.0;
            m_pumpHead = 0.0;
            m_hydraulicPower = 0.0;
            m_effectivePower =0.0;

            using rapidjson::Document;
            Document pumpspec;
            pumpspec.Parse(specification.c_str());

            for (const auto& item : pumpspec.GetObject()) {
                std::string key = item.name.GetString();

                if (key == "PumpEfficiency") {
                    m_pumpEfficiency = item.value.GetDouble();
                }
                if (key == "OutletPressure") {
                    m_pumpSpecType = SpecType::OutletPressure;
                    m_outletPressure = item.value.GetDouble();
                }
                else if (key == "DeltaPressure") {
                    m_pumpSpecType = SpecType::DeltaPressure;
                    m_deltaPressure = item.value.GetDouble();
                }
                else if (key == "HydraulicPower") {
                    m_pumpSpecType = SpecType::HydraulicPower;
                    m_hydraulicPower = item.value.GetDouble();
                }
                else if (key == "EffectivePower") {
                    m_pumpSpecType = SpecType::EffectivePower;
                    m_effectivePower = item.value.GetDouble();
                }
            }

        }

        /**
         * @brief
         * @param streamName
         * @return
         */
        const Stream& outputStream(const std::string& streamName = "")
        {
            return m_outletStream;
        }

        /**
         * @brief
         */
        void compute() {
            using namespace PCProps::Globals;

            m_fluidProps = FluidProperties(m_inletStream->properties());
            calcResults();

            auto temperatureRise = G_ACCL * m_pumpHead * (1.0 / m_pumpEfficiency - 1.0) / computeCp();

            m_outletStream = *m_inletStream;
            m_outletStream.flash("PT", m_outletPressure, m_fluidProps.front().Temperature + temperatureRise);
        }

        /**
         * @brief
         * @return
         */
        std::string results() const {

            using rapidjson::StringBuffer;
            using rapidjson::Writer;

            StringBuffer         s;
            Writer<StringBuffer> writer(s);
            writer.StartObject();
            writer.Key("InletPressure");
            writer.Double(m_fluidProps.front().Pressure);
            writer.Key("OutletPressure");
            writer.Double(m_outletPressure);
            writer.Key("DeltaPressure");
            writer.Double(m_deltaPressure);
            writer.Key("PumpHead");
            writer.Double(m_pumpHead);
            writer.Key("HydraulicPower");
            writer.Double(m_hydraulicPower);
            writer.Key("EffectivePower");
            writer.Double(m_effectivePower);
            writer.EndObject();

            return s.GetString();
        }

    };

    /**
     * @details
     */
    CentrifugalPump::CentrifugalPump(const JSONString& specification) : m_impl(std::make_unique<impl>(specification)) {}

    /**
     * @details
     */
    CentrifugalPump::CentrifugalPump(const CentrifugalPump& other) : m_impl(std::make_unique<impl>(*other.m_impl)) {};

    /**
     * @details
     */
    CentrifugalPump::CentrifugalPump(CentrifugalPump&& other) noexcept = default;

    /**
     * @details
     */
    CentrifugalPump::~CentrifugalPump() = default;

    /**
     * @details
     */
    CentrifugalPump& CentrifugalPump::operator=(const CentrifugalPump& other)
    {
        CentrifugalPump copy = other;
        *this                = std::move(copy);
        return *this;
    }

    /**
     * @details
     */
    CentrifugalPump& CentrifugalPump::operator=(CentrifugalPump&& other) noexcept = default;

    /**
     * @details
     */
//    const Stream& CentrifugalPump::operator()() const
//    {
//        return m_impl->operator()();
//    }

    /**
     * @details
     */
    void CentrifugalPump::setInletStream(const Stream* stream) {
        m_impl->setInletStream(stream);
    }

    /**
     * @details
     */
    void CentrifugalPump::setSpecification(const JSONString& specification) {}

    /**
     * @details
     */
    const Stream& CentrifugalPump::outputStream(const std::string& streamName)
    {
        return m_impl->outputStream();
    }

    /**
     * @details
     */
    void CentrifugalPump::compute() {
        m_impl->compute();
    }

    /**
     * @details
     */
    std::string CentrifugalPump::results() const
    {
        return m_impl->results();
    }


} //  namespace PCProps::UnitOps