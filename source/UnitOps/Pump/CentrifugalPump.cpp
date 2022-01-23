//
// Created by Kenneth Balslev on 22/01/2021.
//

#include "CentrifugalPump.hpp"

#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

#include <FluidProperties.hpp>
#include <Common/Globals.hpp>
#include <numeric>

namespace PCProps::UnitOps {

    using JSONString = std::string;

    /**
     * @brief
     * @todo: Implement viscosity correction
     * @todo: Implement pump curves
     */
    class CentrifugalPump::impl {

        enum class SpecType { OutletPressure, DeltaPressure, HydraulicPower, EffectivePower };

        SpecType m_pumpSpecType;
        std::vector<const Stream*> m_inletStreams;

        mutable std::vector<Stream> m_outletStreams {};
        double m_pumpEfficiency {};

        mutable FluidProperties m_fluidProps;
        mutable double m_outletPressure;
        mutable double m_deltaPressure;
        mutable double m_pumpHead;
        mutable double m_hydraulicPower;
        mutable double m_effectivePower;

    public:

        /**
         * @brief
         * @param specification
         */
        impl(const JSONString& specification) {

            m_outletStreams.reserve(1);

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
         * @return
         */
        double computeDensity() const {
            return 1.0 / (m_fluidProps.mixtureMolarVolume() / m_fluidProps.mixtureMolarWeight() * 1000);
        }

        /**
         * @brief
         * @return
         */
        double computeMassFlow() const {
            return m_fluidProps.mixtureMolarFlow() * m_fluidProps.mixtureMolarWeight() / 1000.0;
        }

        /**
         * @brief
         */
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

        /**
         * @brief
         * @param stream
         */
        void setInletStreams(const Stream* stream) {
            m_inletStreams = {stream};
        }

        /**
         * @brief
         * @param stream
         */
        void setInletStreams(const std::vector<const Stream*> streams) {
            m_inletStreams = streams;
        }

        /**
         * @brief
         * @param specification
         */
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
        const std::vector<Stream>& outletStreams()
        {
            return m_outletStreams;
        }

        /**
         * @brief
         */
        void compute() {
            using namespace PCProps::Globals;

            m_fluidProps = FluidProperties(m_inletStreams.front()->properties());
            calcResults();

            m_outletStreams = { *m_inletStreams.front() };
            m_outletStreams.front().flash("PH", m_outletPressure, m_fluidProps.mixtureEnthalpy() + m_effectivePower / m_fluidProps.mixtureMolarFlow());
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
     * @brief
     * @param streams
     */
    void CentrifugalPump::setInletStreams(const std::vector<const Stream*> streams) {
        m_impl->setInletStreams(streams);
    }

    /**
     * @details
     */
    void CentrifugalPump::setSpecification(const JSONString& specification)
    {
        m_impl->setSpecification(specification);
    }

    /**
     * @details
     */
    const std::vector<Stream>& CentrifugalPump::outletStreams()
    {
        return m_impl->outletStreams();
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