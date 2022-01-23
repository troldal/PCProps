//
// Created by Kenneth Balslev on 23/01/2022.
//

#include "Valve.hpp"

#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

#include <Common/Globals.hpp>
#include <FluidProperties.hpp>
#include <numeric>

namespace PCProps::UnitOps
{

    using JSONString = std::string;

    /**
     * @brief
     * @todo: Implement viscosity correction
     * @todo: Implement pump curves
     */
    class Valve::impl
    {
        enum class SpecType { OutletPressure, DeltaPressure };

        std::vector<const Stream*>  m_inletStreams;
        mutable std::vector<Stream> m_outletStreams {};

        SpecType       m_specification;
        mutable double m_outletPressure;
        mutable double m_deltaPressure;

        mutable FluidProperties m_fluidProps;

    public:
        /**
         * @brief
         * @param specification
         */
        impl(const JSONString& specification)
        {
            m_outletStreams.reserve(1);

            using rapidjson::Document;
            Document pumpspec;
            pumpspec.Parse(specification.c_str());

            for (const auto& item : pumpspec.GetObject()) {
                std::string key = item.name.GetString();

                if (key == "OutletPressure") {
                    m_specification  = SpecType::OutletPressure;
                    m_outletPressure = item.value.GetDouble();
                }
                else if (key == "DeltaPressure") {
                    m_specification = SpecType::DeltaPressure;
                    m_deltaPressure = item.value.GetDouble();
                }
            }
        }

        /**
         * @brief
         * @return
         */
        //        double computeDensity() const {
        //            return 1.0 / (m_fluidProps.mixtureMolarVolume() / m_fluidProps.mixtureMolarWeight() * 1000);
        //        }

        /**
         * @brief
         * @return
         */
        //        double computeMassFlow() const {
        //            return m_fluidProps.mixtureMolarFlow() * m_fluidProps.mixtureMolarWeight() / 1000.0;
        //        }

        /**
         * @brief
         */
        void calcResults() const
        {
            using namespace PCProps::Globals;

            switch (m_specification) {
                case SpecType::OutletPressure: {
                    m_deltaPressure = m_outletPressure - m_fluidProps.front().Pressure;
                } break;
                case SpecType::DeltaPressure: {
                    m_outletPressure = m_fluidProps.front().Pressure + m_deltaPressure;
                } break;
            }
        }

        /**
         * @brief
         * @param stream
         */
        void setInletStreams(const Stream* stream)
        {
            m_inletStreams = { stream };
        }

        /**
         * @brief
         * @param stream
         */
        void setInletStreams(const std::vector<const Stream*> streams)
        {
            m_inletStreams = streams;
        }

        /**
         * @brief
         * @param specification
         */
        void setSpecification(const JSONString& specification)
        {
            m_outletPressure = 0.0;
            m_deltaPressure  = 0.0;

            using rapidjson::Document;
            Document pumpspec;
            pumpspec.Parse(specification.c_str());

            for (const auto& item : pumpspec.GetObject()) {
                std::string key = item.name.GetString();

                if (key == "OutletPressure") {
                    m_specification  = SpecType::OutletPressure;
                    m_outletPressure = item.value.GetDouble();
                }
                else if (key == "DeltaPressure") {
                    m_specification = SpecType::DeltaPressure;
                    m_deltaPressure = item.value.GetDouble();
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
        void compute()
        {
            using namespace PCProps::Globals;

            m_fluidProps = FluidProperties(m_inletStreams.front()->properties());
            calcResults();

            m_outletStreams = { *m_inletStreams.front() };
            m_outletStreams.front().flash("PH", m_outletPressure, m_fluidProps.mixtureEnthalpy());
        }

        /**
         * @brief
         * @return
         */
        std::string results() const
        {
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
            writer.EndObject();

            return s.GetString();
        }
    };

    /**
     * @details
     */
    Valve::Valve(const Valve::JSONString& specification) : m_impl(std::make_unique<impl>(specification)) {}

    /**
     * @details
     */
    Valve::Valve(const Valve& other) : m_impl(std::make_unique<impl>(*other.m_impl)) {};

    /**
     * @details
     */
    Valve::Valve(Valve&& other) noexcept = default;

    /**
     * @details
     */
    Valve::~Valve() = default;

    /**
     * @details
     */
    Valve& Valve::operator=(const Valve& other)
    {
        Valve copy = other;
        *this      = std::move(copy);
        return *this;
    }

    /**
     * @details
     */
    Valve& Valve::operator=(Valve&& other) noexcept = default;

    /**
     * @details
     */
    void Valve::setInletStreams(const std::vector<const Stream*> streams)
    {
        m_impl->setInletStreams(streams);
    }

    /**
     * @details
     */
    void Valve::setSpecification(const Valve::JSONString& specification)
    {
        m_impl->setSpecification(specification);
    }

    /**
     * @details
     */
    const std::vector<Stream>& Valve::outletStreams()
    {
        return m_impl->outletStreams();
    }

    /**
     * @details
     */
    void Valve::compute()
    {
        m_impl->compute();
    }

    /**
     * @details
     */
    Valve::JSONString Valve::results() const
    {
        return m_impl->results();
    }

}    // namespace PCProps::UnitOps