//
// Created by Kenneth Balslev on 22/01/2022.
//

#include "StreamMixer.hpp"

#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

#include <Common/Globals.hpp>
#include <FluidProperties.hpp>

#include <algorithm>
#include <numeric>

namespace PCProps::UnitOps
{

    using JSONString = std::string;

    class StreamMixer::impl
    {
        enum class SpecType { OutletPressure, DeltaPressure, HydraulicPower, EffectivePower };

        std::vector<const Stream*>  m_inletStreams;
        mutable std::vector<Stream> m_outletStreams {};

    public:
        /**
         * @brief
         * @param specification
         */
        impl(const JSONString& specification)
        {
            m_outletStreams.reserve(1);
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
        void setSpecification(const JSONString& specification) {}

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

//            auto molarFlow = std::accumulate(m_inletStreams.begin(), m_inletStreams.end(), [](const Stream& stream) {
//                auto properties = FluidProperties(stream.properties());
//                return std::accumulate(properties.begin(), properties.end(), [](const PhaseProperties& phase) { return phase.MolarFlow; });
//            } );
//
//            auto pressure = std::min_element(m_inletStreams.begin(), m_inletStreams.end(), [](const Stream* s1, const Stream* s2) {
//                return FluidProperties(s1->properties()).front().Pressure < FluidProperties(s2->properties()).front().Pressure;
//            });

            Stream result = *m_inletStreams.front();


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

            writer.EndObject();

            return s.GetString();
        }
    };

    /**
     * @details
     */
    StreamMixer::StreamMixer(const JSONString& specification) : m_impl(std::make_unique<impl>(specification)) {}

    /**
     * @details
     */
    StreamMixer::StreamMixer(const StreamMixer& other) : m_impl(std::make_unique<impl>(*other.m_impl)) {};

    /**
     * @details
     */
    StreamMixer::StreamMixer(StreamMixer&& other) noexcept = default;

    /**
     * @details
     */
    StreamMixer::~StreamMixer() = default;

    /**
     * @details
     */
    StreamMixer& StreamMixer::operator=(const StreamMixer& other)
    {
        StreamMixer copy = other;
        *this                = std::move(copy);
        return *this;
    }

    /**
     * @details
     */
    StreamMixer& StreamMixer::operator=(StreamMixer&& other) noexcept = default;

    /**
     * @details
     */
    void StreamMixer::setInletStreams(const std::vector<const Stream*> streams) {
        m_impl->setInletStreams(streams);
    }

    /**
     * @details
     */
    void StreamMixer::setSpecification(const JSONString& specification)
    {
        m_impl->setSpecification(specification);
    }

    /**
     * @details
     */
    const std::vector<Stream>& StreamMixer::outletStreams()
    {
        return m_impl->outletStreams();
    }

    /**
     * @details
     */
    void StreamMixer::compute() {
        m_impl->compute();
    }

    /**
     * @details
     */
    JSONString StreamMixer::results() const
    {
        return m_impl->results();
    }

}    // namespace PCProps::UnitOps