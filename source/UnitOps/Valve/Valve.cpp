/*

888    d8P  8888888b.
888   d8P   888   Y88b
888  d8P    888    888
888d88K     888   d88P 888d888 .d88b.  88888b.  .d8888b
8888888b    8888888P"  888P"  d88""88b 888 "88b 88K
888  Y88b   888        888    888  888 888  888 "Y8888b.
888   Y88b  888        888    Y88..88P 888 d88P      X88
888    Y88b 888        888     "Y88P"  88888P"   88888P'
                                       888
                                       888
                                       888

Copyright (c) 2022 Kenneth Troldal Balslev

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

// ===== KProps Headers ===== //
#include "Valve.hpp"
#include <FluidProperties.hpp>

// ===== External Headers ===== //
#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

// ===== Standard Library Headers ===== //
#include <numeric>

namespace PCProps::UnitOps
{

    using JSONString = std::string;

    /**
     * @brief The Valve::impl class is a nested (hidden) implementation class for the UnitOps::Valve class.
     */
    class Valve::impl
    {
        // ===== Specification type enum
        enum class SpecType { OutletPressure, DeltaPressure };

        // ===== Equipment specification data
        SpecType       m_specification;
        mutable double m_outletPressure;
        mutable double m_deltaPressure;

        // ===== Inlet/outlet streams
        std::vector<const Stream*>  m_inletStreams;
        mutable std::vector<Stream> m_outletStreams {};

        /**
         * @brief Helper function for getting the fluid properties from the inlet stream.
         * @return A FluidProperties object with the fluid properties.
         */
        FluidProperties inletFluidProps() const {
            return FluidProperties(m_inletStreams.front()->results());
        }

        /**
         * @brief Calculates the equipment data (not the outlet stream). I.e., if the outlet pressure
         * is given, calculate the delta pressure, or vice versa.
         */
        void calcResults() const
        {
            switch (m_specification) {
                case SpecType::OutletPressure: {
                    m_deltaPressure = m_outletPressure - inletFluidProps().front().Pressure;
                } break;
                case SpecType::DeltaPressure: {
                    m_outletPressure = inletFluidProps().front().Pressure + m_deltaPressure;
                } break;
            }
        }

    public:

        /**
         * @brief Constructor taking a JSON string as an argument.
         * @param specification A JSON string with the specification for the valve.
         * @details The specification can either be "OutletPressure" or "DeltaPressure" and must be given in pascal.
         */
        impl(const JSONString& specification)
        {
            m_outletStreams.reserve(1);
            setSpecification(specification);
        }

        /**
         * @brief Set the inlet streams to the equipment item.
         * @param streams A std::vector with const Stream* pointers to the input streams.
         * @note The Valve class only supports a single output stream. If more than one is provided,
         * only the first one is considered; the rest are ignored.
         */
        void setInletStreams(const std::vector<const Stream*> streams)
        {
            m_inletStreams = streams;
        }

        /**
         * @brief Sets the equipment specification, given the JSON string with the specification.
         * @param specification specification A JSON string with the specification for the valve.
         * @details The specification can either be "OutletPressure" or "DeltaPressure" and must be given in pascal.
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

            using rapidjson::StringBuffer;
            using rapidjson::Writer;

            StringBuffer         s;
            Writer<StringBuffer> writer(s);
            writer.StartObject();
            writer.Key("FlashSpecification");
            writer.String("PH");
            writer.Key("Pressure");
            writer.Double(m_outletPressure);
            writer.Key("MolarEnthalpy");
            writer.Double(inletFluidProps().mixtureEnthalpy());
            writer.EndObject();

            calcResults();
            m_outletStreams = { *m_inletStreams.front() };
            m_outletStreams.front().setSpecification(s.GetString());
            m_outletStreams.front().compute();
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
            writer.Double(inletFluidProps().front().Pressure);
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
    Valve::Valve()  : m_impl(std::make_unique<impl>("{\"DeltaPressure\": 0.0}")) {}

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