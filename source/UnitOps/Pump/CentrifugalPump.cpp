//
// Created by Kenneth Balslev on 22/01/2021.
//

#include "CentrifugalPump.hpp"

#include <rapidjson/document.h>

namespace PCProps::UnitOps {

    using JSONString = std::string;

    class CentrifugalPump::impl {

        enum class SpecType { OutletPressure, DeltaPressure, HydraulicPower, EffectivePower };

        SpecType m_pumpSpecType;
        double m_pumpSpecValue;
        const Stream* m_inletStream;
        mutable Stream m_outletStream;
        double m_dp;
        double m_pumpEfficiency {};

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

                if (key == "OutletPressure") m_pumpSpecType = SpecType::OutletPressure;
                else if (key == "DeltaPressure") m_pumpSpecType = SpecType::DeltaPressure;
                else if (key == "HydraulicPower") m_pumpSpecType = SpecType::HydraulicPower;
                else if (key == "EffectivePower") m_pumpSpecType = SpecType::EffectivePower;

                m_pumpSpecValue = item.value.GetDouble();
            }

        }

        const Stream& operator()() const {

            auto molarFlow = m_inletStream->properties()[0][PCMolarFlow];
            auto molarVol  = m_inletStream->properties()[0][PCMolarVolume];
            auto volFlow = molarFlow * molarVol;

            auto p = m_inletStream->properties()[0][PCPressure] + m_dp;
            auto t = m_inletStream->properties()[0][PCTemperature];
            m_outletStream = *m_inletStream;
            m_outletStream.flashPT(p, t);
            return m_outletStream;
        }

        std::string results() const {
            return std::string{};
        }

        void setInletStream(Stream* stream) {
            m_inletStream = stream;
        }

        void setDifferentialPressure(double dp) {
            m_dp = dp;
        }

        double computeDensity() const {
            auto density = 0.0;


        }

        double computeLiquidHead(double pressure) const {
//            auto density = 0.0;
//            auto phases = m_inletStream->properties();
//            for (auto& phase : phases) {
//                volFlow +=
//            }
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
    const Stream& CentrifugalPump::operator()() const
    {
        return m_impl->operator()();
    }

    /**
     * @details
     */
    std::string CentrifugalPump::results() const
    {
        return m_impl->results();
    }

    /**
     * @details
     */
    void CentrifugalPump::setInletStream(Stream* stream) {
        m_impl->setInletStream(stream);
    }

    /**
     * @details
     */
    void CentrifugalPump::setSpecification(const JSONString& specification) {}

} //  namespace PCProps::UnitOps