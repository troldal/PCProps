//
// Created by Kenneth Balslev on 22/01/2021.
//

#include "CentrifugalPump.hpp"

#include <json/json.hpp>

namespace PCProps::UnitOps {

    class CentrifugalPump::impl {

        const Stream* m_inletStream;
        mutable Stream m_outletStream;
        double m_dp;



    public:



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



    };

    CentrifugalPump::CentrifugalPump() : m_impl(std::make_unique<impl>()) {}

    CentrifugalPump::CentrifugalPump(const CentrifugalPump& other) : m_impl(std::make_unique<impl>(*other.m_impl)) {};

    CentrifugalPump::CentrifugalPump(CentrifugalPump&& other) noexcept = default;

    CentrifugalPump::~CentrifugalPump() = default;

    CentrifugalPump& CentrifugalPump::operator=(const CentrifugalPump& other)
    {
        CentrifugalPump copy = other;
        *this                = std::move(copy);
        return *this;
    }

    CentrifugalPump& CentrifugalPump::operator=(CentrifugalPump&& other) noexcept = default;

    const Stream& CentrifugalPump::operator()() const
    {
        return m_impl->operator()();
    }

    std::string CentrifugalPump::results() const
    {
        return m_impl->results();
    }

    void CentrifugalPump::setDifferentialPressure(double pressure) {
        m_impl->setDifferentialPressure(pressure);
    }

    void CentrifugalPump::setInletStream(Stream* stream) {
        m_impl->setInletStream(stream);
    }

} //  namespace PCProps::UnitOps