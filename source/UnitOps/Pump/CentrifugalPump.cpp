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
        double m_eff;



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

        impl(double eff) : m_eff(eff) {}

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

    CentrifugalPump::CentrifugalPump(double eff) : m_impl(std::make_unique<impl>(eff)) {}

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