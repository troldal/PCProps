//
// Created by Kenneth Balslev on 22/01/2021.
//

#ifndef PCPROPS_CENTRIFUGALPUMP_HPP
#define PCPROPS_CENTRIFUGALPUMP_HPP

#include <IPropertyPackage.hpp>
#include <Stream/Stream.hpp>

#include <array>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace PCProps::UnitOps
{

    class CentrifugalPump
    {
    public:

        CentrifugalPump(double eff = 1.0);

        CentrifugalPump(const CentrifugalPump& other);

        CentrifugalPump(CentrifugalPump&& other) noexcept;

        ~CentrifugalPump();

        CentrifugalPump& operator=(const CentrifugalPump& other);

        CentrifugalPump& operator=(CentrifugalPump&& other) noexcept;

//        CentrifugalPump(const IPropertyPackage& inletFluid);

        // Inlet + Outlet P
        // Inlet + dP
        // Inlet + Power

        const Stream& operator()() const;

        std::string results() const;

        void setDifferentialPressure(double pressure);

        void setInletStream(Stream* stream);


    private:

        class impl;
        std::unique_ptr<impl> m_impl;

    };
}    // namespace PCProps::UnitOps

#endif    // PCPROPS_CENTRIFUGALPUMP_HPP
