//
// Created by Kenneth Balslev on 24/01/2021.
//

#ifndef PCPROPS_STREAM_HPP
#define PCPROPS_STREAM_HPP

#include <Interfaces/IFluid.hpp>

#include <memory>

namespace PCProps::UnitOps
{
    class Stream
    {
    public:

        Stream();

        Stream(const IFluid& fluid, double quantity);

        Stream(const Stream& other);

        Stream(Stream&& other) noexcept;

        ~Stream();

        Stream& operator=(const Stream& other);

        Stream& operator=(Stream&& other) noexcept;

        PCPhases flashPT(double pressure, double temperature) const;

        PCPhases flashPx(double pressure, double vaporFraction) const;

        PCPhases flashTx(double temperature, double vaporFraction) const;

        PCPhases flashPH(double pressure, double enthalpy) const;

        PCPhases flashPS(double pressure, double entropy) const;

        PCPhases flashTV(double temperature, double volume) const;

        PCPhases properties() const;

    private:
        class impl;
        std::unique_ptr<impl> m_impl;

    };
}    // namespace PCProps::UnitOps

#endif    // PCPROPS_STREAM_HPP
