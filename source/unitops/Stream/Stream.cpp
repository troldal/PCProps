//
// Created by Kenneth Balslev on 24/01/2021.
//

#include "Stream.hpp"

namespace PCProps::UnitOps {

    class Stream::impl {

        IFluid m_fluid;
        double m_quantity;

    public:

        impl(const IFluid& fluid, double quantity) : m_fluid{fluid}, m_quantity{quantity} {}

        PCPhases flashPT(double pressure, double temperature) const {
            m_fluid.flashPT(pressure, temperature), m_quantity;
            return properties();
        }

        PCPhases flashPx(double pressure, double vaporFraction) const {
            m_fluid.flashPx(pressure, vaporFraction), m_quantity;
            return properties();
        }

        PCPhases flashTx(double temperature, double vaporFraction) const {
            m_fluid.flashTx(temperature, vaporFraction), m_quantity;
            return properties();
        }

        PCPhases flashPH(double pressure, double enthalpy) const {
            m_fluid.flashPH(pressure, enthalpy), m_quantity;
            return properties();
        }

        PCPhases flashPS(double pressure, double entropy) const {
            m_fluid.flashPS(pressure, entropy), m_quantity;
            return properties();
        }

        PCPhases flashTV(double temperature, double volume) const {
            m_fluid.flashPS(temperature, volume), m_quantity;
            return properties();
        }

        PCPhases properties() const {
            auto phases = m_fluid.properties();
            for(auto& item : phases) item[PCMolarFlow] *= m_quantity;
            return phases;
        }
    };

    Stream::Stream() = default;

    Stream::Stream(const IFluid& fluid, double quantity) : m_impl(std::make_unique<impl>(fluid, quantity)) {}

    Stream::Stream(const Stream& other) : m_impl(std::make_unique<impl>(*other.m_impl)) {};

    Stream::Stream(Stream&& other) noexcept = default;

    Stream::~Stream() = default;

    Stream& Stream::operator=(const Stream& other)
    {
        Stream copy = other;
        *this       = std::move(copy);
        return *this;
    }

    Stream& Stream::operator=(Stream&& other) noexcept = default;

    PCPhases Stream::flashPT(double pressure, double temperature) const
    {
        return m_impl->flashPT(pressure, temperature);
    }

    PCPhases Stream::flashPx(double pressure, double vaporFraction) const
    {
        return m_impl->flashPx(pressure, vaporFraction);
    }

    PCPhases Stream::flashTx(double temperature, double vaporFraction) const
    {
        return m_impl->flashTx(temperature, vaporFraction);
    }

    PCPhases Stream::flashPH(double pressure, double enthalpy) const
    {
        return m_impl->flashPH(pressure, enthalpy);
    }

    PCPhases Stream::flashPS(double pressure, double entropy) const
    {
        return m_impl->flashPS(pressure, entropy);
    }

    PCPhases Stream::flashTV(double temperature, double volume) const
    {
        return m_impl->flashTV(temperature, volume);
    }

    PCPhases Stream::properties() const
    {
        return m_impl->properties();
    }

} // namespace PCProps::UnitOps