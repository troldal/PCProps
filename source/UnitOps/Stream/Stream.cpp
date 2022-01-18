//
// Created by Kenneth Balslev on 24/01/2021.
//

#include "Stream.hpp"

#include <FluidProperties.hpp>

using JSONString = std::string;

namespace PCProps::UnitOps {

    /**
     * @brief
     */
    class Stream::impl {
        IPropertyPackage m_fluid;
        double m_quantity;

    public:

        /**
         * @brief
         * @param fluid
         * @param quantity
         */
        impl(const IPropertyPackage& fluid, double quantity) : m_fluid{fluid}, m_quantity{quantity} {}

        /**
         * @brief
         * @param pressure
         * @param temperature
         * @return
         */
        JSONString flashPT(double pressure, double temperature) const {
            m_fluid.flashPT(pressure, temperature), m_quantity;
            return properties();
        }

        /**
         * @brief
         * @param pressure
         * @param vaporFraction
         * @return
         */
        JSONString flashPx(double pressure, double vaporFraction) const {
            m_fluid.flashPx(pressure, vaporFraction), m_quantity;
            return properties();
        }

        /**
         * @brief
         * @param temperature
         * @param vaporFraction
         * @return
         */
        JSONString flashTx(double temperature, double vaporFraction) const {
            m_fluid.flashTx(temperature, vaporFraction), m_quantity;
            return properties();
        }

        /**
         * @brief
         * @param pressure
         * @param enthalpy
         * @return
         */
        JSONString flashPH(double pressure, double enthalpy) const {
            m_fluid.flashPH(pressure, enthalpy), m_quantity;
            return properties();
        }

        /**
         * @brief
         * @param pressure
         * @param entropy
         * @return
         */
        JSONString flashPS(double pressure, double entropy) const {
            m_fluid.flashPS(pressure, entropy), m_quantity;
            return properties();
        }

        /**
         * @brief
         * @param temperature
         * @param volume
         * @return
         */
        JSONString flashTV(double temperature, double volume) const {
            m_fluid.flashPS(temperature, volume), m_quantity;
            return properties();
        }

        /**
         * @brief
         * @return
         */
        JSONString properties() const {
            auto fluid = FluidProperties(m_fluid.properties());
            for (auto& phase : fluid)
                phase.MolarFlow *= m_quantity;

            return fluid.asJSON();
//            auto phases = m_fluid.properties();
//            for(auto& item : phases) item[PCMolarFlow] *= m_quantity;
//            return phases;
        }
    };

    /**
     * @details
     */
    Stream::Stream() = default;

    /**
     * @details
     */
    Stream::Stream(const IPropertyPackage& fluid, double quantity) : m_impl(std::make_unique<impl>(fluid, quantity)) {}

    /**
     * @details
     */
    Stream::Stream(const Stream& other) : m_impl(std::make_unique<impl>(*other.m_impl)) {};

    /**
     * @details
     */
    Stream::Stream(Stream&& other) noexcept = default;

    /**
     * @details
     */
    Stream::~Stream() = default;

    /**
     * @details
     */
    Stream& Stream::operator=(const Stream& other)
    {
        Stream copy = other;
        *this       = std::move(copy);
        return *this;
    }

    /**
     * @details
     */
    Stream& Stream::operator=(Stream&& other) noexcept = default;

    /**
     * @details
     */
    JSONString Stream::flashPT(double pressure, double temperature) const
    {
        return m_impl->flashPT(pressure, temperature);
    }

    /**
     * @details
     */
    JSONString Stream::flashPx(double pressure, double vaporFraction) const
    {
        return m_impl->flashPx(pressure, vaporFraction);
    }

    /**
     * @details
     */
    JSONString Stream::flashTx(double temperature, double vaporFraction) const
    {
        return m_impl->flashTx(temperature, vaporFraction);
    }

    /**
     * @details
     */
    JSONString Stream::flashPH(double pressure, double enthalpy) const
    {
        return m_impl->flashPH(pressure, enthalpy);
    }

    /**
     * @details
     */
    JSONString Stream::flashPS(double pressure, double entropy) const
    {
        return m_impl->flashPS(pressure, entropy);
    }

    /**
     * @details
     */
    JSONString Stream::flashTV(double temperature, double volume) const
    {
        return m_impl->flashTV(temperature, volume);
    }

    /**
     * @details
     */
    JSONString Stream::properties() const
    {
        return m_impl->properties();
    }



} // namespace PCProps::UnitOps