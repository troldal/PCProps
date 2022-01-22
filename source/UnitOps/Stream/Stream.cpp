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
        FluidProperties m_fluidProps;

    public:

        /**
         * @brief
         * @param fluid
         * @param quantity
         */
        impl(const IPropertyPackage& fluid, double quantity) : m_fluid{fluid}, m_quantity{quantity} {}

        void setQuantity() {
            for (auto& phase : m_fluidProps)
                phase.MolarFlow = phase.MolarFraction * m_quantity;
        }

        /**
         * @brief
         * @param pressure
         * @param temperature
         * @return
         */
        JSONString flashPT(double pressure, double temperature) {
            m_fluidProps = FluidProperties(m_fluid.flash("PT", pressure, temperature));
            setQuantity();
            return properties();
        }

        /**
         * @brief
         * @param pressure
         * @param vaporFraction
         * @return
         */
        JSONString flashPx(double pressure, double vaporFraction) {
            m_fluidProps = FluidProperties(m_fluid.flash("Px", pressure, vaporFraction));
            setQuantity();
            return properties();
        }

        /**
         * @brief
         * @param temperature
         * @param vaporFraction
         * @return
         */
        JSONString flashTx(double temperature, double vaporFraction) {
            m_fluidProps = FluidProperties(m_fluid.flash("Tx", temperature, vaporFraction));
            setQuantity();
            return properties();
        }

        /**
         * @brief
         * @param pressure
         * @param enthalpy
         * @return
         */
        JSONString flashPH(double pressure, double enthalpy) {
            m_fluidProps = FluidProperties(m_fluid.flash("PH", pressure, enthalpy));
            setQuantity();
            return properties();
        }

        /**
         * @brief
         * @param pressure
         * @param entropy
         * @return
         */
        JSONString flashPS(double pressure, double entropy) {
            m_fluidProps = FluidProperties(m_fluid.flash("PS", pressure, entropy));
            setQuantity();
            return properties();
        }

        /**
         * @brief
         * @param temperature
         * @param volume
         * @return
         */
        JSONString flashTV(double temperature, double volume) {
            m_fluidProps = FluidProperties(m_fluid.flash("TV", temperature, volume));
            setQuantity();
            return properties();
        }

        /**
         * @brief
         * @return
         */
        JSONString properties() const {
            return m_fluidProps.asJSON();
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
    JSONString Stream::flash(const std::string& spec, double s1, double s2)
    {
        if (spec == "PT")
            return m_impl->flashPT(s1, s2);
        else if (spec == "Px")
            return m_impl->flashPx(s1, s2);
        else if (spec == "Tx")
            return m_impl->flashTx(s1, s2);
        else if (spec == "PH")
            return m_impl->flashPH(s1, s2);
        else if (spec == "PS")
            return m_impl->flashPS(s1, s2);
        else if (spec == "TV")
            return m_impl->flashTV(s1, s2);
    }

    /**
     * @details
     */
    JSONString Stream::properties() const
    {
        return m_impl->properties();
    }

} // namespace PCProps::UnitOps