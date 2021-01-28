//
// Created by Kenneth Balslev on 20/01/2021.
//

#ifndef PCPROPS_IFLUID_HPP
#define PCPROPS_IFLUID_HPP

#include <iostream>

#include <common/PropertyData.hpp>

namespace PCProps::UnitOps
{
    class IFluid
    {

    public:
        /**
         * @brief Default constructor
         */
        IFluid() : m_fluid() {}

        /**
         * @brief Constructor, taking the target object as an argument.
         * @tparam T The type of the target object (will be auto deducted)
         * @param x The target object
         */
        template<typename T>
        IFluid(const T& x) : m_fluid { std::make_unique<FluidModel<T>>(x) }
        {}

        /**
         * @brief
         * @param other
         */
        IFluid(const IFluid& other) : m_fluid(other.m_fluid ? other.m_fluid->clone() : nullptr) {}

        /**
         * @brief
         * @param other
         */
        IFluid(IFluid&& other) noexcept {}

        /**
         * @brief
         */
        ~IFluid() {}

        /**
         * @brief
         * @tparam T
         * @param x
         * @return
         */
        template<typename T>
        IFluid& operator=(const T& x)
        {
            std::cout << "IFluid (template) copy assignment operator called" << std::endl;
            m_fluid = std::make_unique<FluidModel<T>>(x);
            return *this;
        }

        /**
         * @brief
         * @param other
         * @return
         */
        IFluid& operator=(const IFluid& other)
        {
            std::cout << "IFluid copy assignment operator called" << std::endl;
            IFluid copy(other);
            *this = std::move(copy);
            return *this;
        }

        /**
         * @brief
         * @param other
         * @return
         */
        IFluid& operator=(IFluid&& other) noexcept {
            std::cout << "IFluid move assignment operator called" << std::endl;
            m_fluid = std::move(other.m_fluid);
            return *this;
        }

        /**
         * @brief
         * @return
         */
        explicit operator bool() const
        {
            return m_fluid != nullptr;
        }

        const PCPhases& flashPT(double pressure, double temperature) const {
            return m_fluid->flashPT(pressure, temperature);
        }

        const PCPhases& flashPx(double pressure, double vaporFraction) const {
            return m_fluid->flashPx(pressure, vaporFraction);
        }

        const PCPhases& flashTx(double temperature, double vaporFraction) const {
            return m_fluid->flashTx(temperature, vaporFraction);
        }

        const PCPhases& flashPH(double pressure, double enthalpy) const {
            return m_fluid->flashPH(pressure, enthalpy);
        }

        const PCPhases& flashPS(double pressure, double entropy) const {
            return m_fluid->flashPS(pressure, entropy);
        }

        const PCPhases& flashTV(double temperature, double volume) const {
            return m_fluid->flashPS(temperature, volume);
        }

        const PCPhases& properties() const {
            return m_fluid->properties();
        }

    private:
        /**
         * @brief
         */
        struct FluidConcept
        {
        public:
            /**
             * @brief
             */
            FluidConcept() = default;

            /**
             * @brief
             */
            FluidConcept(const FluidConcept&) = default;

            /**
             * @brief
             */
            FluidConcept(FluidConcept&&) noexcept = default;

            /**
             * @brief
             */
            virtual ~FluidConcept() = default;

            /**
             * @brief
             * @return
             */
            FluidConcept& operator=(const FluidConcept&) = default;

            /**
             * @brief
             * @return
             */
            FluidConcept& operator=(FluidConcept&&) noexcept = default;

            /**
             * @brief
             * @return
             */
            virtual std::unique_ptr<FluidConcept> clone() const = 0;

            /**
             * @brief
             * @param pressure
             * @param temperature
             * @return
             */
            virtual const PCPhases& flashPT(double pressure, double temperature) const = 0;

            /**
             * @brief
             * @param pressure
             * @param vaporFraction
             * @return
             */
            virtual const PCPhases& flashPx(double pressure, double vaporFraction) const = 0;

            /**
             * @brief
             * @param temperature
             * @param vaporFraction
             * @return
             */
            virtual const PCPhases& flashTx(double temperature, double vaporFraction) const = 0;

            /**
             * @brief
             * @param pressure
             * @param enthalpy
             * @return
             */
            virtual const PCPhases& flashPH(double pressure, double enthalpy) const = 0;

            /**
             * @brief
             * @param pressure
             * @param entropy
             * @return
             */
            virtual const PCPhases& flashPS(double pressure, double entropy) const = 0;

            /**
             * @brief
             * @param temperature
             * @param volume
             * @return
             */
            virtual const PCPhases& flashTV(double temperature, double volume) const = 0;

            virtual const PCPhases& properties() const = 0;

        };

        /**
         * @brief
         * @tparam T
         */
        template<typename T>
        struct FluidModel : FluidConcept
        {
        public:
            /**
             * @brief
             * @param x
             */
            explicit FluidModel(const T& x) : FluidType(x) {}

            /**
             * @brief
             * @param other
             */
            FluidModel(const FluidModel& other) = default;

            /**
             * @brief
             * @param other
             */
            FluidModel(FluidModel&& other) noexcept = default;

            /**
             * @brief
             */
            ~FluidModel() override = default;

            /**
             * @brief
             * @param other
             * @return
             */
            FluidModel& operator=(const FluidModel& other) = default;

            /**
             * @brief
             * @param other
             * @return
             */
            FluidModel& operator=(FluidModel&& other) noexcept = default;

            /**
             * @brief
             * @return
             */
            std::unique_ptr<FluidConcept> clone() const override
            {
                return std::make_unique<FluidModel<T>>(FluidType);
            }

            const PCPhases& flashPT(double pressure, double temperature) const override
            {
                return FluidType.flashPT(pressure, temperature);
            }

            const PCPhases& flashPx(double pressure, double vaporFraction) const override
            {
                return FluidType.flashPx(pressure, vaporFraction);
            }

            const PCPhases& flashTx(double temperature, double vaporFraction) const override
            {
                return FluidType.flashTx(temperature, vaporFraction);
            }

            const PCPhases& flashPH(double pressure, double enthalpy) const override
            {
                return FluidType.flashPH(pressure, enthalpy);
            }

            const PCPhases& flashPS(double pressure, double entropy) const override
            {
                return FluidType.flashPS(pressure, entropy);
            }

            const PCPhases& flashTV(double temperature, double volume) const override
            {
                return FluidType.flashTV(temperature, volume);
            }

            const PCPhases& properties() const override
            {
                return FluidType.properties();
            }

        private:
            T FluidType;
        };

        std::unique_ptr<FluidConcept> m_fluid;

    };
}    // namespace PCProps::UnitOps

#endif    // PCPROPS_IFLUID_HPP
