//
// Created by Kenneth Balslev on 20/01/2021.
//

#ifndef PCPROPS_IFLUID_HPP
#define PCPROPS_IFLUID_HPP

#include <iostream>
#include <string>

namespace PCProps::UnitOps
{
    class IFluid
    {
        using JSONString = std::string;
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
        inline IFluid& operator=(const T& x)
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
        inline IFluid& operator=(const IFluid& other)
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
        inline IFluid& operator=(IFluid&& other) noexcept {
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

        inline JSONString flashPT(double pressure, double temperature) const {
            return m_fluid->flashPT(pressure, temperature);
        }

        inline JSONString flashPx(double pressure, double vaporFraction) const {
            return m_fluid->flashPx(pressure, vaporFraction);
        }

        inline JSONString flashTx(double temperature, double vaporFraction) const {
            return m_fluid->flashTx(temperature, vaporFraction);
        }

        inline JSONString flashPH(double pressure, double enthalpy) const {
            return m_fluid->flashPH(pressure, enthalpy);
        }

        inline JSONString flashPS(double pressure, double entropy) const {
            return m_fluid->flashPS(pressure, entropy);
        }

        inline JSONString flashTV(double temperature, double volume) const {
            return m_fluid->flashPS(temperature, volume);
        }

        inline JSONString properties() const {
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
            inline FluidConcept& operator=(const FluidConcept&) = default;

            /**
             * @brief
             * @return
             */
            inline FluidConcept& operator=(FluidConcept&&) noexcept = default;

            /**
             * @brief
             * @return
             */
            inline virtual std::unique_ptr<FluidConcept> clone() const = 0;

            /**
             * @brief
             * @param pressure
             * @param temperature
             * @return
             */
            inline virtual JSONString flashPT(double pressure, double temperature) const = 0;

            /**
             * @brief
             * @param pressure
             * @param vaporFraction
             * @return
             */
            inline virtual JSONString flashPx(double pressure, double vaporFraction) const = 0;

            /**
             * @brief
             * @param temperature
             * @param vaporFraction
             * @return
             */
            inline virtual JSONString flashTx(double temperature, double vaporFraction) const = 0;

            /**
             * @brief
             * @param pressure
             * @param enthalpy
             * @return
             */
            inline virtual JSONString flashPH(double pressure, double enthalpy) const = 0;

            /**
             * @brief
             * @param pressure
             * @param entropy
             * @return
             */
            inline virtual JSONString flashPS(double pressure, double entropy) const = 0;

            /**
             * @brief
             * @param temperature
             * @param volume
             * @return
             */
            inline virtual JSONString flashTV(double temperature, double volume) const = 0;

            inline virtual JSONString properties() const = 0;

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
            inline FluidModel& operator=(const FluidModel& other) = default;

            /**
             * @brief
             * @param other
             * @return
             */
            inline FluidModel& operator=(FluidModel&& other) noexcept = default;

            /**
             * @brief
             * @return
             */
            inline std::unique_ptr<FluidConcept> clone() const override
            {
                return std::make_unique<FluidModel<T>>(FluidType);
            }

            inline JSONString flashPT(double pressure, double temperature) const override
            {
                return FluidType.flashPT(pressure, temperature);
            }

            inline JSONString flashPx(double pressure, double vaporFraction) const override
            {
                return FluidType.flashPx(pressure, vaporFraction);
            }

            inline JSONString flashTx(double temperature, double vaporFraction) const override
            {
                return FluidType.flashTx(temperature, vaporFraction);
            }

            inline JSONString flashPH(double pressure, double enthalpy) const override
            {
                return FluidType.flashPH(pressure, enthalpy);
            }

            inline JSONString flashPS(double pressure, double entropy) const override
            {
                return FluidType.flashPS(pressure, entropy);
            }

            inline JSONString flashTV(double temperature, double volume) const override
            {
                return FluidType.flashTV(temperature, volume);
            }

            inline JSONString properties() const override
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
