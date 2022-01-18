//
// Created by Kenneth Balslev on 20/01/2021.
//

#ifndef PCPROPS_IPROPERTYPACKAGE_HPP
#define PCPROPS_IPROPERTYPACKAGE_HPP

#include <iostream>
#include <memory>
#include <string>

namespace PCProps::UnitOps
{
    class IPropertyPackage
    {
        using JSONString = std::string;
    public:
        /**
         * @brief Default constructor
         */
        IPropertyPackage() : m_fluid() {}

        /**
         * @brief Constructor, taking the target object as an argument.
         * @tparam T The type of the target object (will be auto deducted)
         * @param x The target object
         */
        template<typename T>
        IPropertyPackage(const T& x) : m_fluid { std::make_unique<Model<T>>(x) }
        {}

        /**
         * @brief
         * @param other
         */
        IPropertyPackage(const IPropertyPackage& other) : m_fluid(other.m_fluid ? other.m_fluid->clone() : nullptr) {}

        /**
         * @brief
         * @param other
         */
        IPropertyPackage(IPropertyPackage&& other) noexcept {}

        /**
         * @brief
         */
        ~IPropertyPackage() {}

        /**
         * @brief
         * @tparam T
         * @param x
         * @return
         */
        template<typename T>
        inline IPropertyPackage& operator=(const T& x)
        {
            std::cout << "IPropertyPackage (template) copy assignment operator called" << std::endl;
            m_fluid = std::make_unique<Model<T>>(x);
            return *this;
        }

        /**
         * @brief
         * @param other
         * @return
         */
        inline IPropertyPackage& operator=(const IPropertyPackage& other)
        {
            std::cout << "IPropertyPackage copy assignment operator called" << std::endl;
            IPropertyPackage copy(other);
            *this = std::move(copy);
            return *this;
        }

        /**
         * @brief
         * @param other
         * @return
         */
        inline IPropertyPackage& operator=(IPropertyPackage&& other) noexcept {
            std::cout << "IPropertyPackage move assignment operator called" << std::endl;
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
        struct Concept
        {
        public:
            /**
             * @brief
             */
            Concept() = default;

            /**
             * @brief
             */
            Concept(const Concept&) = default;

            /**
             * @brief
             */
            Concept(Concept&&) noexcept = default;

            /**
             * @brief
             */
            virtual ~Concept() = default;

            /**
             * @brief
             * @return
             */
            inline Concept& operator=(const Concept&) = default;

            /**
             * @brief
             * @return
             */
            inline Concept& operator=(Concept&&) noexcept = default;

            /**
             * @brief
             * @return
             */
            inline virtual std::unique_ptr<Concept> clone() const = 0;

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
        struct Model : Concept
        {
        public:
            /**
             * @brief
             * @param x
             */
            explicit Model(const T& x) : FluidType(x) {}

            /**
             * @brief
             * @param other
             */
            Model(const Model& other) = default;

            /**
             * @brief
             * @param other
             */
            Model(Model&& other) noexcept = default;

            /**
             * @brief
             */
            ~Model() override = default;

            /**
             * @brief
             * @param other
             * @return
             */
            inline Model& operator=(const Model& other) = default;

            /**
             * @brief
             * @param other
             * @return
             */
            inline Model& operator=(Model&& other) noexcept = default;

            /**
             * @brief
             * @return
             */
            inline std::unique_ptr<Concept> clone() const override
            {
                return std::make_unique<Model<T>>(FluidType);
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

        std::unique_ptr<Concept> m_fluid;

    };
}    // namespace PCProps::UnitOps

#endif    // PCPROPS_IPROPERTYPACKAGE_HPP
