//
// Created by Kenneth Balslev on 18/01/2021.
//

#ifndef PCPROPS_IPURECOMPONENT_HPP
#define PCPROPS_IPURECOMPONENT_HPP

#include <functional>
#include <memory>

namespace PCProps {

    /**
 * @brief This class functions as a wrapper around any class that provides the necessary functionality for
 * cubic equations of state.
 * @details This class works by applying 'type erasure'. This enables the use of objects of any class, the only
 * requirement being that it provides the right interface. No inheritance from a base class is needed.
 */
    class IPureComponent
    {
    public:
        /**
         * @brief Default constructor
         */
        IPureComponent() : m_pureComponent() {}

        /**
         * @brief Constructor, taking the target object as an argument.
         * @tparam T The type of the target object (will be auto deducted)
         * @param x The target object
         */
        template<typename T>
        IPureComponent(const T& x) : m_pureComponent { std::make_unique<PCModel<T>>(x) }
        {}

        /**
         * @brief
         * @param other
         */
        IPureComponent(const IPureComponent& other) : m_pureComponent(other.m_pureComponent ? other.m_pureComponent->clone() : nullptr) {}

        /**
         * @brief
         * @param other
         */
        IPureComponent(IPureComponent&& other) noexcept = default;

        /**
         * @brief
         */
        ~IPureComponent() = default;

        /**
         * @brief
         * @tparam T
         * @param x
         * @return
         */
        template<typename T>
        IPureComponent& operator=(const T& x)
        {
            m_pureComponent = std::make_unique<PCModel<T>>(x);
            return *this;
        }

        /**
         * @brief
         * @param other
         * @return
         */
        IPureComponent& operator=(const IPureComponent& other)
        {
            IPureComponent copy(other);
            *this = std::move(copy);
            return *this;
        }

        /**
         * @brief
         * @param other
         * @return
         */
        IPureComponent& operator=(IPureComponent&& other) noexcept = default;

        /**
         * @brief
         * @return
         */
        explicit operator bool() const
        {
            return m_pureComponent != nullptr;
        }

        double molarWeight() const
        {
            return m_pureComponent->molarWeight();
        }

        double criticalTemperature() const
        {
            return m_pureComponent->criticalTemperature();
        }

        double criticalPressure() const
        {
            return m_pureComponent->criticalPressure();
        }

        double criticalVolume() const
        {
            return m_pureComponent->criticalVolume();
        }

        double criticalCompressibility() const
        {
            return m_pureComponent->criticalCompressibility();
        }

        double acentricFactor() const
        {
            return m_pureComponent->acentricFactor();
        }

        double dipoleMoment() const
        {
            return m_pureComponent->dipoleMoment();
        }

        double idealGasCp(double temperature) const
        {
            return m_pureComponent->idealGasCp(temperature);
        }

        double satLiquidVolume(double temperature) const
        {
            return m_pureComponent->satLiquidVolume(temperature);
        }

        double satVaporViscosity(double temperature) const
        {
            return m_pureComponent->satVaporViscosity(temperature);
        }

        double satLiquidViscosity(double temperature) const
        {
            return m_pureComponent->satLiquidViscosity(temperature);
        }


    private:
        /**
         * @brief
         */
        struct PCConcept
        {
        public:
            /**
             * @brief
             */
            PCConcept() = default;

            /**
             * @brief
             */
            PCConcept(const PCConcept&) = default;

            /**
             * @brief
             */
            PCConcept(PCConcept&&) noexcept = default;

            /**
             * @brief
             */
            virtual ~PCConcept() = default;

            /**
             * @brief
             * @return
             */
            PCConcept& operator=(const PCConcept&) = default;

            /**
             * @brief
             * @return
             */
            PCConcept& operator=(PCConcept&&) noexcept = default;

            /**
             * @brief
             * @return
             */
            virtual std::unique_ptr<PCConcept> clone() const = 0;

            virtual double molarWeight() const = 0;

            virtual double criticalTemperature() const = 0;

            virtual double criticalPressure() const = 0;

            virtual double criticalVolume() const = 0;

            virtual double criticalCompressibility() const = 0;

            virtual double acentricFactor() const = 0;

            virtual double dipoleMoment() const = 0;



            virtual double idealGasCp(double temperature) const = 0;

            virtual double satLiquidVolume(double temperature) const = 0;

            virtual double satVaporViscosity(double temperature) const = 0;

            virtual double satLiquidViscosity(double temperature) const = 0;

        };

        /**
         * @brief
         * @tparam T
         */
        template<typename T>
        struct PCModel : PCConcept
        {
        public:
            /**
             * @brief
             * @param x
             */
            explicit PCModel(const T& x) : PCType(x) {}

            /**
             * @brief
             * @param other
             */
            PCModel(const PCModel& other) = default;

            /**
             * @brief
             * @param other
             */
            PCModel(PCModel&& other) noexcept = default;

            /**
             * @brief
             */
            ~PCModel() override = default;

            /**
             * @brief
             * @param other
             * @return
             */
            PCModel& operator=(const PCModel& other) = default;

            /**
             * @brief
             * @param other
             * @return
             */
            PCModel& operator=(PCModel&& other) noexcept = default;

            /**
             * @brief
             * @return
             */
            std::unique_ptr<PCConcept> clone() const override
            {
                return std::make_unique<PCModel<T>>(PCType);
            }

            double molarWeight() const override
            {
                return PCType.molarWeight();
            }

            double criticalTemperature() const override
            {
                return PCType.criticalTemperature();
            }

            double criticalPressure() const override
            {
                return PCType.criticalPressure();
            }

            double criticalVolume() const override
            {
                return PCType.criticalVolume();
            }

            double criticalCompressibility() const override
            {
                return PCType.criticalCompressibility();
            }

            double acentricFactor() const override
            {
                return PCType.acentricFactor();
            }

            double dipoleMoment() const override
            {
                return PCType.dipoleMoment();
            }

            double idealGasCp(double temperature) const override
            {
                return PCType.idealGasCp(temperature);
            }

            double satLiquidVolume(double temperature) const override
            {
                return PCType.satLiquidVolume(temperature);
            }

            double satVaporViscosity(double temperature) const override
            {
                return PCType.satVaporViscosity(temperature);
            }

            double satLiquidViscosity(double temperature) const override
            {
                return PCType.satLiquidViscosity(temperature);
            }

        private:
            T PCType;
        };

        std::unique_ptr<PCConcept> m_pureComponent;
    };

}

#endif    // PCPROPS_IPURECOMPONENT_HPP
