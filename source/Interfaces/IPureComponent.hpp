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
        IPureComponent() : m_pureComponent() {} // NOLINT

        /**
         * @brief Constructor, taking the target object as an argument.
         * @tparam T The type of the target object (will be auto deducted)
         * @param x The target object
         */
        template<typename T>
        IPureComponent(const T& x) : m_pureComponent { std::make_unique<PCModel<T>>(x) } // NOLINT
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

        /**
         *
         * @param ID
         * @return
         */
        double property(const std::string& ID) const {
            return m_pureComponent->property(ID);
        }

        /**
         *
         * @param ID
         * @param temperature
         * @return
         */
        double property(const std::string& ID, double temperature) const {
            return m_pureComponent->property(ID, temperature);
        }

        /**
         *
         * @param ID
         * @param parameters
         * @return
         */
        double property(const std::string& ID, const std::vector<double>& parameters) const {
            return m_pureComponent->property(ID, parameters);
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

            /**
             *
             * @param ID
             * @return
             */
            virtual double property(const std::string& ID) const = 0;

            /**
             *
             * @param ID
             * @param temperature
             * @return
             */
            virtual double property(const std::string& ID, double temperature) const = 0;

            /**
             *
             * @param ID
             * @param parameters
             * @return
             */
            virtual double property(const std::string& ID, const std::vector<double>& parameters) const = 0;

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

            /**
             *
             * @param ID
             * @return
             */
            double property(const std::string& ID) const override {
                return PCType.property(ID);
            }

            /**
             *
             * @param ID
             * @param temperature
             * @return
             */
            double property(const std::string& ID, double temperature) const override {
                return PCType.property(ID, temperature);
            }

            /**
             *
             * @param ID
             * @param parameters
             * @return
             */
            double property(const std::string& ID, const std::vector<double>& parameters) const override {
                return PCType.property(ID, parameters);
            }

        private:
            T PCType; /**< */
        };

        std::unique_ptr<PCConcept> m_pureComponent; /**< */
    };

} // namespace PCProps

#endif    // PCPROPS_IPURECOMPONENT_HPP
