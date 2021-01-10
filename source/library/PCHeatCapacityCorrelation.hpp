//
// Created by Kenneth Balslev on 05/01/2021.
//

#ifndef PCPROPS_PCHEATCAPACITYCORRELATION_HPP
#define PCPROPS_PCHEATCAPACITYCORRELATION_HPP

namespace PCProps
{
    class PCHeatCapacityCorrelation
    {
    public:
        /**
         * @brief Default constructor
         */
        PCHeatCapacityCorrelation() : m_heatCapacityObject() {};

        /**
         * @brief Constructor, taking the target object as an argument.
         * @tparam T The type of the target object (will be auto deducted)
         * @param x The target object
         */
        template<typename T>
        explicit PCHeatCapacityCorrelation(const T& x) : m_heatCapacityObject { std::make_unique<CPModel<T>>(x) }
        {}

        /**
         * @brief
         * @param other
         */
        PCHeatCapacityCorrelation(const PCHeatCapacityCorrelation& other) : m_heatCapacityObject(other.m_heatCapacityObject ? other.m_heatCapacityObject->clone() : nullptr) {}

        /**
         * @brief
         * @param other
         */
        PCHeatCapacityCorrelation(PCHeatCapacityCorrelation&& other) noexcept = default;

        /**
         * @brief
         */
        ~PCHeatCapacityCorrelation() = default;

        /**
         * @brief
         * @tparam T
         * @param x
         * @return
         */
        template<typename T>
        PCHeatCapacityCorrelation& operator=(const T& x)
        {
            m_heatCapacityObject = std::make_unique<CPModel<T>>(x);
            return *this;
        }

        /**
         * @brief
         * @param other
         * @return
         */
        PCHeatCapacityCorrelation& operator=(const PCHeatCapacityCorrelation& other)
        {
            PCHeatCapacityCorrelation copy(other);
            *this = std::move(copy);
            return *this;
        }

        /**
         * @brief
         * @param other
         * @return
         */
        PCHeatCapacityCorrelation& operator=(PCHeatCapacityCorrelation&& other) noexcept = default;

        explicit operator bool() const
        {
            return m_heatCapacityObject != nullptr;
        }

        double evaluateCp(double temperature) const
        {
            return m_heatCapacityObject->evaluateCp(temperature);
        }

        double derivativeOfCp(double temperature) const
        {
            return m_heatCapacityObject->derivativeOfCp(temperature);
        }

        double integralOfCp(double temperature) const
        {
            return m_heatCapacityObject->integralOfCp(temperature);
        }

        double integralOfCpOverT(double temperature) const
        {
            return m_heatCapacityObject->integralOfCpOverT(temperature);
        }

    private:
        /**
         * @brief
         */
        struct CPConcept
        {
        public:
            /**
             * @brief
             */
            CPConcept() = default;

            /**
             * @brief
             */
            CPConcept(const CPConcept&) = default;

            /**
             * @brief
             */
            CPConcept(CPConcept&&) noexcept = default;

            /**
             * @brief
             */
            virtual ~CPConcept() = default;

            /**
             * @brief
             * @return
             */
            CPConcept& operator=(const CPConcept&) = default;

            /**
             * @brief
             * @return
             */
            CPConcept& operator=(CPConcept&&) noexcept = default;

            /**
             * @brief
             * @return
             */
            virtual std::unique_ptr<CPConcept> clone() const = 0;

            virtual double evaluateCp(double temperature) const = 0;

            virtual double derivativeOfCp(double temperature) const = 0;

            virtual double integralOfCp(double temperature) const = 0;

            virtual double integralOfCpOverT(double temperature) const = 0;
        };

        /**
         * @brief
         * @tparam T
         */
        template<typename T>
        struct CPModel : CPConcept
        {
        public:
            /**
             * @brief
             * @param x
             */
            explicit CPModel(const T& x) : CPType(x) {}

            /**
             * @brief
             * @param other
             */
            CPModel(const CPModel& other) = default;

            /**
             * @brief
             * @param other
             */
            CPModel(CPModel&& other) noexcept = default;

            /**
             * @brief
             */
            ~CPModel() override = default;

            /**
             * @brief
             * @param other
             * @return
             */
            CPModel& operator=(const CPModel& other) = default;

            /**
             * @brief
             * @param other
             * @return
             */
            CPModel& operator=(CPModel&& other) noexcept = default;

            /**
             * @brief
             * @return
             */
            std::unique_ptr<CPConcept> clone() const override
            {
                return std::make_unique<CPModel<T>>(CPType);
            }

            double evaluateCp(double temperature) const override
            {
                return CPType.evaluateCp(temperature);
            }

            double derivativeOfCp(double temperature) const override
            {
                return CPType.derivativeOfCp(temperature);
            }

            double integralOfCp(double temperature) const override
            {
                return CPType.integralOfCp(temperature);
            }

            double integralOfCpOverT(double temperature) const override
            {
                return CPType.integralOfCpOverT(temperature);
            }

        private:
            T CPType;
        };

        std::unique_ptr<CPConcept> m_heatCapacityObject {};
    };
}    // namespace PCProps

#endif    // PCPROPS_PCHEATCAPACITYCORRELATION_HPP