/*

8888888b.   .d8888b.  8888888b.
888   Y88b d88P  Y88b 888   Y88b
888    888 888    888 888    888
888   d88P 888        888   d88P 888d888 .d88b.  88888b.  .d8888b
8888888P"  888        8888888P"  888P"  d88""88b 888 "88b 88K
888        888    888 888        888    888  888 888  888 "Y8888b.
888        Y88b  d88P 888        888    Y88..88P 888 d88P      X88
888         "Y8888P"  888        888     "Y88P"  88888P"   88888P'
                                                 888
                                                 888
                                                 888

Copyright (c) 2020 Kenneth Troldal Balslev

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#ifndef PCPROPS_PCEQUATIONOFSTATE_HPP
#define PCPROPS_PCEQUATIONOFSTATE_HPP

#include <memory>

#include <EquationOfState/EOSUtilities.hpp>

using PCProps::EquationOfState::Phases;

namespace PCProps
{
    /**
     * @brief This class functions as a wrapper around any class that provides the necessary functionality for
     * cubic equations of state.
     * @details This class works by applying 'type erasure'. This enables the use of objects of any class, the only
     * requirement being that it provides the right interface. No inheritance from a base class is needed.
     */
    class PCEquationOfState
    {
    public:
        /**
         * @brief Default constructor
         */
        PCEquationOfState() : m_equationOfState() {}

        /**
         * @brief Constructor, taking the target object as an argument.
         * @tparam T The type of the target object (will be auto deducted)
         * @param x The target object
         */
        template<typename T>
        explicit PCEquationOfState(const T& x) : m_equationOfState { std::make_unique<EOSModel<T>>(x) }
        {}

        /**
         * @brief
         * @param other
         */
        PCEquationOfState(const PCEquationOfState& other) : m_equationOfState(other.m_equationOfState ? other.m_equationOfState->clone() : nullptr) {}

        /**
         * @brief
         * @param other
         */
        PCEquationOfState(PCEquationOfState&& other) noexcept = default;

        /**
         * @brief
         */
        ~PCEquationOfState() = default;

        /**
         * @brief
         * @tparam T
         * @param x
         * @return
         */
        template<typename T>
        PCEquationOfState& operator=(const T& x)
        {
            m_equationOfState = std::make_unique<EOSModel<T>>(x);
            return *this;
        }

        /**
         * @brief
         * @param other
         * @return
         */
        PCEquationOfState& operator=(const PCEquationOfState& other)
        {
            PCEquationOfState copy(other);
            *this = std::move(copy);
            return *this;
        }

        /**
         * @brief
         * @param other
         * @return
         */
        PCEquationOfState& operator=(PCEquationOfState&& other) noexcept = default;

        explicit operator bool() const
        {
            return m_equationOfState != nullptr;
        }

        void setProperties(double criticalTemperature, double criticalPressure, double acentricFactor)
        {
            m_equationOfState->setProperties(criticalTemperature, criticalPressure, acentricFactor);
        }

        void setVaporPressureFunction(const std::function<double(double)>& vaporPressureFunction)
        {
            m_equationOfState->setVaporPressureFunction(vaporPressureFunction);
        }

        void setIdealGasCpFunction(const std::function<double(double)>& idealGasCpFunction)
        {
            m_equationOfState->setIdealGasCpFunction(idealGasCpFunction);
        }

        void setIdealGasCpDerivativeFunction(const std::function<double(double)>& idealGasCpDerivativeFunction)
        {
            m_equationOfState->setIdealGasCpDerivativeFunction(idealGasCpDerivativeFunction);
        }

        void setIdealGasCpIntegralFunction(const std::function<double(double)>& idealGasCpIntegralFunction)
        {
            m_equationOfState->setIdealGasCpIntegralFunction(idealGasCpIntegralFunction);
        }

        void setIdealGasCpOverTIntegralFunction(const std::function<double(double)>& idealGasOverTIntegralFunction)
        {
            m_equationOfState->setIdealGasCpOverTIntegralFunction(idealGasOverTIntegralFunction);
        }

        /**
         * @brief
         * @param pressure
         * @param temperature
         * @return
         */
        Phases flashPT(double pressure, double temperature) const
        {
            return m_equationOfState->flashPT(pressure, temperature);
        }

        /**
         * @brief
         * @param pressure
         * @param vaporFraction
         * @return
         */
        Phases flashPx(double pressure, double vaporFraction) const
        {
            return m_equationOfState->flashPx(pressure, vaporFraction);
        }

        /**
         * @brief
         * @param temperature
         * @param vaporFraction
         * @return
         */
        Phases flashTx(double temperature, double vaporFraction) const
        {
            return m_equationOfState->flashTx(temperature, vaporFraction);
        }

        /**
         * @brief
         * @param pressure
         * @param enthalpy
         * @return
         */
        Phases flashPH(double pressure, double enthalpy) const
        {
            return m_equationOfState->flashPH(pressure, enthalpy);
        }

        /**
         * @brief
         * @param pressure
         * @param entropy
         * @return
         */
        Phases flashPS(double pressure, double entropy) const
        {
            return m_equationOfState->flashPS(pressure, entropy);
        }

        double saturationPressure(double temperature) const
        {
            return m_equationOfState->saturationPressure(temperature);
        }

        double saturationTemperature(double pressure) const
        {
            return m_equationOfState->saturationTemperature(pressure);
        }

    private:
        /**
         * @brief
         */
        struct EOSConcept
        {
        public:
            /**
             * @brief
             */
            EOSConcept() = default;

            /**
             * @brief
             */
            EOSConcept(const EOSConcept&) = default;

            /**
             * @brief
             */
            EOSConcept(EOSConcept&&) noexcept = default;

            /**
             * @brief
             */
            virtual ~EOSConcept() = default;

            /**
             * @brief
             * @return
             */
            EOSConcept& operator=(const EOSConcept&) = default;

            /**
             * @brief
             * @return
             */
            EOSConcept& operator=(EOSConcept&&) noexcept = default;

            /**
             * @brief
             * @return
             */
            virtual std::unique_ptr<EOSConcept> clone() const = 0;

            /**
             * @brief
             * @param pressure
             * @param temperature
             * @return
             */
            virtual Phases flashPT(double pressure, double temperature) const = 0;

            /**
             * @brief
             * @param pressure
             * @param vaporFraction
             * @return
             */
            virtual Phases flashPx(double pressure, double vaporFraction) const = 0;

            /**
             * @brief
             * @param temperature
             * @param vaporFraction
             * @return
             */
            virtual Phases flashTx(double temperature, double vaporFraction) const = 0;

            /**
             * @brief
             * @param pressure
             * @param enthalpy
             * @return
             */
            virtual Phases flashPH(double pressure, double enthalpy) const = 0;

            /**
             * @brief
             * @param pressure
             * @param entropy
             * @return
             */
            virtual Phases flashPS(double pressure, double entropy) const = 0;

            /**
             * @brief
             * @param temperature
             * @return
             */
            virtual double saturationPressure(double temperature) const = 0;

            /**
             * @brief
             * @param pressure
             * @return
             */
            virtual double saturationTemperature(double pressure) const = 0;

            /**
             * @brief
             * @param criticalTemperature
             * @param criticalPressure
             * @param acentricFactor
             */
            virtual void setProperties(double criticalTemperature, double criticalPressure, double acentricFactor) = 0;

            /**
             * @brief
             * @param vaporPressureFunction
             */
            virtual void setVaporPressureFunction(const std::function<double(double)>& vaporPressureFunction) = 0;

            /**
             * @brief
             * @param idealGasCpFunction
             */
            virtual void setIdealGasCpFunction(const std::function<double(double)>& idealGasCpFunction) = 0;

            /**
             * @brief
             * @param idealGasCpDerivativeFunction
             */
            virtual void setIdealGasCpDerivativeFunction(const std::function<double(double)>& idealGasCpDerivativeFunction) = 0;

            /**
             * @brief
             * @param idealGasCpIntegralFunction
             */
            virtual void setIdealGasCpIntegralFunction(const std::function<double(double)>& idealGasCpIntegralFunction) = 0;

            /**
             * @brief
             * @param idealGasOverTIntegralFunction
             */
            virtual void setIdealGasCpOverTIntegralFunction(const std::function<double(double)>& idealGasOverTIntegralFunction) = 0;
        };

        /**
         * @brief
         * @tparam T
         */
        template<typename T>
        struct EOSModel : EOSConcept
        {
        public:
            /**
             * @brief
             * @param x
             */
            explicit EOSModel(const T& x) : EOSType(x) {}

            /**
             * @brief
             * @param other
             */
            EOSModel(const EOSModel& other) = default;

            /**
             * @brief
             * @param other
             */
            EOSModel(EOSModel&& other) noexcept = default;

            /**
             * @brief
             */
            ~EOSModel() override = default;

            /**
             * @brief
             * @param other
             * @return
             */
            EOSModel& operator=(const EOSModel& other) = default;

            /**
             * @brief
             * @param other
             * @return
             */
            EOSModel& operator=(EOSModel&& other) noexcept = default;

            /**
             * @brief
             * @return
             */
            std::unique_ptr<EOSConcept> clone() const override
            {
                return std::make_unique<EOSModel<T>>(EOSType);
            }

            /**
             * @brief
             * @param pressure
             * @param temperature
             * @return
             */
            Phases flashPT(double pressure, double temperature) const override
            {
                return EOSType.flashPT(pressure, temperature);
            }

            /**
             * @brief
             * @param pressure
             * @param vaporFraction
             * @return
             */
            Phases flashPx(double pressure, double vaporFraction) const override
            {
                return EOSType.flashPx(pressure, vaporFraction);
            }

            /**
             * @brief
             * @param temperature
             * @param vaporFraction
             * @return
             */
            Phases flashTx(double temperature, double vaporFraction) const override
            {
                return EOSType.flashTx(temperature, vaporFraction);
            }

            /**
             * @brief
             * @param pressure
             * @param enthalpy
             * @return
             */
            Phases flashPH(double pressure, double enthalpy) const override
            {
                return EOSType.flashPH(pressure, enthalpy);
            }

            /**
             * @brief
             * @param pressure
             * @param entropy
             * @return
             */
            Phases flashPS(double pressure, double entropy) const override
            {
                return EOSType.flashPS(pressure, entropy);
            }

            /**
             * @brief
             * @param temperature
             * @return
             */
            double saturationPressure(double temperature) const override
            {
                return EOSType.saturationPressure(temperature);
            }

            /**
             * @brief
             * @param pressure
             * @return
             */
            double saturationTemperature(double pressure) const override
            {
                return EOSType.saturationTemperature(pressure);
            }

            void setProperties(double criticalTemperature, double criticalPressure, double acentricFactor) override
            {
                EOSType.setProperties(criticalTemperature, criticalPressure, acentricFactor);
            }

            void setVaporPressureFunction(const std::function<double(double)>& vaporPressureFunction) override
            {
                EOSType.setVaporPressureFunction(vaporPressureFunction);
            }

            void setIdealGasCpFunction(const std::function<double(double)>& idealGasCpFunction) override
            {
                EOSType.setIdealGasCpFunction(idealGasCpFunction);
            }

            void setIdealGasCpDerivativeFunction(const std::function<double(double)>& idealGasCpDerivativeFunction) override
            {
                EOSType.setIdealGasCpDerivativeFunction(idealGasCpDerivativeFunction);
            }

            void setIdealGasCpIntegralFunction(const std::function<double(double)>& idealGasCpIntegralFunction) override
            {
                EOSType.setIdealGasCpIntegralFunction(idealGasCpIntegralFunction);
            }

            void setIdealGasCpOverTIntegralFunction(const std::function<double(double)>& idealGasOverTIntegralFunction) override
            {
                EOSType.setIdealGasCpOverTIntegralFunction(idealGasOverTIntegralFunction);
            }

        private:
            T EOSType;
        };

        std::unique_ptr<EOSConcept> m_equationOfState;
    };

}    // namespace PCProps

#endif    // PCPROPS_PCEQUATIONOFSTATE_HPP