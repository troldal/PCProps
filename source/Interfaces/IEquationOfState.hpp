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

#ifndef PCPROPS_IEQUATIONOFSTATE_HPP
#define PCPROPS_IEQUATIONOFSTATE_HPP

#include <functional>
#include <memory>

#include <common/PropertyData.hpp>
#include <IPureComponent.hpp>

namespace PCProps
{
    /**
     * @brief This class functions as a wrapper around any class that provides the necessary functionality for
     * equations of state.
     * @details This class works by applying 'type erasure'. This enables the use of objects of any class, the only
     * requirement being that it provides the right interface. No inheritance from a base class is needed.
     */
    class IEquationOfState
    {
        using JSONString = std::string;
    public:
        /**
         * @brief Default constructor
         */
        IEquationOfState() : m_equationOfState() {} // NOLINT

        /**
         * @brief Constructor, taking the target object as an argument.
         * @tparam T The type of the target object (will be auto deducted)
         * @param x The target object
         * @note This methed is deliberately not marked 'explicit', because as a templated constructor, it should be able
         * to take any type as an argument. However, only objects that satisfy the required interface can be used.
         */
        template<typename T>
        IEquationOfState(const T& x) : m_equationOfState { std::make_unique<EOSModel<T>>(x) } {} // NOLINT

        /**
         * @brief Copy constructor
         * @param other
         */
        IEquationOfState(const IEquationOfState& other) : m_equationOfState(other.m_equationOfState ? other.m_equationOfState->clone() : nullptr) {}

        /**
         * @brief Move constructor
         * @param other
         */
        IEquationOfState(IEquationOfState&& other) noexcept = default;

        /**
         * @brief Destructor
         */
        ~IEquationOfState() = default;

        /**
         * @brief
         * @tparam T
         * @param x
         * @return
         */
        template<typename T>
        inline IEquationOfState& operator=(const T& x)
        {
            m_equationOfState = std::make_unique<EOSModel<T>>(x);
            return *this;
        }

        /**
         * @brief
         * @param other
         * @return
         */
        inline IEquationOfState& operator=(const IEquationOfState& other)
        {
            IEquationOfState copy(other);
            *this = std::move(copy);
            return *this;
        }

        /**
         * @brief
         * @param other
         * @return
         */
        inline IEquationOfState& operator=(IEquationOfState&& other) noexcept = default;

        /**
         * @brief
         * @return
         */
        inline explicit operator bool() const
        {
            return m_equationOfState != nullptr;
        }

        /**
         * @brief
         * @param jsonString
         */
        inline void init(const std::tuple<double, double, double, std::function<double(double)> >& properties) {
            m_equationOfState->init(properties);
        }

        /**
         * @brief
         * @param pressure
         * @param temperature
         * @return
         */
        inline JSONString flashPT(double pressure, double temperature) const
        {
            return m_equationOfState->flashPT(pressure, temperature);
        }

        /**
         * @brief
         * @param pressure
         * @param vaporFraction
         * @return
         */
        inline JSONString flashPx(double pressure, double vaporFraction) const
        {
            return m_equationOfState->flashPx(pressure, vaporFraction);
        }

        /**
         * @brief
         * @param temperature
         * @param vaporFraction
         * @return
         */
        inline JSONString flashTx(double temperature, double vaporFraction) const
        {
            return m_equationOfState->flashTx(temperature, vaporFraction);
        }

        /**
         * @brief
         * @param pressure
         * @param enthalpy
         * @return
         */
        inline JSONString flashPH(double pressure, double enthalpy) const
        {
            return m_equationOfState->flashPH(pressure, enthalpy);
        }

        /**
         * @brief
         * @param pressure
         * @param entropy
         * @return
         */
        inline JSONString flashPS(double pressure, double entropy) const
        {
            return m_equationOfState->flashPS(pressure, entropy);
        }

        /**
         * @brief
         * @param temperature
         * @param volume
         * @return
         */
        inline JSONString flashTV(double temperature, double volume) const
        {
            return m_equationOfState->flashTV(temperature, volume);
        }

        /**
         * @brief
         * @param temperature
         * @return
         */
        inline double saturationPressure(double temperature) const
        {
            return m_equationOfState->saturationPressure(temperature);
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
            inline EOSConcept& operator=(const EOSConcept&) = default;

            /**
             * @brief
             * @return
             */
            inline EOSConcept& operator=(EOSConcept&&) noexcept = default;

            /**
             * @brief
             * @return
             */
            inline virtual std::unique_ptr<EOSConcept> clone() const = 0;

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

            /**
             * @brief
             * @param temperature
             * @return
             */
            inline virtual double saturationPressure(double temperature) const = 0;

            /**
             * @brief
             * @param jsonString
             */
            inline virtual void init(const std::tuple<double, double, double, std::function<double(double)> >& properties) = 0;

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
            inline EOSModel& operator=(const EOSModel& other) = default;

            /**
             * @brief
             * @param other
             * @return
             */
            inline EOSModel& operator=(EOSModel&& other) noexcept = default;

            /**
             * @brief
             * @return
             */
            inline std::unique_ptr<EOSConcept> clone() const override
            {
                return std::make_unique<EOSModel<T>>(EOSType);
            }

            /**
             * @brief
             * @param pressure
             * @param temperature
             * @return
             */
            inline JSONString flashPT(double pressure, double temperature) const override
            {
                return EOSType.flashPT(pressure, temperature);
            }

            /**
             * @brief
             * @param pressure
             * @param vaporFraction
             * @return
             */
            inline JSONString flashPx(double pressure, double vaporFraction) const override
            {
                return EOSType.flashPx(pressure, vaporFraction);
            }

            /**
             * @brief
             * @param temperature
             * @param vaporFraction
             * @return
             */
            inline JSONString flashTx(double temperature, double vaporFraction) const override
            {
                return EOSType.flashTx(temperature, vaporFraction);
            }

            /**
             * @brief
             * @param pressure
             * @param enthalpy
             * @return
             */
            inline JSONString flashPH(double pressure, double enthalpy) const override
            {
                return EOSType.flashPH(pressure, enthalpy);
            }

            /**
             * @brief
             * @param pressure
             * @param entropy
             * @return
             */
            inline JSONString flashPS(double pressure, double entropy) const override
            {
                return EOSType.flashPS(pressure, entropy);
            }

            /**
             * @brief
             * @param temperature
             * @param volume
             * @return
             */
            inline JSONString flashTV(double temperature, double volume) const override
            {
                return EOSType.flashTV(temperature, volume);
            }

            /**
             * @brief
             * @param temperature
             * @return
             */
            inline double saturationPressure(double temperature) const override
            {
                return EOSType.saturationPressure(temperature);
            }

            /**
             * @brief
             * @param jsonString
             */
            inline void init(const std::tuple<double, double, double, std::function<double(double)> >& properties) override {
                EOSType.init(properties);
            }

        private:
            T EOSType;
        };

        std::unique_ptr<EOSConcept> m_equationOfState;
    };

}    // namespace PCProps

#endif    // PCPROPS_IEQUATIONOFSTATE_HPP
