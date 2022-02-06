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
#include <string>

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
        IEquationOfState(const T& eos) : m_equationOfState { std::make_unique<Model<T>>(eos) } {} // NOLINT

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
        inline IEquationOfState& operator=(const T& eos)
        {
            m_equationOfState = std::make_unique<Model<T>>(eos);
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
         * @param pureComponent
         */
        inline void init(const IPureComponent& pureComponent)   {
            m_equationOfState->init(pureComponent);
        }

        /**
         * @brief
         * @param pressure
         * @param temperature
         * @return
         */
        inline JSONString computePropertiesPT(double pressure, double temperature) const
        {
            return m_equationOfState->computePropertiesPT(pressure, temperature);
        }

        /**
         * @brief
         * @param temperature
         * @param molarVolume
         * @return
         */
        inline JSONString computePropertiesTV(double temperature, double molarVolume) const
        {
            return m_equationOfState->computePropertiesTV(temperature, molarVolume);
        }

        /**
         * @brief
         * @param temperature
         * @param molarVolume
         * @return
         */
        inline JSONString computePropertiesPV(double pressure, double molarVolume) const
        {
            return m_equationOfState->computePropertiesPV(pressure, molarVolume);
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

        /**
         * @brief
         * @param pressure
         * @return
         */
        inline double saturationTemperature(double pressure) const
        {
            return m_equationOfState->saturationTemperature(pressure);
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
            inline virtual JSONString computePropertiesPT(double pressure, double temperature) const = 0;

            /**
             * @brief
             * @param temperature
             * @param molarVolume
             * @return
             */
            inline virtual JSONString computePropertiesTV(double temperature, double molarVolume) const = 0;

            /**
             * @brief
             * @param temperature
             * @param molarVolume
             * @return
             */
            inline virtual JSONString computePropertiesPV(double pressure, double molarVolume) const = 0;

            /**
             * @brief
             * @param temperature
             * @return
             */
            inline virtual double saturationPressure(double temperature) const = 0;

            /**
             * @brief
             * @param pressure
             * @return
             */
            inline virtual double saturationTemperature(double pressure) const = 0;

            /**
             * @brief
             * @param pureComponent
             */
            inline virtual void init(const IPureComponent& pureComponent) = 0;

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
            explicit Model(const T& x) : EOSType(x) {}

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
                return std::make_unique<Model<T>>(EOSType);
            }

            /**
             * @brief
             * @param pressure
             * @param temperature
             * @return
             */
            inline JSONString computePropertiesPT(double pressure, double temperature) const override
            {
                return EOSType.computePropertiesPT(pressure, temperature);
            }

            /**
             * @brief
             * @param temperature
             * @param molarVolume
             * @return
             */
            inline JSONString computePropertiesTV(double temperature, double molarVolume) const override
            {
                return EOSType.computePropertiesTV(temperature, molarVolume);
            }

            /**
             * @brief
             * @param temperature
             * @param molarVolume
             * @return
             */
            inline JSONString computePropertiesPV(double pressure, double molarVolume) const override
            {
                return EOSType.computePropertiesPV(pressure, molarVolume);
            }

            /**
             * @brief
             * @param temperature
             * @return
             */
            inline double saturationPressure(double temperature) const override
            {
                return EOSType.computePSat(temperature);
            }

            /**
             * @brief
             * @param pressure
             * @return
             */
            inline double saturationTemperature(double pressure) const override
            {
                return EOSType.computeTSat(pressure);
            }

            /**
             *
             * @param pureComponent
             */
            inline void init(const IPureComponent& pureComponent) override {
                EOSType.init(pureComponent);
            }


        private:
            T EOSType;
        };

        std::unique_ptr<Concept> m_equationOfState;
    };

}    // namespace PCProps

#endif    // PCPROPS_IEQUATIONOFSTATE_HPP
