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

#ifndef PCPROPS_IUNITOPERATION_HPP
#define PCPROPS_IUNITOPERATION_HPP

#include <functional>
#include <memory>
#include <string>

namespace PCProps
{
    /**
     * @brief This class functions as a wrapper around any class that provides the necessary functionality for
     * equations of state.
     * @details This class works by applying 'type erasure'. This enables the use of objects of any class, the only
     * requirement being that it provides the right interface. No inheritance from a base class is needed.
     */
    class IUnitOperation
    {
        using JSONString = std::string;

    public:
        /**
         * @brief Default constructor
         */
        IUnitOperation() : m_unitOperation() {}    // NOLINT

        /**
         * @brief Constructor, taking the target object as an argument.
         * @tparam T The type of the target object (will be auto deducted)
         * @param x The target object
         * @note This methed is deliberately not marked 'explicit', because as a templated constructor, it should be able
         * to take any type as an argument. However, only objects that satisfy the required interface can be used.
         */
        template<typename T>
        IUnitOperation(const T& unitOperation) : m_unitOperation { std::make_unique<Model<T>>(unitOperation) }
        {}    // NOLINT

        /**
         * @brief Copy constructor
         * @param other
         */
        IUnitOperation(const IUnitOperation& other) : m_unitOperation(other.m_unitOperation ? other.m_unitOperation->clone() : nullptr) {}

        /**
         * @brief Move constructor
         * @param other
         */
        IUnitOperation(IUnitOperation&& other) noexcept = default;

        /**
         * @brief Destructor
         */
        ~IUnitOperation() = default;

        /**
         * @brief
         * @tparam T
         * @param x
         * @return
         */
        template<typename T>
        inline IUnitOperation& operator=(const T& unitOperation)
        {
            m_unitOperation = std::make_unique<Model<T>>(unitOperation);
            return *this;
        }

        /**
         * @brief
         * @param other
         * @return
         */
        inline IUnitOperation& operator=(const IUnitOperation& other)
        {
            IUnitOperation copy(other);
            *this = std::move(copy);
            return *this;
        }

        /**
         * @brief
         * @param other
         * @return
         */
        inline IUnitOperation& operator=(IUnitOperation&& other) noexcept = default;

        /**
         * @brief
         * @return
         */
        inline explicit operator bool() const { return m_unitOperation != nullptr; }

        inline void setSpecification(const JSONString& specification) {
            m_unitOperation->setSpecification(specification);
        }

        inline void compute() {
            m_unitOperation->compute();
        }

        inline JSONString results() const  {
            return m_unitOperation->results();
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


            inline virtual void setSpecification(const JSONString& specification) = 0;

            inline virtual void compute() = 0;

            inline virtual JSONString results() const = 0;

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
            explicit Model(const T& x) : ConcreteType(x) {}

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
            inline std::unique_ptr<Concept> clone() const override { return std::make_unique<Model<T>>(ConcreteType); }

            /**
             * @brief
             * @param specification
             */
            inline void setSpecification(const JSONString& specification) override {
                ConcreteType.setSpecification(specification);
            }

            /**
             * @brief
             */
            inline void compute() override {
                ConcreteType.compute();
            }

            /**
             * @brief
             * @return
             */
            inline JSONString results() const override {
                return ConcreteType.results();
            }

        private:
            T ConcreteType;
        };

        std::unique_ptr<Concept> m_unitOperation;
    };

}    // namespace PCProps

#endif    // PCPROPS_IUNITOPERATION_HPP
