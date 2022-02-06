/*

888    d8P  8888888b.
888   d8P   888   Y88b
888  d8P    888    888
888d88K     888   d88P 888d888 .d88b.  88888b.  .d8888b
8888888b    8888888P"  888P"  d88""88b 888 "88b 88K
888  Y88b   888        888    888  888 888  888 "Y8888b.
888   Y88b  888        888    Y88..88P 888 d88P      X88
888    Y88b 888        888     "Y88P"  88888P"   88888P'
                                       888
                                       888
                                       888

Copyright (c) 2022 Kenneth Troldal Balslev

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

#ifndef PCPROPS_VALVE_HPP
#define PCPROPS_VALVE_HPP

// ===== KProps Headers ===== //
#include <Stream/Stream.hpp>

// ===== Standard Library Headers ===== //
#include <string>
#include <vector>

namespace PCProps::UnitOps
{

    /**
     * @brief The Valve class is a part of the UnitOps sub-namespace, and encapsulates
     * the concept of a valve unit operation.
     */
    class Valve
    {
        using JSONString = std::string;

    public:

        /**
         * @brief Default constructor. The object is constructed with dP = 0.
         */
        Valve();

        /**
         * @brief Constructor taking a JSON string with the specification.
         * @param specification The JSON string input. Keys can either be "DeltaPressure" or "OutletPressure".
         */
        Valve(const JSONString& specification);

        /**
         * @brief Copy constructor.
         * @param other Object to be copied.
         */
        Valve(const Valve& other);

        /**
         * @brief Move constructor.
         * @param other Object to be moved.
         */
        Valve(Valve&& other) noexcept;

        /**
         * @brief Destructor.
         */
        ~Valve();

        /**
         * @brief Copy assignment operator.
         * @param other Object to be copied.
         */
        Valve& operator=(const Valve& other);

        /**
         * @brief Move assignment operator.
         * @param other Object to be moved.
         */
        Valve& operator=(Valve&& other) noexcept;

        /**
         * @brief Set the inlet streams to the valve.
         * @param stream A std::vector holding const Stream* pointers.
         * @note Only one inlet stream should be provided. If the input vector holds more than
         * one stream object, only the first one will be used; the rest will be ignored.
         */
        void setInletStreams(const std::vector<const Stream*> streams);

        /**
         * @brief Set the specification for the valve. This overwrites the specification given in the constructor.
         * @param specification The JSON string input. Keys can either be "DeltaPressure" or "OutletPressure".
         */
        void setSpecification(const JSONString& specification);

        /**
         * @brief Retreive the outlet streams from the valve.
         * @return A const reference to the outlet streams.
         * @note Only one outlet stream is provided. Hence, the size of the vector will always be exactly one.
         */
        const std::vector<Stream>& outletStreams();

        /**
         * @brief Calculate the outlet stream properties, based on the inlet streams and equipment specifications.
         */
        void compute();

        /**
         * @brief Get the calculation results (for the equipment item, not the outlet streams).
         * @return A JSON string with the calculation results.
         */
        JSONString results() const;

    private:

        class impl;
        std::unique_ptr<impl> m_impl;

    };
}    // namespace PCProps::UnitOps

#endif    // PCPROPS_VALVE_HPP
