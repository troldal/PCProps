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

#ifndef PCPROPS_STREAM_HPP
#define PCPROPS_STREAM_HPP

#include <IPropertyPackage.hpp>

#include <memory>
#include <string>

namespace PCProps::UnitOps
{
    /**
     * @brief
     */
    class Stream
    {
        using JSONString = std::string;

    public:

        /**
         * @brief
         */
        Stream();

        /**
         * @brief
         * @param fluid
         * @param quantity
         */
        Stream(const IPropertyPackage& fluid, double quantity);

        /**
         * @brief
         * @param other
         */
        Stream(const Stream& other);

        /**
         * @brief
         * @param other
         */
        Stream(Stream&& other) noexcept;

        /**
         * @brief
         */
        ~Stream();

        /**
         * @brief
         * @param other
         * @return
         */
        Stream& operator=(const Stream& other);

        /**
         * @brief
         * @param other
         * @return
         */
        Stream& operator=(Stream&& other) noexcept;


        // ===== Generic Interface ===== //

        /**
         * @brief
         * @param specification
         */
        void setSpecification(const JSONString& specification);

        /**
         * @brief Calculate the stream properties, based on the specifications.
         */
        void compute();

        /**
         * @brief Get the calculation results for the stream.
         * @return A JSON string with the calculation results.
         */
        JSONString results() const;


        // ===== Specific Interface ===== //

        /**
         * @brief
         * @param spec
         * @param s1
         * @param s2
         * @return
         */
        JSONString flash(const std::string& flashSpec, double specOne, double specTwo);

    private:
        class impl;
        std::unique_ptr<impl> m_impl;

    };
}    // namespace PCProps::UnitOps

#endif    // PCPROPS_STREAM_HPP
