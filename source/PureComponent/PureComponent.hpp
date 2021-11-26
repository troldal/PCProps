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

#ifndef PCPROPS_PURECOMPONENT_HPP
#define PCPROPS_PURECOMPONENT_HPP

#include <functional>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>
#include <variant>

namespace PCProps
{
    using ComponentDataItem = std::variant<std::string, double, std::function<double(double)>, std::function<double(std::vector<double>)> >;
    using ComponentDataRecord = std::tuple<std::string, ComponentDataItem, std::string, std::string>;

    /**
     * @brief A PCComponent object holds an PCComponentData object, as well as all required access functions.
     * @details In order for a PCComponent object to work properly, it should be constructed from a PCComponentData object.
     * The reason for dividing the definition of a pure component into two different classes is that the PCComponent does not
     * have any setter functions, an therefore cannot be modified after construction. The purpose is to avoid accidental
     * modification of pure component data after construction. If it is absolutely necessary to modify an existing PCComponent
     * object after construction, this should be done by copying the data object, modify the data and construct a new object
     * in the same place as the old object.
     */
    class PureComponent
    {
        // ===== Private Data Members
        std::vector<ComponentDataRecord> m_dataItems {};

        ComponentDataItem find(const std::string& ID) const {
            return std::get<1>(*std::find_if(m_dataItems.begin(),
                                             m_dataItems.end(),
                                             [&](const auto& item) {return std::get<0>(item) == ID;}));
        }


    public:

        // ===== Constructors & Assignment Operators ===== //

        /**
         * @brief Default constructor.
         */
        PureComponent() = default;

        /**
         * @brief Copy constructor.
         */
        PureComponent(const PureComponent& other) = default;

        /**
         * @brief Move constructor.
         */
        PureComponent(PureComponent&& other) noexcept = default;

        /**
         * @brief Destructor.
         */
        ~PureComponent() = default;

        /**
         * @brief Copy assignment operator.
         */
        PureComponent& operator=(const PureComponent& other) = default;

        /**
         * @brief Move assignment operator.
         */
        PureComponent& operator=(PureComponent&& other) noexcept = default;

        // ===== Accessors (Constants) ===== //

        void addDataItem(const std::string& ID, const ComponentDataItem& dataItem, const std::string& description = "", const std::string& unit = "") {
            m_dataItems.emplace_back(std::make_tuple(ID, dataItem, description, unit));
        }

        std::string metaData(const std::string& ID) {
            auto item = find(ID);
            if (std::holds_alternative<std::string>(item)) return std::get<std::string>(item);
            // TODO: else, throw exception
        }

        double property(const std::string& ID) const {
            auto item = find(ID);
            if (std::holds_alternative<double>(item)) return std::get<double>(item);
            // TODO: else, throw exception
        }

        double correlation(const std::string& ID, double temperature) const {
            auto item = find(ID);
            if (std::holds_alternative<std::function<double(double)>>(item))
                return std::get<std::function<double(double)>>(item)(temperature);
            // TODO: else, throw exception
        }

        double correlation(const std::string& ID, const std::vector<double>& parameters) const {
            auto item = find(ID);
            if (std::holds_alternative<std::function<double(std::vector<double>)>>(item))
                return std::get<std::function<double(std::vector<double>)>>(item)(parameters);
            // TODO: else, throw exception
        }

    };
}    // namespace PCProps

#endif    // PCPROPS_PURECOMPONENT_HPP
