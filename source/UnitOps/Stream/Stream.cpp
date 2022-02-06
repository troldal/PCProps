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

#include "Stream.hpp"

#include <FluidProperties.hpp>

// ===== External Headers ===== //
#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

using JSONString = std::string;

namespace PCProps::UnitOps
{

    /**
     * @brief
     */
    class Stream::impl
    {
        enum class FlashSpecification { PT, Px, Tx, PH, PS, TV };

        FlashSpecification m_flashSpec { FlashSpecification::PT };

        IPropertyPackage m_fluid {};
        FluidProperties  m_fluidProps {};

        double m_specOne {};
        double m_specTwo {};
        double m_quantity {};

    public:

        /**
         * @brief
         * @param fluid
         * @param quantity
         */
        impl(const IPropertyPackage& fluid, double quantity) : m_fluid { fluid }, m_quantity { quantity } {}

        /**
         * @brief
         */
        void setQuantity()
        {
            for (auto& phase : m_fluidProps) phase.MolarFlow = phase.MolarFraction * m_quantity;
        }

        /**
         * @brief
         * @param pressure
         * @param temperature
         * @return
         */
        void flashPT(double pressure, double temperature) {
            m_fluidProps = FluidProperties(m_fluid.flash("PT", pressure, temperature));
            setQuantity();
        }

        /**
         * @brief
         * @param pressure
         * @param vaporFraction
         * @return
         */
        void flashPx(double pressure, double vaporFraction) {
            m_fluidProps = FluidProperties(m_fluid.flash("Px", pressure, vaporFraction));
            setQuantity();
        }

        /**
         * @brief
         * @param temperature
         * @param vaporFraction
         * @return
         */
        void flashTx(double temperature, double vaporFraction){
            m_fluidProps = FluidProperties(m_fluid.flash("Tx", temperature, vaporFraction));
            setQuantity();
        }

        /**
         * @brief
         * @param pressure
         * @param enthalpy
         * @return
         */
        void flashPH(double pressure, double enthalpy) {
            m_fluidProps = FluidProperties(m_fluid.flash("PH", pressure, enthalpy));
            setQuantity();
        }

        /**
         * @brief
         * @param pressure
         * @param entropy
         * @return
         */
        void flashPS(double pressure, double entropy){
            m_fluidProps = FluidProperties(m_fluid.flash("PS", pressure, entropy));
            setQuantity();
        }

        /**
         * @brief
         * @param temperature
         * @param volume
         * @return
         */
        void flashTV(double temperature, double volume)
        {
            m_fluidProps = FluidProperties(m_fluid.flash("TV", temperature, volume));
            setQuantity();
        }

        /**
         * @brief
         * @param specification
         */
        void setSpecification(const JSONString& specification)
        {
            double pressure {};
            double temperature {};
            double vaporFraction {};
            double molarEnthalpy {};
            double molarEntropy {};
            double molarVolume {};

            using rapidjson::Document;
            Document pumpspec;
            pumpspec.Parse(specification.c_str());

            for (const auto& item : pumpspec.GetObject()) {
                std::string key = item.name.GetString();

                if (key == "FlashSpecification") {
                    std::string value = item.value.GetString();

                    if (value == "PT")
                        m_flashSpec = FlashSpecification::PT;
                    else if (value == "Px")
                        m_flashSpec = FlashSpecification::Px;
                    else if (value == "Tx")
                        m_flashSpec = FlashSpecification::Tx;
                    else if (value == "PH")
                        m_flashSpec = FlashSpecification::PH;
                    else if (value == "PS")
                        m_flashSpec = FlashSpecification::PS;
                    else if (value == "TV")
                        m_flashSpec = FlashSpecification::TV;
                }
                else if (key == "Pressure")
                    pressure = item.value.GetDouble();
                else if (key == "Temperature")
                    temperature = item.value.GetDouble();
                else if (key == "VaporFraction")
                    vaporFraction = item.value.GetDouble();
                else if (key == "MolarEnthalpy")
                    molarEnthalpy = item.value.GetDouble();
                else if (key == "MolarEntropy")
                    molarEntropy = item.value.GetDouble();
                else if (key == "MolarVolume")
                    molarVolume = item.value.GetDouble();

                switch (m_flashSpec) {
                    case FlashSpecification::PT:
                        m_specOne = pressure;
                        m_specTwo = temperature;
                        break;
                    case FlashSpecification::Px:
                        m_specOne = pressure;
                        m_specTwo = vaporFraction;
                        break;
                    case FlashSpecification::Tx:
                        m_specOne = temperature;
                        m_specTwo = vaporFraction;
                        break;
                    case FlashSpecification::PH:
                        m_specOne = pressure;
                        m_specTwo = molarEnthalpy;
                        break;
                    case FlashSpecification::PS:
                        m_specOne = pressure;
                        m_specTwo = molarEntropy;
                        break;
                    case FlashSpecification::TV:
                        m_specOne = temperature;
                        m_specTwo = molarVolume;
                        break;
                }
            }
        }

        /**
         * @brief
         */
        void compute()
        {
            switch (m_flashSpec) {
                case FlashSpecification::PT:
                    flashPT(m_specOne, m_specTwo);
                    break;
                case FlashSpecification::Px:
                    flashPx(m_specOne, m_specTwo);
                    break;
                case FlashSpecification::Tx:
                    flashTx(m_specOne, m_specTwo);
                    break;
                case FlashSpecification::PH:
                    flashPH(m_specOne, m_specTwo);
                    break;
                case FlashSpecification::PS:
                    flashPS(m_specOne, m_specTwo);
                    break;
                case FlashSpecification::TV:
                    flashTV(m_specOne, m_specTwo);
                    break;
            }
        }

        /**
         * @brief
         * @return
         */
        JSONString results() const { return m_fluidProps.asJSON(); }
    };

    /**
     * @details
     */
    Stream::Stream() = default;

    /**
     * @details
     */
    Stream::Stream(const IPropertyPackage& fluid, double quantity) : m_impl(std::make_unique<impl>(fluid, quantity)) {}

    /**
     * @details
     */
    Stream::Stream(const Stream& other) : m_impl(std::make_unique<impl>(*other.m_impl)) {};

    /**
     * @details
     */
    Stream::Stream(Stream&& other) noexcept = default;

    /**
     * @details
     */
    Stream::~Stream() = default;

    /**
     * @details
     */
    Stream& Stream::operator=(const Stream& other)
    {
        Stream copy = other;
        *this       = std::move(copy);
        return *this;
    }

    /**
     * @details
     */
    Stream& Stream::operator=(Stream&& other) noexcept = default;

    /**
     * @details
     */
    JSONString Stream::flash(const std::string& flashSpec, double specOne, double specTwo)
    {
        if (flashSpec == "PT")
            m_impl->flashPT(specOne, specTwo);
        else if (flashSpec == "Px")
            m_impl->flashPx(specOne, specTwo);
        else if (flashSpec == "Tx")
            m_impl->flashTx(specOne, specTwo);
        else if (flashSpec == "PH")
            m_impl->flashPH(specOne, specTwo);
        else if (flashSpec == "PS")
            m_impl->flashPS(specOne, specTwo);
        else if (flashSpec == "TV")
            m_impl->flashTV(specOne, specTwo);

        return m_impl->results();
    }

    /**
     * @details
     */
    void Stream::setSpecification(const JSONString& specification) { m_impl->setSpecification(specification); }

    /**
     * @details
     */
    void Stream::compute() { m_impl->compute(); }

    /**
     * @details
     */
    JSONString Stream::results() const { return m_impl->results(); }

}    // namespace PCProps::UnitOps