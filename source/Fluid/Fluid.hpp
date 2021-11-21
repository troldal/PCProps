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

#ifndef PCPROPS_FLUID_HPP
#define PCPROPS_FLUID_HPP

#include <IEquationOfState.hpp>
#include <IPureComponent.hpp>

#include <memory>

namespace PCProps
{

    /**
     * @brief The Fluid class encapsulates the concept of a fluid (liquid and/or gas). It consists of a object
     * representing a pure component (Note: in the future it may be multiple components) and optionally an object
     * representing an equation of state. Equation of state can be of any type (cubic, Excess-G, high-procision, custom)
     * as long as it supports the required interface.
     */
    class Fluid
    {
        using JSONString = std::string;
    public:
        // =====================================================================
        // CONSTRUCTORS & ASSIGNMENT OPERATORS
        // =====================================================================

        /**
         * @brief Constructor, default
         */
        Fluid();

        /**
         * @brief Constructor, taking a pure component object (any object supporting the IPureComponent interface), and
         * an optional equation of state object (any object supporting the IEquationOfState interface).
         * @param pureComponent Object representing a pure component.
         * @param eos Object representing an equation of state.
         */
        explicit Fluid(const IPureComponent& pureComponent, const IEquationOfState& eos = {});

        /**
         * @brief Copy constructor
         */
        Fluid(const Fluid& other);

        /**
         * @brief Move constructor
         */
        Fluid(Fluid&& other) noexcept;

        /**
         * @brief Destructor
         */
        ~Fluid();

        /**
         * @brief Copy assignment operator
         */
        Fluid& operator=(const Fluid& other);

        /**
         * @brief Move assignment operator
         */
        Fluid& operator=(Fluid&& other) noexcept;

        // =====================================================================
        // FLASH ALGORITHMS
        // =====================================================================

        /**
         * @brief Flash the fluid at the specified pressure and temperature.
         * @param pressure The pressure [Pa]
         * @param temperature The temperature [K]
         * @return A JSONString object with the phase properties.
         */
        JSONString flashPT(double pressure, double temperature) const;

        /**
         * @brief Flash the fluid at the specified pressure and vapor fraction.
         * @param pressure The pressure [Pa]
         * @param vaporFraction The vapor fraction [-]
         * @return A JSONString object with the phase properties.
         */
        JSONString flashPx(double pressure, double vaporFraction) const;

        /**
         * @brief Flash the fluid at the specified temperature and vapor fraction.
         * @param temperature The temperature [K]
         * @param vaporFraction The vapor fraction [-]
         * @return A JSONString object with the phase properties.
         */
        JSONString flashTx(double temperature, double vaporFraction) const;

        /**
         * @brief Flash the fluid at the specified pressure and enthalpy.
         * @param pressure The pressure [Pa]
         * @param enthalpy The enthalpy [J/mol]
         * @return A JSONString object with the phase properties.
         */
        JSONString flashPH(double pressure, double enthalpy) const;

        /**
         * @brief Flash the fluid at the specified pressure and entropy.
         * @param pressure The pressure [Pa]
         * @param entropy The entropy [J/mol-K]
         * @return A JSONString object with the phase properties.
         */
        JSONString flashPS(double pressure, double entropy) const;

        /**
         * @brief Flash the fluid at the specified temperature and molar volume.
         * @param temperature The temperature [K]
         * @param volume The molar volume [m3/mol]
         * @return A JSONString object with the phase properties.
         */
        JSONString flashTV(double temperature, double volume) const;

        /**
         * @brief Get the fluid properties in its current state.
         * @return A JSONString object with the phase properties.
         */
        JSONString properties() const;

    private:
        class impl;
        std::unique_ptr<impl> m_impl;
    };
}    // namespace PCProps

#endif    // PCPROPS_FLUID_HPP
