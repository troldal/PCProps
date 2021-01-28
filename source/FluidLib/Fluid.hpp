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

#include <interfaces/IEquationOfState.hpp>
#include <interfaces/IPureComponent.hpp>

#include <memory>

namespace PCProps
{
    class Fluid
    {
    public:
        // =====================================================================
        // CONSTRUCTORS & ASSIGNMENT OPERATORS
        // =====================================================================


        Fluid();

        Fluid(const IPureComponent& pc, const IEquationOfState& eos);

        Fluid(const Fluid& other);

        Fluid(Fluid&& other) noexcept;

        ~Fluid();

        Fluid& operator=(const Fluid& other);

        Fluid& operator=(Fluid&& other) noexcept;

        // =====================================================================
        // FLASH ALGORITHMS
        // =====================================================================

        /**
         * @brief
         * @param pressure
         * @param temperature
         * @return
         */
        const PCPhases& flashPT(double pressure, double temperature) const;

        /**
         * @brief
         * @param pressure
         * @param vaporFraction
         * @return
         */
        const PCPhases& flashPx(double pressure, double vaporFraction) const;

        /**
         * @brief
         * @param temperature
         * @param vaporFraction
         * @return
         */
        const PCPhases& flashTx(double temperature, double vaporFraction) const;

        /**
         * @brief
         * @param pressure
         * @param enthalpy
         * @return
         */
        const PCPhases& flashPH(double pressure, double enthalpy) const;

        /**
         * @brief
         * @param pressure
         * @param entropy
         * @return
         */
        const PCPhases& flashPS(double pressure, double entropy) const;

        /**
         * @brief
         * @param temperature
         * @param volume
         * @return
         */
        const PCPhases& flashTV(double temperature, double volume) const;

        /**
         * @brief
         * @return
         */
        const PCPhases& properties() const;

    private:
        class impl;
        std::unique_ptr<impl> m_impl;
    };
}    // namespace PCProps

#endif    // PCPROPS_FLUID_HPP
