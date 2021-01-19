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
#include <types/types.hpp>

#include <memory>

namespace PCProps
{
    using Pressure      = types::NamedType<double, struct PressureTag>;
    using Temperature   = types::NamedType<double, struct TemperatureTag>;
    using VaporFraction = types::NamedType<double, struct VaporFractionTag>;
    using Enthalpy      = types::NamedType<double, struct EnthalpyTag>;
    using Entropy       = types::NamedType<double, struct EntropyTag>;
    using MolarVolume   = types::NamedType<double, struct MolarVolumeTag>;

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

        PCPhases flash(Pressure pressure, Temperature temperature) const;

        PCPhases flash(Pressure pressure, VaporFraction vaporFraction) const;

        PCPhases flash(Temperature temperature, VaporFraction vaporFraction) const;

        PCPhases flash(Pressure pressure, Enthalpy enthalpy) const;

        PCPhases flash(Pressure pressure, Entropy entropy) const;

        PCPhases flash(Temperature temperature, MolarVolume volume) const;

    private:
        class impl;
        std::unique_ptr<impl> m_impl;
    };
}    // namespace PCProps

#endif    // PCPROPS_FLUID_HPP
