//
// Created by Kenneth Balslev on 18/01/2021.
//

#ifndef PCPROPS_FLUID_HPP
#define PCPROPS_FLUID_HPP

#include <interfaces/IEquationOfState.hpp>
#include <interfaces/IPureComponent.hpp>
#include <types/types.hpp>

namespace PCProps
{
    using Pressure = types::NamedType<double, struct PressureTag>;
    using Temperature = types::NamedType<double, struct TemperatureTag>;
    using VaporFraction = types::NamedType<double, struct VaporFractionTag>;
    using Enthalpy = types::NamedType<double, struct EnthalpyTag>;
    using Entropy = types::NamedType<double, struct EntropyTag>;
    using MolarVolume = types::NamedType<double, struct MolarVolumeTag>;


    class Fluid
    {
        IPureComponent   m_pureComponent {};
        IEquationOfState m_equationOfState {};

    public:
        Fluid();

        Fluid(const IPureComponent& pc, const IEquationOfState& eos);

        Fluid(const Fluid& other);

        Fluid(Fluid&& other) noexcept;

        ~Fluid();

        Fluid& operator=(const Fluid& other);

        Fluid& operator=(Fluid&& other) noexcept;

        PCPhases flash(Pressure pressure, Temperature temperature) const;

        PCPhases flash(Pressure pressure, VaporFraction vaporFraction) const;

        PCPhases flash(Temperature temperature, VaporFraction vaporFraction) const;

        PCPhases flash(Pressure pressure, Enthalpy enthalpy) const;

        PCPhases flash(Pressure pressure, Entropy entropy) const;

        PCPhases flash(Temperature temperature, MolarVolume volume) const;


    };
}    // namespace PCProps

#endif    // PCPROPS_FLUID_HPP
