//
// Created by Kenneth Balslev on 18/01/2021.
//

#ifndef PCPROPS_PCFLUID_HPP
#define PCPROPS_PCFLUID_HPP

#include <interfaces/PCEquationOfState.hpp>
#include <interfaces/PCPureComponent.hpp>
#include <types/types.hpp>

namespace PCProps
{
    using Pressure = types::NamedType<double, struct PressureTag>;
    using Temperature = types::NamedType<double, struct TemperatureTag>;
    using VaporFraction = types::NamedType<double, struct VaporFractionTag>;
    using Enthalpy = types::NamedType<double, struct EnthalpyTag>;
    using Entropy = types::NamedType<double, struct EntropyTag>;
    using MolarVolume = types::NamedType<double, struct MolarVolumeTag>;


    class PCFluid
    {

        PCPureComponent m_pureComponent {};
        PCEquationOfState m_equationOfState {};

    public:

        PCFluid();

        PCFluid(const PCPureComponent& pc, const PCEquationOfState& eos);

        PCFluid(const PCFluid& other);

        PCFluid(PCFluid&& other) noexcept;

        ~PCFluid();

        PCFluid& operator=(const PCFluid& other);

        PCFluid& operator=(PCFluid&& other) noexcept;

        PCPhases flash(Pressure pressure, Temperature temperature) const;

        PCPhases flash(Pressure pressure, VaporFraction vaporFraction) const;

        PCPhases flash(Temperature temperature, VaporFraction vaporFraction) const;

        PCPhases flash(Pressure pressure, Enthalpy enthalpy) const;

        PCPhases flash(Pressure pressure, Entropy entropy) const;

        PCPhases flash(Temperature temperature, MolarVolume volume) const;


    };
}    // namespace PCProps

#endif    // PCPROPS_PCFLUID_HPP
