//
// Created by Kenneth Balslev on 04/01/2021.
//

#ifndef PCPROPS_NAMEDTYPES_HPP
#define PCPROPS_NAMEDTYPES_HPP

namespace PCProps::Utilities
{
    template<typename T, typename Parameter>
    class NamedType
    {
    public:
        explicit NamedType(T const& value) : value_(value) {}
        T& get()
        {
            return value_;
        }
        T const& get() const
        {
            return value_;
        }

    private:
        T value_;
    };

    using Temperature   = PCProps::Utilities::NamedType<double, struct TemperatureTag>;
    using Pressure      = PCProps::Utilities::NamedType<double, struct PressureTag>;
    using Volume        = PCProps::Utilities::NamedType<double, struct VolumeTag>;
    using Enthalpy      = PCProps::Utilities::NamedType<double, struct EnthalpyTag>;
    using Entropy       = PCProps::Utilities::NamedType<double, struct EntropyTag>;
    using VaporFraction = PCProps::Utilities::NamedType<double, struct VaporFractionTag>;

}    // namespace PCProps::Utilities

#endif    // PCPROPS_NAMEDTYPES_HPP
