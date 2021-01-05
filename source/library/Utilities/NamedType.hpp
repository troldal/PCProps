//
// Created by Kenneth Balslev on 04/01/2021.
//

#ifndef PCPROPS_NAMEDTYPE_HPP
#define PCPROPS_NAMEDTYPE_HPP

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
}    // namespace PCProps::Utilities

#endif    // PCPROPS_NAMEDTYPE_HPP
