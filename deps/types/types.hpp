//
// Created by Kenneth Balslev on 16/01/2021.
//

#ifndef PCPROPS_TYPES_HPP
#define PCPROPS_TYPES_HPP

namespace types
{

    /**
     * @brief
     * @tparam T
     * @tparam Parameter
     */
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

} // namespace types

#endif    // PCPROPS_TYPES_HPP
