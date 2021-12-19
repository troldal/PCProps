//
// Created by Kenneth Balslev on 16/01/2021.
//

#ifndef PCPROPS_DIFFERENTIATION_HPP
#define PCPROPS_DIFFERENTIATION_HPP

#include <functional>
#include <cmath>

namespace numeric::impl {

    /**
     * @brief
     * @param func
     * @param x
     * @return
     */
    inline double diff_central(const std::function<double(double)>& func, double x)
    {
        using std::sqrt;

        auto h = sqrt(std::numeric_limits<double>::epsilon());
        auto result = (func(x - 2*h) - 8*func(x - h) + 8*func(x + h) - func(x + 2*h)) / (12 * h);
        if (std::isnan(result)) result = (func(x + h) - func(x - h)) / (2 * h);
        if (std::isnan(result)) throw std::runtime_error("Derivative could not be computed.");

        return result;

    }

    /**
     * @brief
     * @param func
     * @param x
     * @return
     */
    inline double diff_forward(const std::function<double(double)>& func, double x)
    {
        using std::sqrt;

        auto h = sqrt(std::numeric_limits<double>::epsilon());
        auto result = (-25*func(x) + 48*func(x + h) - 36*func(x + 2*h) + 16*func(x + 3*h) - 3*func(x + 4*h)) / (12 * h);
        if (std::isnan(result)) result = (-3 * func(x) + 4 * func(x + h) - func(x + 2*h)) / (2 * h);
        if (std::isnan(result)) throw std::runtime_error("Derivative could not be computed.");

        return result;
    }

    /**
     * @brief
     * @param func
     * @param x
     * @return
     */
    inline double diff_backward(const std::function<double(double)>& func, double x)
    {
        using std::sqrt;

        auto h = -sqrt(std::numeric_limits<double>::epsilon());
        auto result = (-25*func(x) + 48*func(x + h) - 36*func(x + 2*h) + 16*func(x + 3*h) - 3*func(x + 4*h)) / (12 * h);
        if (std::isnan(result)) result = (-3 * func(x) + 4 * func(x + h) - func(x + 2*h)) / (2 * h);
        if (std::isnan(result)) throw std::runtime_error("Derivative could not be computed.");

        return result;
    }

} // namespace numeric

#endif    // PCPROPS_DIFFERENTIATION_HPP
