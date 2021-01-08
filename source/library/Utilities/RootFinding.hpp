//
// Created by Kenneth Balslev on 29/12/2020.
//

#ifndef PCPROPS_ROOTFINDING_HPP
#define PCPROPS_ROOTFINDING_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>

#include "Calculus.hpp"

namespace PCProps::Numerics
{
    double newton(const std::function<double(double)>& func, double x, double eps = 1E-6, int maxiter = 10)
    {
        using std::abs;
        int counter = 0;
        while (true) {
            double x1 = x - (func(x) / diff_central(func, x));
            if (abs(x - x1) < eps) return x1;
            if (counter > maxiter) return NAN;
            x = x1;
            ++counter;
        }
    }

    double ridders(const std::function<double(double)>& objective, double x1, double x2, double eps = 1.0E-6, double max_iter = 100)
    {
        using std::abs;
        using std::pow;
        using std::sqrt;

        int                   sign = 0;
        std::array<double, 4> values {};
        values[0] = x1;
        values[1] = x2;

        int counter = 0;

        while (true) {
            if (counter > max_iter) break;
            values[2] = (values[0] + values[1]) / 2.0;
            sign      = ((objective(values[0]) - objective(values[1])) < 0.0 ? -1 : 1);
            values[3] =
                values[2] + (values[2] - values[0]) * ((sign * objective(values[2])) /
                                                       sqrt(pow(objective(values[2]), 2) - objective(values[0]) * objective(values[1])));

            if (abs(objective(values[3])) < eps) break;

            std::sort(values.begin(), values.end());
            for (std::size_t i = 0; i <= 2; ++i)
                if (objective(values[i]) * objective(values[i + 1]) < 0.0) {
                    values[0] = values[i];
                    values[1] = values[i + 1];
                }
            ++counter;
        }

        return values[3];
    }

}    // namespace PCProps::Numerics

#endif    // PCPROPS_ROOTFINDING_HPP
