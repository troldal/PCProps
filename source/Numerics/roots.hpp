//
// Created by Kenneth Balslev on 16/01/2021.
//

#ifndef PCPROPS_ROOTS_HPP
#define PCPROPS_ROOTS_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "differentiation.hpp"

namespace numeric::impl
{

    /**
     * @brief
     * @param func
     * @param x
     * @param eps
     * @param maxiter
     * @return
     */
    inline double newton(const std::function<double(double)>& func, double x, double eps = 1E-6, int maxiter = 20)
    {
        using std::abs;
        int counter = 0;
        while (true) {
            double x1 = x - (func(x) / numeric::diff_central(func, x));
            if (abs(x - x1) < eps) return x1;
            if (counter > maxiter) return static_cast<double>(NAN);
            x = x1;
            ++counter;
        }
    }

    /**
     * @brief
     * @param objective
     * @param x1
     * @param x2
     * @param eps
     * @param max_iter
     * @return
     */
    inline double ridders(const std::function<double(double)>& objective, double x1, double x2, double eps = 1.0E-6, double max_iter = 100)
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
            values[3] = values[2] + (values[2] - values[0]) * ((sign * objective(values[2])) / sqrt(pow(objective(values[2]), 2) - objective(values[0]) * objective(values[1])));

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

    /**
     * @brief
     * @param objective
     * @param lower
     * @param upper
     * @return
     */
    inline std::pair<double, double> bracket_search_up(const std::function<double(double)>& objective, double lower, double upper, double max_iter = 100)
    {
        if (upper <= lower) throw std::logic_error("Upper value must be higher than the lower value!");

        auto diff = upper - lower;
        for (int i = 0; i < max_iter; ++i) {
            if (objective(lower) * objective(upper) < 0.0) return std::make_pair(lower, upper);
            lower = upper;
            upper += diff;
        }

        throw std::logic_error("Bracket not found!");
    }

    inline std::vector<double> solve_cubic(double coeff_a, double coeff_b, double coeff_c)
    {
        using std::abs;
        using std::acos;
        using std::cos;
        using std::pow;
        using std::sqrt;

        std::vector<double> result;

        auto q = (pow(coeff_a, 2) - 3 * coeff_b);
        auto r = (2 * pow(coeff_a, 3) - 9 * coeff_a * coeff_b + 27 * coeff_c);

        auto Q = q / 9;
        auto R = r / 54;

        auto Q3 = pow(Q, 3);
        auto R2 = pow(R, 2);

        //        auto CR2 = 729 * r * r;
        //        auto CQ3 = 2916 * q * q * q;

        // ===== Three real roots
        //        if (R == 0 && Q == 0)
        //            result = { -coeff_a / 3, -coeff_a / 3, -coeff_a / 3 };
        //
        //        else if (CR2 == CQ3) {
        //            double sqrtQ = sqrt(Q);
        //
        //            if (R > 0)
        //                result = { -2 * sqrtQ - coeff_a / 3, sqrtQ - coeff_a / 3, sqrtQ - coeff_a / 3 };
        //
        //            else
        //                result = { -sqrtQ - coeff_a / 3, -sqrtQ - coeff_a / 3, 2 * sqrtQ - coeff_a / 3 };
        //        }
        //        else
        if (R2 < Q3) {
            auto sgnR  = (R >= 0 ? 1 : -1);
            auto ratio = sgnR * sqrt(R2 / Q3);
            auto theta = acos(ratio);
            auto norm  = -2 * sqrt(Q);

            result = { norm * cos(theta / 3) - coeff_a / 3, norm * cos((theta + 2.0 * M_PI) / 3) - coeff_a / 3, norm * cos((theta - 2.0 * M_PI) / 3) - coeff_a / 3 };
        }

        // ===== One real root
        else {
            auto sgnR = (R >= 0 ? 1 : -1);
            auto A    = -sgnR * pow(abs(R) + sqrt(R2 - Q3), 1.0 / 3.0);
            auto B    = Q / A;
            result    = { A + B - coeff_a / 3 };
        }

        std::sort(result.begin(), result.end());
        return result;



        //            // ===== Compute the constants required for an analytic solution.
        //            auto p = (1.0 / 3.0) * (3.0 * a_1 - pow(a_2, 2.0));
        //            auto q = (1.0 / 27.0) * (2.0 * pow(a_2, 3.0) - 9.0 * a_2 * a_1 + 27.0 * a_0);
        //            auto R = (pow(q, 2.0) / 4.0) + (pow(p, 3.0) / 27.0);
        //
        //            // ===== If R <= 0, there are three real roots
        //            if (R <= 0.0) {
        //                auto m     = 2.0 * sqrt(-p / 3.0);
        //                auto theta = acos(3 * q / (p * m)) / 3.0;
        //
        //                std::vector<double> roots { m * cos(theta) - a_2 / 3.0, m * cos(theta + 2.0 * PI / 3.0) - a_2 / 3.0, m * cos(theta + 4.0 * PI / 3.0) - a_2 / 3.0 };
        //
        //                std::sort(roots.begin(), roots.end());
        //                roots.erase(roots.begin() + 1);
        //
        //                // ===== If any of the roots are negative, delete them.
        //                if (roots[1] <= 0.0) roots.pop_back();
        //                if (roots[0] <= 0.0) roots.erase(roots.begin());
        //
        //                // ===== If the two roots are equal, delete the last one.
        //                if (roots.size() == 2 && roots[0] == roots[1]) roots.pop_back();
        //
        //                return roots;
        //            }
        //
        //            // ===== If R > 0, there is one real root
        //            auto P = cbrt(-q / 2.0 + sqrt(R));
        //            auto Q = cbrt(-q / 2.0 - sqrt(R));
        //
        //            return { P + Q - a_2 / 3.0 };
    }
}    // namespace numeric

#endif    // PCPROPS_ROOTS_HPP
