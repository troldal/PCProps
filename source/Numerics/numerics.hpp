//
// Created by Kenneth Balslev on 28/11/2021.
//

#ifndef PCPROPS_NUMERICS_HPP
#define PCPROPS_NUMERICS_HPP

#include <functional>

namespace numeric {

    double newton(const std::function<double(double)>& func, double x, double eps = 1E-6, int maxiter = 20);

    double ridders(const std::function<double(double)>& objective, double x1, double x2, double eps = 1.0E-6, int max_iter = 100);

    std::pair<double, double> bracket_search_up(const std::function<double(double)>& objective, double lower, double upper, double max_iter = 100);

    std::vector<double> solve_cubic(double coeff_a, double coeff_b, double coeff_c);

    double integrate(const std::function<double(double)>& func, double x1, double x2, double precision = 1E-6);

    double diff_central(const std::function<double(double)>& func, double x);

    double diff_backward(const std::function<double(double)>& func, double x);

}

#endif    // PCPROPS_NUMERICS_HPP
