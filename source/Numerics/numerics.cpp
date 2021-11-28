//
// Created by Kenneth Balslev on 28/11/2021.
//

#include "numerics.hpp"
#include "roots.hpp"
#include "integration.hpp"
#include "differentiation.hpp"

#include <vector>
#include <functional>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

namespace {
    double function_wrapper_f(double arg, void * params) {
        return (*reinterpret_cast<std::function<double(double)>*>(params))(arg);
    }

    double function_wrapper_df(double arg, void * params) {
        auto f = (*reinterpret_cast<std::function<double(double)>*>(params));
        return numeric::diff_central(f, arg);
    }

    void function_wrapper_fdf(double arg, void * params, double * f, double * df) {
        *f = function_wrapper_f(arg, params);
        *df = function_wrapper_df(arg, params);
    }
}

namespace numeric {

    double newton(const std::function<double(double)>& func, double arg, double eps, int maxiter)
    {
//        auto function = func;
//        int status;
//        int iter = 0, max_iter = 100;
//        const gsl_root_fdfsolver_type *T;
//        gsl_root_fdfsolver *s;
//        double x0, x = arg;
//        gsl_function_fdf FDF;
//
//        FDF.f = &function_wrapper_f;
//        FDF.df = &function_wrapper_df;
//        FDF.fdf = &function_wrapper_fdf;
//        FDF.params = &function;
//
//        T = gsl_root_fdfsolver_newton;
//        s = gsl_root_fdfsolver_alloc (T);
//        gsl_root_fdfsolver_set (s, &FDF, x);
//
//        do
//        {
//            iter++;
//            status = gsl_root_fdfsolver_iterate (s);
//            x0 = x;
//            x = gsl_root_fdfsolver_root (s);
//            status = gsl_root_test_delta (x, x0, 0, 1e-3);
//        }
//        while (status == GSL_CONTINUE && iter < max_iter);
//
//        gsl_root_fdfsolver_free (s);
//        return x;



        return impl::newton(func, arg, eps, maxiter);
    }
    double ridders(const std::function<double(double)>& objective, double x1, double x2, double eps, int max_it)
    {
//        auto function = objective;
//        int status;
//        int cur_iter = 0;
//        int max_iter = max_it;
//
//        const gsl_root_fsolver_type *T;
//        gsl_root_fsolver *s;
//        gsl_function F;
//
//        double r;
//        double x_lo = x1;
//        double x_hi = x2;
//
//        F.function = &function_wrapper_f;
//        F.params = &function;
//
//        T = gsl_root_fsolver_brent;
//        s = gsl_root_fsolver_alloc (T);
//        gsl_root_fsolver_set (s, &F, x_lo, x_hi);
//
//        do
//        {
//            cur_iter++;
//            status = gsl_root_fsolver_iterate (s);
//            r = gsl_root_fsolver_root (s);
//            x_lo = gsl_root_fsolver_x_lower (s);
//            x_hi = gsl_root_fsolver_x_upper (s);
//            status = gsl_root_test_interval (x_lo, x_hi,0, eps);
//        }
//        while (status == GSL_CONTINUE && cur_iter < max_iter);
//
//        gsl_root_fsolver_free (s);
//
//        return r;

        return impl::ridders(objective, x1, x2, eps, max_it);
    }

    std::pair<double, double> bracket_search_up(const std::function<double(double)>& objective, double lower, double upper, double max_iter)
    {
        return impl::bracket_search_up(objective, lower, upper, max_iter);
    }

    std::vector<double> solve_cubic(double coeff_a, double coeff_b, double coeff_c)
    {
//        std::vector<double> result {-1.0, -1.0, -1.0};
//        gsl_poly_solve_cubic(coeff_a, coeff_b, coeff_c, &result[0], &result[1], &result[2]);
//        return result;

        return impl::solve_cubic(coeff_a, coeff_b, coeff_c);
    }



    double integrate(const std::function<double(double)>& func, double x1, double x2, double precision)
    {
//        auto temp = func;
//        gsl_integration_workspace * w
//            = gsl_integration_workspace_alloc (1000);
//
//        double result, error;
//
//        gsl_function F;
//        F.function = &function_wrapper_f;
//        F.params = &temp;
//
//        gsl_integration_qags (&F, x1, x2, 0, 1e-7, 1000,
//                             w, &result, &error);
//
//        gsl_integration_workspace_free (w);
//        return result;
        return impl::integrate(func, x1, x2, precision);
    }

    double diff_central(const std::function<double(double)>& func, double x)
    {
//        auto temp = func;
//        gsl_function F;
//        double result, abserr;
//
//        F.function = &function_wrapper_f;
//        F.params = &temp;
//        gsl_deriv_central (&F, x, 1e-8, &result, &abserr);
//
//        return result;

        return impl::diff_central(func, x);
    }

    double diff_backward(const std::function<double(double)>& func, double x)
    {
//        auto temp = func;
//        gsl_function F;
//        double result, abserr;
//
//        F.function = &function_wrapper_f;
//        F.params = &temp;
//        gsl_deriv_backward(&F, x, 1e-8, &result, &abserr);
//
//        return result;

        return impl::diff_backward(func, x);
    }

}