/*

8888888b.   .d8888b.  8888888b.
888   Y88b d88P  Y88b 888   Y88b
888    888 888    888 888    888
888   d88P 888        888   d88P 888d888 .d88b.  88888b.  .d8888b
8888888P"  888        8888888P"  888P"  d88""88b 888 "88b 88K
888        888    888 888        888    888  888 888  888 "Y8888b.
888        Y88b  d88P 888        888    Y88..88P 888 d88P      X88
888         "Y8888P"  888        888     "Y88P"  88888P"   88888P'
                                                 888
                                                 888
                                                 888

Copyright (c) 2020 Kenneth Troldal Balslev

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "EOSPengRobinson.hpp"
#include <PCConfig.hpp>
#include <PCPropsException.hpp>
#include <Utilities/Integration.hpp>
#include <Utilities/RootFinding.hpp>

namespace
{
    // =========================================================================
    // ===== Enthalpy functions
    // =========================================================================

    /**
     * @brief Calculates the ideal gas enthalpy.
     * @details The calculation is done by integrating Cp (heat capacity at constant pressure) from the standard state
     * conditions (usually 298.15 K, but this can be changed by modifying the STANDARD_T global constant) to the given
     * temperature. This overload evaluates the integral directly, by using the supplied function object.
     * @param t The temperature [K] at which to evaluate the ideal gas enthalpy.
     * @param igCpIntegral A function object for evaluating the integral of the ideal gas heat capacity.
     * @return The ideal gas enthalpy [J/mol].
     */
    double idealGasEnthalpy(double t, const std::function<double(double, double)>& igCpIntegral)
    {
        return igCpIntegral(PCProps::Globals::STANDARD_T, t);
    }

    /**
     * @brief Calculates the ideal gas enthalpy.
     * @details The calculation is done by integrating Cp (heat capacity at constant pressure) from the standard state
     * conditions (usually 298.15 K, but this can be changed by modifying the STANDARD_T global constant) to the given
     * temperature. This overload evaluates the integral numerically, and should only be used if the integral
     * cannot be evaluated directly.
     * @param t The temperature [K] at which to evaluate the ideal gas enthalpy.
     * @param igCp A function object for evaluating the ideal gas heat capacity [J/mol-K].
     * @return The ideal gas enthalpy [J/mol].
     * @note This overload evaluates the Cp integral numerically. It should only be used if the integral cannot
     * be evaluated directly.
     */
    double idealGasEnthalpy(double t, const std::function<double(double)>& igCp)
    {
        using PCProps::Numerics::integrate;
        return integrate(igCp, PCProps::Globals::STANDARD_T, t);
    }

    /**
     * @brief Computes the enthalpy departure (H-H^ig) at constant T and P.
     * @param z The compressibility factor of the fluid.
     * @param t The temperature [K].
     * @param tc The critical temperature [K].
     * @param A The A coefficient evaluated at the given temperature.
     * @param B The B coefficient evaluated at the given temperature.
     * @param alpha The alpha coefficient evaluated at the given temperature.
     * @param kappa The kappa coefficient.
     * @return The enthalpy departure (H-H^ig) [J/mol]
     */
    double enthalpyDeparture(double z, double t, double tc, double A, double B, double alpha, double kappa)
    {
        using std::log;
        using std::sqrt;

        return (z - 1.0 -
                A / (B * sqrt(8)) * (1 + kappa * sqrt(t / tc) / sqrt(alpha)) * log((z + (1 + sqrt(2)) * B) / (z + (1 - sqrt(2)) * B))) *
               PCProps::Globals::R_CONST * t;
    }

    /**
     * @brief Computes the enthalpy, relative to the standard state (298.15 K by default, but can be modified)
     * @param t The temperature at which to evaluate the enthalpy [K].
     * @param tc The critical temperature [K].
     * @param z The compressibility factor [-].
     * @param A The A coefficient evaluated at the given temperature.
     * @param B The B coefficient evaluated at the given temperature.
     * @param alpha The alpha coefficient evaluated at the given temperature.
     * @param kappa The kappa coefficient.
     * @param igCp A function object for evaluating the ideal gas heat capacity.
     * @return The enthalpy relative to the standard state [J/mol]
     */
    double
        enthalpy(double t, double tc, double z, double A, double B, double alpha, double kappa, const std::function<double(double)>& igCp)
    {
        return idealGasEnthalpy(t, igCp) + enthalpyDeparture(z, t, tc, A, B, alpha, kappa);
    }

    // =========================================================================
    // ===== Entropy functions
    // =========================================================================

    /**
     * @brief Computes the ideal gas entropy.
     * @details The calculation is done by integrating Cp/T (heat capacity at constant pressure divided by T) from the s
     * tandard state conditions (usually 298.15 K, but this can be changed by modifying the STANDARD_T global constant)
     * to the given temperature. This overload evaluates the integral directly, by using the supplied function object.
     * @param t The temperature [K]
     * @param p The pressure [Pa]
     * @param igCpOverTIntegral A function object for evaluating the integral of Cp/T directly.
     * @return The ideal gas entropy [J/mol-K].
     */
    double idealGasEntropy(double t, double p, const std::function<double(double, double)>& igCpOverTIntegral)
    {
        return igCpOverTIntegral(PCProps::Globals::STANDARD_T, t) - PCProps::Globals::R_CONST * log(p / PCProps::Globals::STANDARD_P);
    }

    /**
     * @brief Computes the ideal gas entropy.
     * @details The calculation is done by integrating Cp/T (heat capacity at constant pressure divided by T) from the s
     * tandard state conditions (usually 298.15 K, but this can be changed by modifying the STANDARD_T global constant)
     * to the given temperature. This overload evaluates the integral numerically, and should only be used if the integral
     * cannot be evaluated directly.
     * @param t The temperature [K]
     * @param p The pressure [Pa]
     * @param igCp A function object for evaluating the ideal gas heat capacity [J/mol-K].
     * @return The ideal gas entropy [J/mol-K].
     * @note This overload evaluates the Cp/T integral numerically. It should only be used if the integral cannot
     * be evaluated directly.
     */
    double idealGasEntropy(double t, double p, const std::function<double(double)>& igCp)
    {
        using PCProps::Numerics::integrate;
        return integrate([&](double temp) { return igCp(temp) / temp; }, PCProps::Globals::STANDARD_T, t) -
               PCProps::Globals::R_CONST * log(p / PCProps::Globals::STANDARD_P);
    }

    /**
     * @brief omputes the entropy departure (S-S^ig) at constant T and P.
     * @param z The compressibility factor [-]
     * @param t The temperature [K]
     * @param tc The critical temperature [K]
     * @param A The A coefficient evaluated at the given temperature.
     * @param B The B coefficient evaluated at the given temperature.
     * @param alpha The alpha coefficient evaluated at the given temperature.
     * @param kappa The kappa coefficient.
     * @return The entropy departure (S-S^ig) [J/mol-K]
     */
    double entropyDeparture(double z, double t, double tc, double A, double B, double alpha, double kappa)
    {
        using std::log;
        using std::sqrt;

        return (log(z - B) -
                A / (B * sqrt(8)) * kappa * sqrt(t / tc) / sqrt(alpha) * log((z + (1 + sqrt(2)) * B) / (z + (1 - sqrt(2)) * B))) *
               PCProps::Globals::R_CONST;
    }

    /**
     * @brief
     * @param t
     * @param tc
     * @param p
     * @param z
     * @param A
     * @param B
     * @param alpha
     * @param kappa
     * @param igCp
     * @return
     */
    double entropy(
        double                               t,
        double                               tc,
        double                               p,
        double                               z,
        double                               A,
        double                               B,
        double                               alpha,
        double                               kappa,
        const std::function<double(double)>& igCp)
    {
        return idealGasEntropy(t, p, igCp) + entropyDeparture(z, t, tc, A, B, alpha, kappa);
    }

    // =========================================================================
    // ===== Utility functions
    // =========================================================================

    /**
     * @brief Calculates the compressibility factor(s)
     * @param A
     * @param B
     * @return
     */
    std::vector<double> computeCompressibilityFactors(double A, double B)
    {
        using std::acos;
        using std::cbrt;
        using std::cos;
        using std::sqrt;

        // ===== Compute the coefficients for solving Peng Robinson with respect to Z.
        auto a_0 = -(A * B - pow(B, 2) - pow(B, 3));
        auto a_1 = (A - 3 * pow(B, 2) - 2 * B);
        auto a_2 = -(1 - B);

        // ===== Compute the constants required for an analytic solution.
        auto p = (1.0 / 3.0) * (3 * a_1 - pow(a_2, 2));
        auto q = (1.0 / 27.0) * (2 * pow(a_2, 3) - 9 * a_2 * a_1 + 27 * a_0);
        auto R = (pow(q, 2) / 4.0) + (pow(p, 3) / 27.0);

        // ===== If R <= 0, there are three real roots
        if (R <= 0.0) {
            auto m     = 2 * sqrt(-p / 3);
            auto theta = acos(3 * q / (p * m)) / 3.0;

            std::vector<double> all_roots { m * cos(theta) - a_2 / 3,
                                            m * cos(theta + 2 * PCProps::Globals::PI / 3) - a_2 / 3,
                                            m * cos(theta + 4 * PCProps::Globals::PI / 3) - a_2 / 3 };

            std::vector<double> roots { *std::max_element(all_roots.begin(), all_roots.end()),
                                        *std::min_element(all_roots.begin(), all_roots.end()) };

            if (roots[0] == roots[1]) roots.erase(roots.begin());

            if (roots[1] <= 0.0) roots.pop_back();

            if (roots[0] <= 0.0) roots.pop_back();

            return roots;
        }

        // ===== If R > 0, there is one real root
        auto P = cbrt(-q / 2.0 + sqrt(R));
        auto Q = cbrt(-q / 2.0 - sqrt(R));

        return { P + Q - a_2 / 3.0 };
    }

    /**
     * @brief
     * @param A
     * @param B
     * @param Z
     * @return
     */
    inline double computeFugacityCoefficient(double A, double B, double Z)
    {
        using std::exp;
        using std::log;
        using std::sqrt;

        return exp(Z - 1 - log(Z - B) - A / (B * sqrt(8)) * log((Z + (1 + sqrt(2)) * B) / (Z + (1 - sqrt(2)) * B)));
    }

    /**
     * @brief
     * @param A
     * @param B
     * @return
     */
    std::vector<std::pair<double, double>> computeCompressibilityAndFugacity(double A, double B)
    {
        auto                                   zs = computeCompressibilityFactors(A, B);
        std::vector<std::pair<double, double>> result;
        for (const auto& z : zs) result.emplace_back(std::make_pair(z, computeFugacityCoefficient(A, B, z)));

        return result;
    }

    /**
     * @brief
     * @param temperature
     * @param fA
     * @param fB
     * @param fPSat
     * @return
     */
    double computeSaturationPressure(
        double                                       temperature,
        const std::function<double(double, double)>& fA,
        const std::function<double(double, double)>& fB,
        const std::function<double(double)>&         fPSat)
    {
        using std::get;
        auto f = [&](double p) {
            auto A     = fA(temperature, p);
            auto B     = fB(temperature, p);
            auto phi   = computeCompressibilityAndFugacity(A, B);
            auto phi_v = get<1>(phi[0]);
            auto phi_l = get<1>(phi[1]);
            return (phi_l - phi_v) * p;
        };

        return PCProps::Numerics::ridders(f, fPSat(temperature) * 0.8, fPSat(temperature) * 1.2);
    }

    /**
     * @brief
     * @param pressure
     * @param criticalTemperature
     * @param fA
     * @param fB
     * @param fPSat
     * @return
     */
    double computeSaturationTemperature(
        double                                       pressure,
        double                                       criticalTemperature,
        const std::function<double(double, double)>& fA,
        const std::function<double(double, double)>& fB,
        const std::function<double(double)>&         fPSat)
    {
        using std::get;
        auto f = [&](double t) {
            auto A     = fA(t, pressure);
            auto B     = fB(t, pressure);
            auto phi   = computeCompressibilityAndFugacity(A, B);
            auto phi_v = get<1>(phi[0]);
            auto phi_l = get<1>(phi[1]);
            return (phi_l - phi_v) * pressure;
        };

        auto g = [&](double t) { return fPSat(t) - pressure; };

        auto guess = PCProps::Numerics::ridders(g, 0, criticalTemperature);
        return PCProps::Numerics::ridders(f, guess * 0.8, guess * 1.2);
    }

}    // namespace

namespace PCProps::EquationOfState
{
    // ===== Constructor, default
    EOSPengRobinson::EOSPengRobinson() = default;

    // ===== Constructor
    EOSPengRobinson::EOSPengRobinson(
        double                               criticalTemperature,
        double                               criticalPressure,
        double                               acentricFactor,
        double                               molecularWeight,
        const std::function<double(double)>& vaporPressureFunction,
        const std::function<double(double)>& idealGasCpFunction)
        : m_criticalTemperature(criticalTemperature),
          m_criticalPressure(criticalPressure),
          m_molecularWeight(molecularWeight),
          m_vaporPressureFunction(vaporPressureFunction),
          m_idealGasCpFunction(idealGasCpFunction)
    {
        using std::pow;
        using std::sqrt;

        m_kappa   = 0.37464 + 1.54226 * acentricFactor - 0.26992 * pow(acentricFactor, 2);
        m_alpha   = [&](double t) { return pow(1 + m_kappa * (1 - sqrt(t / m_criticalTemperature)), 2); };
        double ac = 0.45723553 * pow(8.31446261815324, 2) * pow(criticalTemperature, 2) / criticalPressure;

        m_a = [=](double temperature) -> double { return ac * pow(1 + m_kappa * (1 - sqrt(temperature / criticalTemperature)), 2); };

        m_A = [=](double temperature, double pressure) -> double {
            return m_a(temperature) * pressure / pow(8.31446261815324, 2) / pow(temperature, 2);
        };

        m_b = 0.07779607 * 8.31446261815324 * criticalTemperature / criticalPressure;

        m_B = [=](double temperature, double pressure) -> double { return m_b * pressure / (8.31446261815324 * temperature); };
    }

    // ===== Copy constructor
    EOSPengRobinson::EOSPengRobinson(const EOSPengRobinson& other) = default;

    // ===== Move constructor
    EOSPengRobinson::EOSPengRobinson(EOSPengRobinson&& other) noexcept = default;

    // ===== Destructor
    EOSPengRobinson::~EOSPengRobinson() = default;

    // ===== Copy assignment operator
    EOSPengRobinson& EOSPengRobinson::operator=(const EOSPengRobinson& other) = default;

    // ===== Move assignment operator
    EOSPengRobinson& EOSPengRobinson::operator=(EOSPengRobinson&& other) noexcept = default;

    // ===== Helper function for creating output data
    PhaseData EOSPengRobinson::createEOSData(double moles, double temperature, double pressure, double z, double phi) const
    {
        using std::get;
        PhaseData result;

        auto A     = m_A(temperature, pressure);
        auto B     = m_B(temperature, pressure);
        auto alpha = m_alpha(temperature);

        get<Moles>(result)           = moles;
        get<MolecularWeight>(result) = m_molecularWeight;
        get<Temperature>(result)     = temperature;
        get<Pressure>(result)        = pressure;
        get<Volume>(result)          = z * PCProps::Globals::R_CONST * temperature / pressure;
        get<Fugacity>(result)        = phi * pressure;
        get<Compressibility>(result) = z;
        get<Enthalpy>(result)        = enthalpy(temperature, m_criticalTemperature, z, A, B, alpha, m_kappa, m_idealGasCpFunction);
        get<Entropy>(result)         = entropy(temperature, m_criticalTemperature, pressure, z, A, B, alpha, m_kappa, m_idealGasCpFunction);
        get<InternalEnergy>(result)  = get<Enthalpy>(result) - pressure * get<Volume>(result);
        get<GibbsEnergy>(result)     = get<Enthalpy>(result) - temperature * get<Entropy>(result);
        get<HelmholzEnergy>(result)  = get<InternalEnergy>(result) - temperature * get<Entropy>(result);

        return result;
    }

    // ===== P,T Flash
    Phases EOSPengRobinson::flashPT(double pressure, double temperature, double moles) const
    {
        using std::get;

        // ===== Compute compressibility factors and fugacity coefficients at given T and P.
        auto z_phi = computeCompressibilityAndFugacity(m_A(temperature, pressure), m_B(temperature, pressure));

        // ===== Identify the Z/Phi with the lowest fugacity coefficient (which is the most stable)
        auto min = std::min_element(z_phi.begin(), z_phi.end(), [](const auto& a, const auto& b) { return get<1>(a) < get<1>(b); });

        return { createEOSData(moles, temperature, pressure, get<0>(*min), get<1>(*min)) };
    }

    // ===== T,x Flash
    Phases EOSPengRobinson::flashTx(double temperature, double vaporFraction, double moles) const
    {
        using std::get;

        // ===== First, calculate the saturation pressure at the specified pressure.
        auto pressure = computeSaturationPressure(temperature, m_A, m_B, m_vaporPressureFunction);

        // ===== Calculate compressibility factors and fugacity coefficients at saturation conditions.
        auto A = m_A(temperature, pressure);
        auto B = m_B(temperature, pressure);

        auto z_phi = computeCompressibilityAndFugacity(A, B);
        auto z_v   = get<0>(z_phi[0]);
        auto z_l   = get<0>(z_phi[1]);
        auto phi_v = get<1>(z_phi[0]);
        auto phi_l = get<1>(z_phi[1]);

        // ===== If the specified vapor fraction is 1.0 (or higher), the fluid is a saturated vapor.
        if (vaporFraction >= 1.0) return { createEOSData(vaporFraction * moles, temperature, pressure, z_v, phi_v) };

        // ===== If the specified vapor fraction is 0.0 (or lower), the fluid is a saturated liquid.
        if (vaporFraction <= 0.0) return { createEOSData((1 - vaporFraction) * moles, temperature, pressure, z_l, phi_l) };

        // ===== If the vapor fraction is between 0.0 and 1.0, the fluid is two-phase.
        return { createEOSData(vaporFraction * moles, temperature, pressure, z_v, phi_v),
                 createEOSData((1 - vaporFraction) * moles, temperature, pressure, z_l, phi_l) };
    }

    // ===== P,x Flash
    Phases EOSPengRobinson::flashPx(double pressure, double vaporFraction, double moles) const
    {
        using std::get;

        // ===== First, calculate the saturation temperature at the specified pressure.
        auto temperature = computeSaturationTemperature(pressure, m_criticalTemperature, m_A, m_B, m_vaporPressureFunction);

        // ===== Calculate compressibility factors and fugacity coefficients at saturation conditions.
        auto A = m_A(temperature, pressure);
        auto B = m_B(temperature, pressure);

        auto z_phi = computeCompressibilityAndFugacity(A, B);
        auto z_v   = get<0>(z_phi[0]);
        auto z_l   = get<0>(z_phi[1]);
        auto phi_v = get<1>(z_phi[0]);
        auto phi_l = get<1>(z_phi[1]);

        // ===== If the specified vapor fraction is 1.0 (or higher), the fluid is a saturated vapor.
        if (vaporFraction >= 1.0) return { createEOSData(vaporFraction * moles, temperature, pressure, z_v, phi_v) };

        // ===== If the specified vapor fraction is 0.0 (or lower), the fluid is a saturated liquid.
        if (vaporFraction <= 0.0) return { createEOSData((1 - vaporFraction) * moles, temperature, pressure, z_l, phi_l) };

        // ===== If the vapor fraction is between 0.0 and 1.0, the fluid is two-phase.
        return { createEOSData(vaporFraction * moles, temperature, pressure, z_v, phi_v),
                 createEOSData((1 - vaporFraction) * moles, temperature, pressure, z_l, phi_l) };
    }

    // ===== P,H Flash
    Phases EOSPengRobinson::flashPH(double pressure, double enthalpy, double moles) const
    {
        using std::get;

        // ===== First, calculate the saturation temperature at the specified pressure.
        auto temperature = computeSaturationTemperature(pressure, m_criticalTemperature, m_A, m_B, m_vaporPressureFunction);

        // ===== Calculate compressibility factors and fugacity coefficients at saturation conditions.
        auto A     = m_A(temperature, pressure);
        auto B     = m_B(temperature, pressure);
        auto alpha = m_alpha(temperature);

        auto z_phi = computeCompressibilityAndFugacity(A, B);
        auto z_v   = get<0>(z_phi[0]);
        auto z_l   = get<0>(z_phi[1]);
        auto phi_v = get<1>(z_phi[0]);
        auto phi_l = get<1>(z_phi[1]);

        // ===== Then, calculate the enthalpy for the vapor and the liquid at saturation conditions.
        auto h_v = ::enthalpy(temperature, m_criticalTemperature, z_v, A, B, alpha, m_kappa, m_idealGasCpFunction);
        auto h_l = ::enthalpy(temperature, m_criticalTemperature, z_l, A, B, alpha, m_kappa, m_idealGasCpFunction);

        // ===== If the specified enthalpy is lower than the saturated liquid entropy, the fluid is a compressed liquid.
        if (enthalpy < h_l) {
            // ===== Define objective function
            auto f = [&](double t) {
                auto A     = m_A(t, pressure);
                auto B     = m_B(t, pressure);
                auto alpha = m_alpha(t);
                auto z     = computeCompressibilityAndFugacity(A, B);
                auto min   = std::min_element(z.begin(), z.end(), [](const auto& a, const auto& b) { return get<0>(a) < get<0>(b); });
                auto z_l   = get<0>(*min);
                return ::enthalpy(t, m_criticalTemperature, z_l, A, B, alpha, m_kappa, m_idealGasCpFunction) - enthalpy;
            };

            auto temp = PCProps::Numerics::ridders(f, 0, temperature);
            return { flashPT(pressure, temp, moles) };
        }

        // ===== If the specified enthalpy is higher than the saturated vapor entropy, the fluid is superheated vapor.
        if (enthalpy > h_v) {
            // ===== Define objective function
            auto f = [&](double t) {
                auto A     = m_A(t, pressure);
                auto B     = m_B(t, pressure);
                auto alpha = m_alpha(t);
                auto z     = computeCompressibilityAndFugacity(A, B);
                auto max   = std::max_element(z.begin(), z.end(), [](const auto& a, const auto& b) { return get<0>(a) < get<0>(b); });
                auto z_v   = get<0>(*max);
                return ::enthalpy(t, m_criticalTemperature, z_v, A, B, alpha, m_kappa, m_idealGasCpFunction) - enthalpy;
            };

            // ===== Determine the interval in which to find the root.
            auto diff = m_criticalTemperature - temperature;
            auto t1   = temperature;
            auto t2   = temperature + diff;
            while (true) {
                if (f(t1) * f(t2) < 0.0) break;
                t1 = t2;
                t2 = t1 + diff;
            }

            auto temp = PCProps::Numerics::ridders(f, t1, t2);
            return { flashPT(pressure, temp, moles) };
        }

        // ===== If the fluid is not a compressed liquid nor a superheated vapor, the fluid is two-phase.
        auto vaporFraction = (h_l - enthalpy) / (h_l - h_v);
        return { createEOSData(vaporFraction * moles, temperature, pressure, z_v, phi_v),
                 createEOSData((1 - vaporFraction) * moles, temperature, pressure, z_l, phi_l) };
    }

    // ===== P,S Flash
    Phases EOSPengRobinson::flashPS(double pressure, double entropy, double moles) const
    {
        using std::get;

        // ===== First, calculate the saturation temperature at the specified pressure.
        auto temperature = computeSaturationTemperature(pressure, m_criticalTemperature, m_A, m_B, m_vaporPressureFunction);

        // ===== Calculate compressibility factors and fugacity coefficients at saturation conditions.
        auto A     = m_A(temperature, pressure);
        auto B     = m_B(temperature, pressure);
        auto alpha = m_alpha(temperature);

        auto z_phi = computeCompressibilityAndFugacity(A, B);
        auto z_v   = get<0>(z_phi[0]);
        auto z_l   = get<0>(z_phi[1]);
        auto phi_v = get<1>(z_phi[0]);
        auto phi_l = get<1>(z_phi[1]);

        // ===== Then, calculate the entropy for the vapor and the liquid at saturation conditions.
        auto s_v = ::entropy(temperature, m_criticalTemperature, pressure, z_v, A, B, alpha, m_kappa, m_idealGasCpFunction);
        auto s_l = ::entropy(temperature, m_criticalTemperature, pressure, z_l, A, B, alpha, m_kappa, m_idealGasCpFunction);

        // ===== If the specified entropy is lower than the saturated liquid entropy, the fluid is a compressed liquid.
        if (entropy < s_l) {
            // ===== Define objective function
            auto f = [&](double t) {
                auto A     = m_A(t, pressure);
                auto B     = m_B(t, pressure);
                auto alpha = m_alpha(t);
                auto z     = computeCompressibilityAndFugacity(A, B);
                auto min   = std::min_element(z.begin(), z.end(), [](const auto& a, const auto& b) { return get<0>(a) < get<0>(b); });
                auto z_l   = get<0>(*min);
                return ::entropy(t, m_criticalTemperature, pressure, z_l, A, B, alpha, m_kappa, m_idealGasCpFunction) - entropy;
            };

            auto temp = PCProps::Numerics::ridders(f, 1, temperature);
            return { flashPT(pressure, temp, moles) };
        }

        // ===== If the specified entropy is higher than the saturated vapor entropy, the fluid is superheated vapor.
        if (entropy > s_v) {
            // ===== Define objective function
            auto f = [&](double t) {
                auto A     = m_A(t, pressure);
                auto B     = m_B(t, pressure);
                auto alpha = m_alpha(t);
                auto z     = computeCompressibilityAndFugacity(A, B);
                auto max   = std::max_element(z.begin(), z.end(), [](const auto& a, const auto& b) { return get<0>(a) < get<0>(b); });
                auto z_v   = get<0>(*max);
                return ::entropy(t, m_criticalTemperature, pressure, z_v, A, B, alpha, m_kappa, m_idealGasCpFunction) - entropy;
            };

            // ===== Determine the interval in which to find the root.
            auto diff = m_criticalTemperature - temperature;
            auto t1   = temperature;
            auto t2   = temperature + diff;
            while (true) {
                if (f(t1) * f(t2) < 0.0) break;
                t1 = t2;
                t2 = t1 + diff;
            }

            auto temp = PCProps::Numerics::ridders(f, t1, t2);
            return { flashPT(pressure, temp, moles) };
        }

        // ===== If the fluid is not a compressed liquid nor a superheated vapor, the fluid is two-phase.
        auto vaporFraction = (s_l - entropy) / (s_l - s_v);
        return { createEOSData(vaporFraction * moles, temperature, pressure, z_v, phi_v),
                 createEOSData((1 - vaporFraction) * moles, temperature, pressure, z_l, phi_l) };
    }
}    // namespace PCProps::EquationOfState