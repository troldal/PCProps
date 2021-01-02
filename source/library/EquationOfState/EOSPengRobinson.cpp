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


}    // namespace

namespace PCProps::EquationOfState
{
    class EOSPengRobinson::impl
    {
        // ===== Basic fluid properties
        double m_criticalTemperature {};
        double m_criticalPressure {};
        double m_acentricFactor {};
        double m_molecularWeight {};

        // ===== Calculated constants
        double m_b {};
        double m_ac {};

        // ===== User-supplied correlations
        std::function<double(double)>         m_vaporPressureFunction {};
        std::function<double(double)>         m_idealGasCpFunction {};
        std::function<double(double, double)> m_idealGasCpIntegralFunction {};
        std::function<double(double, double)> m_idealGasCpOverTemperatureIntegralFunction {};

        /**
         * @brief
         * @param temperature
         * @return
         */
        inline double a(double temperature) const
        {
            return m_ac * pow(1 + kappa() * (1 - sqrt(temperature / criticalTemperature())), 2);
        }

        /**
         * @brief
         * @return
         */
        inline double b() const
        {
            return m_b;
        }

        /**
         * @brief
         * @param temperature
         * @param pressure
         * @return
         */
        inline double A(double temperature, double pressure) const
        {
            return a(temperature) * pressure / pow(8.31446261815324, 2) / pow(temperature, 2);
        }

        /**
         * @brief
         * @param temperature
         * @param pressure
         * @return
         */
        inline double B(double temperature, double pressure) const
        {
            return b() * pressure / (8.31446261815324 * temperature);
        }

        /**
         * @brief
         * @param temperature
         * @return
         */
        inline double alpha(double temperature) const
        {
            return pow(1 + kappa() * (1 - sqrt(temperature / criticalTemperature())), 2);
        }

        /**
         * @brief
         * @return
         */
        inline double kappa() const
        {
            return 0.37464 + 1.54226 * m_acentricFactor - 0.26992 * pow(m_acentricFactor, 2);
        }

        /**
         * @brief
         * @param temperature
         * @return
         */
        inline double vaporPressure(double temperature) const
        {
            return m_vaporPressureFunction(temperature);
        }

        /**
         * @brief
         * @param temperature
         * @return
         */
        inline double idealGasCp(double temperature) const
        {
            return m_idealGasCpFunction(temperature);
        }

        /**
         * @brief
         * @param t1
         * @param t2
         * @return
         */
        inline double idealGasCpIntegral(double t1, double t2) const
        {
            return m_idealGasCpIntegralFunction(t1, t2);
        }

        /**
         * @brief
         * @param t1
         * @param t2
         * @return
         */
        inline double idealGasCpOverTemperatureIntegral(double t1, double t2)
        {
            return m_idealGasCpOverTemperatureIntegralFunction(t1, t2);
        }

        /**
         * @brief
         * @param temperature
         * @param pressure
         * @param compressibility
         * @return
         */
        inline double computeFugacityCoefficient(double temperature, double pressure, double compressibility) const
        {
            using std::exp;
            using std::log;
            using std::sqrt;

            auto coeffA = A(temperature, pressure);
            auto coeffB = B(temperature, pressure);

            return exp(
                compressibility - 1 - log(compressibility - coeffB) -
                coeffA / (coeffB * sqrt(8)) * log((compressibility + (1 + sqrt(2)) * coeffB) / (compressibility + (1 - sqrt(2)) * coeffB)));
        }

        /**
         * @brief
         * @param temperature
         * @param pressure
         * @return
         */
        std::vector<double> computeCompressibilityFactors(double temperature, double pressure) const
        {
            using std::acos;
            using std::cbrt;
            using std::cos;
            using std::sqrt;

            auto coeffA = A(temperature, pressure);
            auto coeffB = B(temperature, pressure);

            // ===== Compute the coefficients for solving Peng Robinson with respect to Z.
            auto a_0 = -(coeffA * coeffB - pow(coeffB, 2) - pow(coeffB, 3));
            auto a_1 = (coeffA - 3 * pow(coeffB, 2) - 2 * coeffB);
            auto a_2 = -(1 - coeffB);

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
         * @param temperature
         * @return
         */
        inline double idealGasEnthalpy(double temperature) const
        {
            using PCProps::Numerics::integrate;
            return integrate(m_idealGasCpFunction, PCProps::Globals::STANDARD_T, temperature);
        }

        /**
         * @brief
         * @param temperature
         * @param pressure
         * @param compressibility
         * @return
         */
        inline double enthalpyDeparture(double temperature, double pressure, double compressibility) const
        {
            using std::log;
            using std::sqrt;

            auto coeffA = A(temperature, pressure);
            auto coeffB = B(temperature, pressure);

            return (compressibility - 1.0 -
                    coeffA / (coeffB * sqrt(8)) * (1 + kappa() * sqrt(temperature / criticalTemperature()) / sqrt(alpha(temperature))) *
                        log((compressibility + (1 + sqrt(2)) * coeffB) / (compressibility + (1 - sqrt(2)) * coeffB))) *
                   PCProps::Globals::R_CONST * temperature;
        }

        /**
         * @brief
         * @param t
         * @param p
         * @return
         */
        inline double idealGasEntropy(double t, double p) const
        {
            using PCProps::Numerics::integrate;
            return integrate([&](double temp) { return m_idealGasCpFunction(temp) / temp; }, PCProps::Globals::STANDARD_T, t) -
                   PCProps::Globals::R_CONST * log(p / PCProps::Globals::STANDARD_P);
        }

        /**
         * @brief
         * @param t
         * @param p
         * @param z
         * @return
         */
        inline double entropyDeparture(double t, double p, double z) const
        {
            using std::log;
            using std::sqrt;

            auto coeffA = A(t, p);
            auto coeffB = B(t, p);

            return (log(z - coeffB) - coeffA / (coeffB * sqrt(8)) * kappa() * sqrt(t / criticalTemperature()) / sqrt(alpha(t)) *
                                          log((z + (1 + sqrt(2)) * coeffB) / (z + (1 - sqrt(2)) * coeffB))) *
                   PCProps::Globals::R_CONST;
        }

    public:
        /**
         * @brief
         * @param criticalTemperature
         * @param criticalPressure
         * @param acentricFactor
         * @param molecularWeight
         * @param vaporPressureFunction
         * @param idealGasCpFunction
         */
        impl(
            double                               criticalTemperature,
            double                               criticalPressure,
            double                               acentricFactor,
            double                               molecularWeight,
            const std::function<double(double)>& vaporPressureFunction,
            const std::function<double(double)>& idealGasCpFunction)
            : m_criticalTemperature(criticalTemperature),
              m_criticalPressure(criticalPressure),
              m_acentricFactor(acentricFactor),
              m_molecularWeight(molecularWeight),
              m_b(0.07779607 * 8.31446261815324 * criticalTemperature / criticalPressure),
              m_ac(0.45723553 * pow(8.31446261815324, 2) * pow(criticalTemperature, 2) / criticalPressure),
              m_vaporPressureFunction(vaporPressureFunction),
              m_idealGasCpFunction(idealGasCpFunction)
        {}

        /**
         * @brief
         * @return
         */
        inline double criticalTemperature() const
        {
            return m_criticalTemperature;
        }

        /**
         * @brief
         * @return
         */
        inline double criticalPressure() const
        {
            return m_criticalPressure;
        }

        /**
         * @brief
         * @return
         */
        inline double molecularWeight() const
        {
            return m_molecularWeight;
        }

        /**
         * @brief
         * @param moles
         * @param temperature
         * @param pressure
         * @param z
         * @param phi
         * @return
         */
        PhaseData createEOSData(double moles, double temperature, double pressure, double z, double phi) const
        {
            using std::get;
            PhaseData result;

            get<Moles>(result)           = moles;
            get<MolecularWeight>(result) = m_molecularWeight;
            get<Temperature>(result)     = temperature;
            get<Pressure>(result)        = pressure;
            get<Volume>(result)          = z * PCProps::Globals::R_CONST * temperature / pressure;
            get<Fugacity>(result)        = phi * pressure;
            get<Compressibility>(result) = z;
            get<Enthalpy>(result)        = enthalpy(temperature, pressure, z);
            get<Entropy>(result)         = entropy(temperature, pressure, z);
            get<InternalEnergy>(result)  = get<Enthalpy>(result) - pressure * get<Volume>(result);
            get<GibbsEnergy>(result)     = get<Enthalpy>(result) - temperature * get<Entropy>(result);
            get<HelmholzEnergy>(result)  = get<InternalEnergy>(result) - temperature * get<Entropy>(result);

            return result;
        }

        /**
         * @brief
         * @param temperature
         * @param pressure
         * @return
         */
        std::vector<std::pair<double, double>> computeCompressibilityAndFugacity(double temperature, double pressure) const
        {
            auto zs = computeCompressibilityFactors(temperature, pressure);

            std::vector<std::pair<double, double>> result;
            for (const auto& z : zs) result.emplace_back(std::make_pair(z, computeFugacityCoefficient(temperature, pressure, z)));

            return result;
        }

        /**
         * @brief
         * @param temperature
         * @return
         */
        inline double computeSaturationPressure(double temperature) const
        {
            using std::get;
            auto f = [&](double p) {
                auto phi   = computeCompressibilityAndFugacity(A(temperature, p), B(temperature, p));
                auto phi_v = get<1>(phi[0]);
                auto phi_l = get<1>(phi[1]);
                return (phi_l - phi_v) * p;
            };

            return PCProps::Numerics::ridders(f, vaporPressure(temperature) * 0.8, vaporPressure(temperature) * 1.2);
        }

        /**
         * @brief
         * @param pressure
         * @return
         */
        inline double computeSaturationTemperature(double pressure) const
        {
            using std::get;
            auto f = [&](double t) {
                auto phi   = computeCompressibilityAndFugacity(t, pressure);
                auto phi_v = get<1>(phi[0]);
                auto phi_l = get<1>(phi[1]);
                return (phi_l - phi_v) * pressure;
            };

            auto g = [&](double t) { return vaporPressure(t) - pressure; };

            auto guess = PCProps::Numerics::ridders(g, 0, criticalTemperature());
            return PCProps::Numerics::ridders(f, guess * 0.8, guess * 1.2);
        }

        /**
         * @brief
         * @param temperature
         * @param pressure
         * @param compressibility
         * @return
         */
        inline double enthalpy(double temperature, double pressure, double compressibility) const
        {
            return idealGasEnthalpy(temperature) + enthalpyDeparture(temperature, pressure, compressibility);
        }

        /**
         * @brief
         * @param temperature
         * @param pressure
         * @param compressibility
         * @return
         */
        inline double entropy(double temperature, double pressure, double compressibility) const
        {
            return idealGasEntropy(temperature, pressure) + entropyDeparture(temperature, pressure, compressibility);
        }
    };

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
        : m_impl(std::make_unique<
                 impl>(criticalTemperature, criticalPressure, acentricFactor, molecularWeight, vaporPressureFunction, idealGasCpFunction))
    {}

    // ===== Copy constructor
    EOSPengRobinson::EOSPengRobinson(const EOSPengRobinson& other) : m_impl(std::make_unique<impl>(*other.m_impl)) {};

    // ===== Move constructor
    EOSPengRobinson::EOSPengRobinson(EOSPengRobinson&& other) noexcept = default;

    // ===== Destructor
    EOSPengRobinson::~EOSPengRobinson() = default;

    // ===== Copy assignment operator
    EOSPengRobinson& EOSPengRobinson::operator=(const EOSPengRobinson& other)
    {
        EOSPengRobinson copy = other;
        *this                = std::move(copy);
        return *this;
    };

    // ===== Move assignment operator
    EOSPengRobinson& EOSPengRobinson::operator=(EOSPengRobinson&& other) noexcept = default;

    // ===== P,T Flash
    Phases EOSPengRobinson::flashPT(double pressure, double temperature, double moles) const
    {
        using std::get;

        // ===== Compute compressibility factors and fugacity coefficients at given T and P.
        auto z_phi = m_impl->computeCompressibilityAndFugacity(temperature, pressure);

        // ===== Identify the Z/Phi with the lowest fugacity coefficient (which is the most stable)
        auto min = std::min_element(z_phi.begin(), z_phi.end(), [](const auto& a, const auto& b) { return get<1>(a) < get<1>(b); });

        return { m_impl->createEOSData(moles, temperature, pressure, get<0>(*min), get<1>(*min)) };
    }

    // ===== T,x Flash
    Phases EOSPengRobinson::flashTx(double temperature, double vaporFraction, double moles) const
    {
        using std::get;

        // ===== First, calculate the saturation pressure at the specified pressure.
        auto pressure = m_impl->computeSaturationPressure(temperature);

        auto z_phi = m_impl->computeCompressibilityAndFugacity(temperature, pressure);
        auto z_v   = get<0>(z_phi[0]);
        auto z_l   = get<0>(z_phi[1]);
        auto phi_v = get<1>(z_phi[0]);
        auto phi_l = get<1>(z_phi[1]);

        // ===== If the specified vapor fraction is 1.0 (or higher), the fluid is a saturated vapor.
        if (vaporFraction >= 1.0) return { m_impl->createEOSData(vaporFraction * moles, temperature, pressure, z_v, phi_v) };

        // ===== If the specified vapor fraction is 0.0 (or lower), the fluid is a saturated liquid.
        if (vaporFraction <= 0.0) return { m_impl->createEOSData((1 - vaporFraction) * moles, temperature, pressure, z_l, phi_l) };

        // ===== If the vapor fraction is between 0.0 and 1.0, the fluid is two-phase.
        return { m_impl->createEOSData(vaporFraction * moles, temperature, pressure, z_v, phi_v),
                 m_impl->createEOSData((1 - vaporFraction) * moles, temperature, pressure, z_l, phi_l) };
    }

    // ===== P,x Flash
    Phases EOSPengRobinson::flashPx(double pressure, double vaporFraction, double moles) const
    {
        using std::get;

        // ===== First, calculate the saturation temperature at the specified pressure.
        auto temperature = m_impl->computeSaturationTemperature(pressure);

        auto z_phi = m_impl->computeCompressibilityAndFugacity(temperature, pressure);
        auto z_v   = get<0>(z_phi[0]);
        auto z_l   = get<0>(z_phi[1]);
        auto phi_v = get<1>(z_phi[0]);
        auto phi_l = get<1>(z_phi[1]);

        // ===== If the specified vapor fraction is 1.0 (or higher), the fluid is a saturated vapor.
        if (vaporFraction >= 1.0) return { m_impl->createEOSData(vaporFraction * moles, temperature, pressure, z_v, phi_v) };

        // ===== If the specified vapor fraction is 0.0 (or lower), the fluid is a saturated liquid.
        if (vaporFraction <= 0.0) return { m_impl->createEOSData((1 - vaporFraction) * moles, temperature, pressure, z_l, phi_l) };

        // ===== If the vapor fraction is between 0.0 and 1.0, the fluid is two-phase.
        return { m_impl->createEOSData(vaporFraction * moles, temperature, pressure, z_v, phi_v),
                 m_impl->createEOSData((1 - vaporFraction) * moles, temperature, pressure, z_l, phi_l) };
    }

    // ===== P,H Flash
    Phases EOSPengRobinson::flashPH(double pressure, double enthalpy, double moles) const
    {
        using std::get;

        // ===== First, calculate the saturation temperature at the specified pressure.
        auto temperature = m_impl->computeSaturationTemperature(pressure);

        auto z_phi = m_impl->computeCompressibilityAndFugacity(temperature, pressure);
        auto z_v   = get<0>(z_phi[0]);
        auto z_l   = get<0>(z_phi[1]);
        auto phi_v = get<1>(z_phi[0]);
        auto phi_l = get<1>(z_phi[1]);

        // ===== Then, calculate the enthalpy for the vapor and the liquid at saturation conditions.
        auto h_v = m_impl->enthalpy(temperature, pressure, z_v);
        auto h_l = m_impl->enthalpy(temperature, pressure, z_l);

        // ===== If the specified enthalpy is lower than the saturated liquid entropy, the fluid is a compressed liquid.
        if (enthalpy < h_l) {
            // ===== Define objective function
            auto f = [&](double t) {
                auto z   = m_impl->computeCompressibilityAndFugacity(temperature, pressure);
                auto min = std::min_element(z.begin(), z.end(), [](const auto& a, const auto& b) { return get<0>(a) < get<0>(b); });
                auto z_l = get<0>(*min);
                return m_impl->enthalpy(t, pressure, z_l) - enthalpy;
            };

            auto temp = PCProps::Numerics::ridders(f, 0, temperature);
            return { flashPT(pressure, temp, moles) };
        }

        // ===== If the specified enthalpy is higher than the saturated vapor entropy, the fluid is superheated vapor.
        if (enthalpy > h_v) {
            // ===== Define objective function
            auto f = [&](double t) {
                auto z   = m_impl->computeCompressibilityAndFugacity(temperature, pressure);
                auto max = std::max_element(z.begin(), z.end(), [](const auto& a, const auto& b) { return get<0>(a) < get<0>(b); });
                auto z_v   = get<0>(*max);
                return m_impl->enthalpy(t, pressure, z_v) - enthalpy;
            };

            // ===== Determine the interval in which to find the root.
            auto diff = m_impl->criticalTemperature() - temperature;
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
        return { m_impl->createEOSData(vaporFraction * moles, temperature, pressure, z_v, phi_v),
                 m_impl->createEOSData((1 - vaporFraction) * moles, temperature, pressure, z_l, phi_l) };
    }

    // ===== P,S Flash
    Phases EOSPengRobinson::flashPS(double pressure, double entropy, double moles) const
    {
        using std::get;

        // ===== First, calculate the saturation temperature at the specified pressure.
        auto temperature = m_impl->computeSaturationTemperature(pressure);

        auto z_phi = m_impl->computeCompressibilityAndFugacity(temperature, pressure);
        auto z_v   = get<0>(z_phi[0]);
        auto z_l   = get<0>(z_phi[1]);
        auto phi_v = get<1>(z_phi[0]);
        auto phi_l = get<1>(z_phi[1]);

        // ===== Then, calculate the entropy for the vapor and the liquid at saturation conditions.
        auto s_v = m_impl->entropy(temperature, pressure, z_v);
        auto s_l = m_impl->entropy(temperature, pressure, z_l);

        // ===== If the specified entropy is lower than the saturated liquid entropy, the fluid is a compressed liquid.
        if (entropy < s_l) {
            // ===== Define objective function
            auto f = [&](double t) {
                auto z   = m_impl->computeCompressibilityAndFugacity(temperature, pressure);
                auto min = std::min_element(z.begin(), z.end(), [](const auto& a, const auto& b) { return get<0>(a) < get<0>(b); });
                auto z_l = get<0>(*min);
                return m_impl->entropy(t, pressure, z_l) - entropy;
            };

            auto temp = PCProps::Numerics::ridders(f, 1, temperature);
            return { flashPT(pressure, temp, moles) };
        }

        // ===== If the specified entropy is higher than the saturated vapor entropy, the fluid is superheated vapor.
        if (entropy > s_v) {
            // ===== Define objective function
            auto f = [&](double t) {
                auto z   = m_impl->computeCompressibilityAndFugacity(temperature, pressure);
                auto max = std::max_element(z.begin(), z.end(), [](const auto& a, const auto& b) { return get<0>(a) < get<0>(b); });
                auto z_v   = get<0>(*max);
                return m_impl->entropy(t, pressure, z_v) - entropy;
            };

            // ===== Determine the interval in which to find the root.
            auto diff = m_impl->criticalTemperature() - temperature;
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
        return { m_impl->createEOSData(vaporFraction * moles, temperature, pressure, z_v, phi_v),
                 m_impl->createEOSData((1 - vaporFraction) * moles, temperature, pressure, z_l, phi_l) };
    }
}    // namespace PCProps::EquationOfState