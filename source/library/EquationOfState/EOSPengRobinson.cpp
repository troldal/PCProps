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
#include <vector>

#include "EOSPengRobinson.hpp"
#include <PCConfig.hpp>
#include <Utilities/Integration.hpp>
#include <Utilities/RootFinding.hpp>

namespace PCProps::EquationOfState
{
    using PCProps::Globals::PI;
    using PCProps::Globals::R_CONST;
    using PCProps::Globals::STANDARD_P;
    using PCProps::Globals::STANDARD_T;

    class EOSPengRobinson::impl
    {
    private:
        // ===== Basic fluid properties
        double m_criticalTemperature {};
        double m_criticalPressure {};

        // ===== Calculated constants
        double m_ac {};
        double m_b {};
        double m_kappa {};

        // ===== User-supplied correlations
        std::function<double(double)> m_vaporPressureFunction {};
        std::function<double(double)> m_idealGasCpFunction {};
        std::function<double(double)> m_idealGasCpDerivativeFunction {};
        std::function<double(double)> m_idealGasCpIntegralFunction {};
        std::function<double(double)> m_idealGasCpOverTemperatureIntegralFunction {};

        /**
         * @brief Compute the 'a' coefficient for the Peng-robinson EOS.
         * @param temperature The temperature [K].
         * @return The 'a' coefficient.
         */
        inline double a(double temperature) const
        {
            return m_ac * pow(1 + m_kappa * (1 - sqrt(temperature / criticalTemperature())), 2);
        }

        /**
         * @brief Compute the 'A' dimensionless coefficient for the Peng-Robinson EOS.
         * @param temperature The temperature [K].
         * @param pressure The pressure [Pa]
         * @return The 'A' coefficient.
         */
        inline double A(double temperature, double pressure) const
        {
            using PCProps::Globals::R_CONST;
            return a(temperature) * pressure / pow(R_CONST, 2) / pow(temperature, 2);
        }

        /**
         * @brief Compute the 'B' dimensionless coefficient for the Peng-Robinson EOS.
         * @param temperature The temperature [K].
         * @param pressure The pressure [Pa]
         * @return The 'B' coefficient.
         */
        inline double B(double temperature, double pressure) const
        {
            using PCProps::Globals::R_CONST;
            return m_b * pressure / (R_CONST * temperature);
        }

        /**
         * @brief Compute the 'alpha' coefficient for the Peng-Robinson EOS.
         * @param temperature The temperature [K].
         * @return The 'alpha' coefficient.
         */
        inline double alpha(double temperature) const
        {
            return pow(1 + m_kappa * (1 - sqrt(temperature / criticalTemperature())), 2);
        }

        /**
         * @brief Compute the vapor pressure at the given T, using the user supplied correlation.
         * @param temperature The temperature [K].
         * @return The vapor pressure [Pa].
         * @todo This may not be needed. Consider removing it.
         */
        inline double vaporPressure(double temperature) const
        {
            return m_vaporPressureFunction(temperature);
        }

        /**
         * @brief Compute the compressibility factors for all phases of the fluid, at the given T and P.
         * @details The compressibility factors are found by solving Peng-Robinson as a cubic equation, with
         * respect to Z. Only real roots are considered. There may be either one or tree real roots. In the latter
         * case, two roots may be equal. If there are three distinct roots, the middle one is discarded, as it has
         * no physical meaning. The roots are found analytically, and it is made sure that only unique roots
         * larger than 0.0 are returned.
         * @param temperature The temperature [K].
         * @param pressure The pressure [Pa].
         * @return A std::vector with the compressibility factors.
         */
        std::vector<double> computeCompressibilityFactors(double temperature, double pressure) const
        {
            using std::acos;
            using std::cbrt;
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

                // TODO (troldal): The middle root should be removed, instead of copying the min and max to a separate vector.
                std::vector<double> all_roots { m * cos(theta) - a_2 / 3,
                                                m * cos(theta + 2 * PI / 3) - a_2 / 3,
                                                m * cos(theta + 4 * PI / 3) - a_2 / 3 };

                // ===== Discard the middle root, as it has no physical meaning.
                std::vector<double> roots { *std::max_element(all_roots.begin(), all_roots.end()),
                                            *std::min_element(all_roots.begin(), all_roots.end()) };

                // ===== If any of the roots are negative, delete them.
                if (roots[1] <= 0.0) roots.pop_back();
                if (roots[0] <= 0.0) roots.erase(roots.begin());

                // ===== If the two roots are equal, delete the last one.
                if (roots.size() == 2 && roots[0] == roots[1]) roots.pop_back();

                return roots;
            }

            // ===== If R > 0, there is one real root
            auto P = cbrt(-q / 2.0 + sqrt(R));
            auto Q = cbrt(-q / 2.0 - sqrt(R));

            return { P + Q - a_2 / 3.0 };
        }

        /**
         * @brief Compute the fugacity coefficient of the fluid with the given compressibility factor and the given T and P.
         * @param temperature The temperature [K].
         * @param pressure The pressure [Pa].
         * @param compressibility The compressibility [-].
         * @return The fugacity coefficient (phi) at the given P, T and Z.
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
         * @brief Compute the ideal gas enthalpy at the given T, relative the standard state.
         * @param temperature The temperature [K].
         * @return The ideal gas enthalpy [J/mol]
         */
        inline double idealGasEnthalpy(double temperature) const
        {
            using PCProps::Globals::STANDARD_T;
            return (m_idealGasCpIntegralFunction(temperature) - m_idealGasCpIntegralFunction(STANDARD_T));

            //            using PCProps::Numerics::integrate;
            //            return integrate(m_idealGasCpFunction, PCProps::Globals::STANDARD_T, temperature);
        }

        /**
         * @brief Compute the enthalpy departure relative to standard (ideal gas) state, at the given T and P.
         * @param temperature The temperature [K].
         * @param pressure The pressure [Pa].
         * @param compressibility The compressibility factor [-].
         * @return The enthalpy departure (H-H^ig) for the fluid [J/mol]
         */
        inline double enthalpyDeparture(double temperature, double pressure, double compressibility) const
        {
            using std::log;
            using std::log;
            using std::sqrt;

            auto coeffA = A(temperature, pressure);
            auto coeffB = B(temperature, pressure);

            return (compressibility - 1.0 -
                    coeffA / (coeffB * sqrt(8)) * (1 + m_kappa * sqrt(temperature / criticalTemperature()) / sqrt(alpha(temperature))) *
                        log((compressibility + (1 + sqrt(2)) * coeffB) / (compressibility + (1 - sqrt(2)) * coeffB))) *
                   R_CONST * temperature;
        }

        /**
         * @brief Compute the ideal gas entropy at the given T and P, relative the standard state.
         * @param t The temperature [K].
         * @param p The pressure [Pa].
         * @return The ideal gas entropy [J/mol-K]
         */
        inline double idealGasEntropy(double t, double p) const
        {
            using PCProps::Globals::R_CONST;
            using PCProps::Globals::STANDARD_P;
            using PCProps::Globals::STANDARD_T;

            return (m_idealGasCpOverTemperatureIntegralFunction(t) - m_idealGasCpOverTemperatureIntegralFunction(STANDARD_T)) - R_CONST * log(p / STANDARD_P);
        }

        /**
         * @brief Compute the entropy departure relative to standard (ideal gas) state, at the given T and P.
         * @param t The temperature [K].
         * @param p The pressure [Pa]
         * @param z The compressibility factor for the fluid [-]
         * @return The entropy departure (S-S^ig) for the fluid [J/mol-K]
         */
        inline double entropyDeparture(double t, double p, double z) const
        {
            using PCProps::Globals::R_CONST;
            using std::log;
            using std::sqrt;

            auto coeffA = A(t, p);
            auto coeffB = B(t, p);

            return (log(z - coeffB) -
                    coeffA / (coeffB * sqrt(8)) * m_kappa * sqrt(t / criticalTemperature()) / sqrt(alpha(t)) * log((z + (1 + sqrt(2)) * coeffB) / (z + (1 - sqrt(2)) * coeffB))) *
                   R_CONST;
        }

    public:
        /**
         * @brief Constructor.
         * @param criticalTemperature The critical temperature [K].
         * @param criticalPressure The critical pressure [Pa]
         * @param acentricFactor The acentric factor [-]
         * @param molecularWeight The molecular weight [g/mol]
         * @param vaporPressureFunction Function object for computing the vapor pressure as a function of temperature.
         * @param idealGasCpFunction Function object for computing the Cp as a function of temperature.
         */
        impl(double criticalTemperature, double criticalPressure, double acentricFactor)
            : m_criticalTemperature(criticalTemperature),
              m_criticalPressure(criticalPressure),
              m_ac(0.45723553 * pow(PCProps::Globals::R_CONST, 2) * pow(criticalTemperature, 2) / criticalPressure),
              m_b(0.07779607 * PCProps::Globals::R_CONST * criticalTemperature / criticalPressure),
              m_kappa(0.37464 + 1.54226 * acentricFactor - 0.26992 * pow(acentricFactor, 2))
        {}

        /**
         * @brief
         * @param criticalTemperature
         * @param criticalPressure
         * @param acentricFactor
         */
        inline void setProperties(double criticalTemperature, double criticalPressure, double acentricFactor)
        {
            m_criticalTemperature = criticalTemperature;
            m_criticalPressure    = criticalPressure;
            m_ac                  = 0.45723553 * pow(PCProps::Globals::R_CONST, 2) * pow(criticalTemperature, 2) / criticalPressure;
            m_b                   = 0.07779607 * PCProps::Globals::R_CONST * criticalTemperature / criticalPressure;
            m_kappa               = 0.37464 + 1.54226 * acentricFactor - 0.26992 * pow(acentricFactor, 2);
        }

        /**
         * @brief
         * @param vaporPressureFunction
         */
        void setVaporPressureFunction(const std::function<double(double)>& vaporPressureFunction)
        {
            m_vaporPressureFunction = vaporPressureFunction;
        }

        /**
         * @brief
         * @param idealGasCpFunction
         */
        void setIdealGasCpFunction(const std::function<double(double)>& idealGasCpFunction)
        {
            m_idealGasCpFunction = idealGasCpFunction;
        }

        /**
         * @brief
         * @param idealGasCpDerivativeFunction
         */
        void setIdealGasCpDerivativeFunction(const std::function<double(double)>& idealGasCpDerivativeFunction)
        {
            m_idealGasCpDerivativeFunction = idealGasCpDerivativeFunction;
        }

        /**
         * @brief
         * @param idealGasCpIntegralFunction
         */
        void setIdealGasCpIntegralFunction(const std::function<double(double)>& idealGasCpIntegralFunction)
        {
            m_idealGasCpIntegralFunction = idealGasCpIntegralFunction;
        }

        /**
         * @brief
         * @param idealGasOverTIntegralFunction
         */
        void setIdealGasCpOverTIntegralFunction(const std::function<double(double)>& idealGasOverTIntegralFunction)
        {
            m_idealGasCpOverTemperatureIntegralFunction = idealGasOverTIntegralFunction;
        }

        /**
         * @brief Accessor for the critical temperature.
         * @return The critical temperature [K].
         */
        inline double criticalTemperature() const
        {
            return m_criticalTemperature;
        }

        /**
         * @brief Accessor for the critical pressure.
         * @return The critical pressure [Pa].
         */
        inline double criticalPressure() const
        {
            return m_criticalPressure;
        }

        /**
         * @brief Helper function for creating the std::tuple with the output data.
         * @param moleFraction The number of moles.
         * @param temperature The temperature [K].
         * @param pressure The pressure [Pa].
         * @param z The compressibility factor [-]
         * @param phi The fugacity coefficient [-].
         * @return A PhaseData object (aka std::tuple) with the output data.
         */
        PhaseData createEOSData(double moleFraction, double temperature, double pressure, double z, double phi) const
        {
            using PCProps::Globals::R_CONST;
            using std::get;

            PhaseData result;

            get<MolarFraction>(result) = moleFraction;
            //            get<MolecularWeight>(result) = m_molecularWeight;
            get<Temperature>(result)         = temperature;
            get<Pressure>(result)            = pressure;
            get<Volume>(result)              = z * R_CONST * temperature / pressure;
            get<FugacityCoefficient>(result) = phi;
            get<Fugacity>(result)            = phi * pressure;
            get<Compressibility>(result)     = z;
            get<Enthalpy>(result)            = computeEnthalpy(temperature, pressure, z);
            get<Entropy>(result)             = computeEntropy(temperature, pressure, z);
            get<InternalEnergy>(result)      = get<Enthalpy>(result) - pressure * get<Volume>(result);
            get<GibbsEnergy>(result)         = get<Enthalpy>(result) - temperature * get<Entropy>(result);
            get<HelmholzEnergy>(result)  = get<InternalEnergy>(result) - temperature * get<Entropy>(result);

            return result;
        }

        /**
         * @brief Compute the compressibilities and fugacity coefficients for the phases present at the given T and P.
         * @param temperature The temperature [K]
         * @param pressure The pressure [Pa].
         * @return A vector of std::pairs with the compressibility and fugacity coefficient for each phase.
         * @note A flash may yield Z and phi for two phases, even though only one phase is present. The individual
         * flash algorithm will discard the Z/Phi pair that is invalid.
         */
        std::vector<std::pair<double, double>> computeCompressibilityAndFugacity(double temperature, double pressure) const
        {
            // ===== Compute the compressibility factor(s)
            auto zs = computeCompressibilityFactors(temperature, pressure);

            // ===== Compute the fugacity coefficient for each Z and add the pair(s) to the results vector.
            std::vector<std::pair<double, double>> result;
            for (const auto& z : zs) result.emplace_back(std::make_pair(z, computeFugacityCoefficient(temperature, pressure, z)));

            return result;
        }

        /**
         * @brief Compute the fluid saturation pressure at the given T.
         * @param temperature The temperature [K].
         * @return The saturation pressure [Pa].
         */
        inline double computeSaturationPressure(double temperature) const
        {
            // TODO (troldal): I'm not sure this works as intended. What if the fluid is supercritical?
            using std::get;
            auto f = [&](double p) {
                auto phi = computeCompressibilityAndFugacity(temperature, p);
                if (phi.size() == 1) return get<1>(phi[0]) * p;

                auto phi_v = get<1>(phi[0]);
                auto phi_l = get<1>(phi[1]);
                return (phi_l - phi_v) * p;
            };

            return PCProps::Numerics::ridders(f, vaporPressure(temperature) * 0.8, vaporPressure(temperature) * 1.2);
        }

        /**
         * @brief Compute the fluid saturation temperature at the given P.
         * @param pressure The pressure [Pa].
         * @return The saturation temperature [K].
         */
        inline double computeSaturationTemperature(double pressure) const
        {
            // TODO (troldal): I'm not sure this works as intended. What if the fluid is supercritical?
            using std::get;
            auto f = [&](double t) {
                auto phi = computeCompressibilityAndFugacity(t, pressure);
                if (phi.size() == 1) return get<1>(phi[0]) * pressure;

                auto phi_v = get<1>(phi[0]);
                auto phi_l = get<1>(phi[1]);
                return (phi_l - phi_v) * pressure;
            };

            auto g = [&](double t) { return vaporPressure(t) - pressure; };

            auto guess = PCProps::Numerics::ridders(g, 0, criticalTemperature());
            return PCProps::Numerics::ridders(f, guess * 0.8, guess * 1.2);
        }

        /**
         * @brief Compute the enthalpy of the fluid, relative to the standard state.
         * @param temperature The temperature [K].
         * @param pressure The pressure [Pa].
         * @param compressibility The fluid compressibility [-].
         * @return The fluid enthalpy [J/mol].
         */
        inline double computeEnthalpy(double temperature, double pressure, double compressibility) const
        {
            return idealGasEnthalpy(temperature) + enthalpyDeparture(temperature, pressure, compressibility);
        }

        /**
         * @brief Compute the entropy of the fluid, relative to the standard state.
         * @param temperature The temperature [K].
         * @param pressure The pressure [Pa].
         * @param compressibility The fluid compressibility [-].
         * @return The fluid entropy [J/mol-K].
         */
        inline double computeEntropy(double temperature, double pressure, double compressibility) const
        {
            return idealGasEntropy(temperature, pressure) + entropyDeparture(temperature, pressure, compressibility);
        }
    };

    // ===== Constructor, default
    EOSPengRobinson::EOSPengRobinson() : m_impl(std::make_unique<impl>(0.0, 0.0, 0.0)) {};

    // ===== Constructor
    EOSPengRobinson::EOSPengRobinson(double criticalTemperature, double criticalPressure, double acentricFactor)
        : m_impl(std::make_unique<impl>(criticalTemperature, criticalPressure, acentricFactor))
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

    void EOSPengRobinson::setProperties(double criticalTemperature, double criticalPressure, double acentricFactor)
    {
        m_impl->setProperties(criticalTemperature, criticalPressure, acentricFactor);
    }

    void EOSPengRobinson::setVaporPressureFunction(const std::function<double(double)>& vaporPressureFunction)
    {
        m_impl->setVaporPressureFunction(vaporPressureFunction);
    }

    void EOSPengRobinson::setIdealGasCpFunction(const std::function<double(double)>& idealGasCpFunction)
    {
        m_impl->setIdealGasCpFunction(idealGasCpFunction);
    }

    void EOSPengRobinson::setIdealGasCpDerivativeFunction(const std::function<double(double)>& idealGasCpDerivativeFunction)
    {
        m_impl->setIdealGasCpDerivativeFunction(idealGasCpDerivativeFunction);
    }

    void EOSPengRobinson::setIdealGasCpIntegralFunction(const std::function<double(double)>& idealGasCpIntegralFunction)
    {
        m_impl->setIdealGasCpIntegralFunction(idealGasCpIntegralFunction);
    }

    void EOSPengRobinson::setIdealGasCpOverTIntegralFunction(const std::function<double(double)>& idealGasOverTIntegralFunction)
    {
        m_impl->setIdealGasCpOverTIntegralFunction(idealGasOverTIntegralFunction);
    }

    // ===== P,T Flash
    Phases EOSPengRobinson::flashPT(double pressure, double temperature) const
    {
        using std::get;

        // ===== Compute compressibility factors and fugacity coefficients at given T and P.
        auto z_phi = m_impl->computeCompressibilityAndFugacity(temperature, pressure);

        // ===== Identify the Z/Phi with the lowest fugacity coefficient (which is the most stable)
        auto min = std::min_element(z_phi.begin(), z_phi.end(), [](const auto& a, const auto& b) { return get<1>(a) < get<1>(b); });

        return { m_impl->createEOSData(1.0, temperature, pressure, get<0>(*min), get<1>(*min)) };
    }

    // ===== T,x Flash
    Phases EOSPengRobinson::flashTx(double temperature, double vaporFraction) const
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
        if (vaporFraction >= 1.0) return { m_impl->createEOSData(vaporFraction, temperature, pressure, z_v, phi_v) };

        // ===== If the specified vapor fraction is 0.0 (or lower), the fluid is a saturated liquid.
        if (vaporFraction <= 0.0) return { m_impl->createEOSData((1 - vaporFraction), temperature, pressure, z_l, phi_l) };

        // ===== If the vapor fraction is between 0.0 and 1.0, the fluid is two-phase.
        return { m_impl->createEOSData(vaporFraction, temperature, pressure, z_v, phi_v), m_impl->createEOSData((1 - vaporFraction), temperature, pressure, z_l, phi_l) };
    }

    // ===== P,x Flash
    Phases EOSPengRobinson::flashPx(double pressure, double vaporFraction) const
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
        if (vaporFraction >= 1.0) return { m_impl->createEOSData(vaporFraction, temperature, pressure, z_v, phi_v) };

        // ===== If the specified vapor fraction is 0.0 (or lower), the fluid is a saturated liquid.
        if (vaporFraction <= 0.0) return { m_impl->createEOSData((1 - vaporFraction), temperature, pressure, z_l, phi_l) };

        // ===== If the vapor fraction is between 0.0 and 1.0, the fluid is two-phase.
        return { m_impl->createEOSData(vaporFraction, temperature, pressure, z_v, phi_v), m_impl->createEOSData((1 - vaporFraction), temperature, pressure, z_l, phi_l) };
    }

    // ===== P,H Flash
    Phases EOSPengRobinson::flashPH(double pressure, double enthalpy) const
    {
        using std::get;

        if (pressure > m_impl->criticalPressure()) {
            // ===== Define objective function
            auto f = [&](double t) {
                auto z   = m_impl->computeCompressibilityAndFugacity(t, pressure);
                auto max = std::max_element(z.begin(), z.end(), [](const auto& a, const auto& b) { return get<0>(a) < get<0>(b); });
                auto z_v = get<0>(*max);
                return m_impl->computeEnthalpy(t, pressure, z_v) - enthalpy;
            };

            // ===== Determine the interval in which to find the root.
            auto t1 = 1.0;
            auto t2 = 1.0 + m_impl->criticalTemperature();
            while (true) {
                if (f(t1) * f(t2) < 0.0) break;
                t1 = t2;
                t2 = t1 + m_impl->criticalTemperature();
            }

            auto temp = PCProps::Numerics::ridders(f, t1, t2);
            return { flashPT(pressure, temp) };
        }

        // ===== First, calculate the saturation temperature at the specified pressure.
        auto temperature = m_impl->computeSaturationTemperature(pressure);

        auto z_phi = m_impl->computeCompressibilityAndFugacity(temperature, pressure);
        auto z_v   = get<0>(z_phi[0]);
        auto z_l   = get<0>(z_phi[1]);
        auto phi_v = get<1>(z_phi[0]);
        auto phi_l = get<1>(z_phi[1]);

        // ===== Then, calculate the enthalpy for the vapor and the liquid at saturation conditions.
        auto h_v = m_impl->computeEnthalpy(temperature, pressure, z_v);
        auto h_l = m_impl->computeEnthalpy(temperature, pressure, z_l);

        // ===== If the specified enthalpy is lower than the saturated liquid entropy, the fluid is a compressed liquid.
        if (enthalpy < h_l) {
            // ===== Define objective function
            auto f = [&](double t) {
                auto z   = m_impl->computeCompressibilityAndFugacity(t, pressure);
                auto min = std::min_element(z.begin(), z.end(), [](const auto& a, const auto& b) { return get<0>(a) < get<0>(b); });
                auto z_l = get<0>(*min);
                return m_impl->computeEnthalpy(t, pressure, z_l) - enthalpy;
            };

            auto temp = PCProps::Numerics::ridders(f, 0, temperature);
            return { flashPT(pressure, temp) };
        }

        // ===== If the specified enthalpy is higher than the saturated vapor entropy, the fluid is superheated vapor.
        if (enthalpy > h_v) {
            // ===== Define objective function
            auto f = [&](double t) {
                auto z   = m_impl->computeCompressibilityAndFugacity(t, pressure);
                auto max = std::max_element(z.begin(), z.end(), [](const auto& a, const auto& b) { return get<0>(a) < get<0>(b); });
                auto z_v   = get<0>(*max);
                return m_impl->computeEnthalpy(t, pressure, z_v) - enthalpy;
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
            return { flashPT(pressure, temp) };
        }

        // ===== If the fluid is not a compressed liquid nor a superheated vapor, the fluid is two-phase.
        auto vaporFraction = (h_l - enthalpy) / (h_l - h_v);
        return { m_impl->createEOSData(vaporFraction, temperature, pressure, z_v, phi_v), m_impl->createEOSData((1 - vaporFraction), temperature, pressure, z_l, phi_l) };
    }

    // ===== P,S Flash
    Phases EOSPengRobinson::flashPS(double pressure, double entropy) const
    {
        using std::get;

        if (pressure > m_impl->criticalPressure()) {
            // ===== Define objective function
            auto f = [&](double t) {
                auto z   = m_impl->computeCompressibilityAndFugacity(t, pressure);
                auto max = std::max_element(z.begin(), z.end(), [](const auto& a, const auto& b) { return get<0>(a) < get<0>(b); });
                auto z_v = get<0>(*max);
                return m_impl->computeEntropy(t, pressure, z_v) - entropy;
            };

            // ===== Determine the interval in which to find the root.
            auto t1 = 1.0;
            auto t2 = 1.0 + m_impl->criticalTemperature();
            while (true) {
                if (f(t1) * f(t2) < 0.0) break;
                t1 = t2;
                t2 = t1 + m_impl->criticalTemperature();
            }

            auto temp = PCProps::Numerics::ridders(f, t1, t2);
            return { flashPT(pressure, temp) };
        }

        // ===== First, calculate the saturation temperature at the specified pressure.
        auto temperature = m_impl->computeSaturationTemperature(pressure);

        auto z_phi = m_impl->computeCompressibilityAndFugacity(temperature, pressure);
        auto z_v   = get<0>(z_phi[0]);
        auto z_l   = get<0>(z_phi[1]);
        auto phi_v = get<1>(z_phi[0]);
        auto phi_l = get<1>(z_phi[1]);

        // ===== Then, calculate the entropy for the vapor and the liquid at saturation conditions.
        auto s_v = m_impl->computeEntropy(temperature, pressure, z_v);
        auto s_l = m_impl->computeEntropy(temperature, pressure, z_l);

        // ===== If the specified entropy is lower than the saturated liquid entropy, the fluid is a compressed liquid.
        if (entropy < s_l) {
            // ===== Define objective function
            auto f = [&](double t) {
                auto z   = m_impl->computeCompressibilityAndFugacity(t, pressure);
                auto min = std::min_element(z.begin(), z.end(), [](const auto& a, const auto& b) { return get<0>(a) < get<0>(b); });
                auto z_l = get<0>(*min);
                return m_impl->computeEntropy(t, pressure, z_l) - entropy;
            };

            auto temp = PCProps::Numerics::ridders(f, 1, temperature);
            return { flashPT(pressure, temp) };
        }

        // ===== If the specified entropy is higher than the saturated vapor entropy, the fluid is superheated vapor.
        if (entropy > s_v) {
            // ===== Define objective function
            auto f = [&](double t) {
                auto z   = m_impl->computeCompressibilityAndFugacity(t, pressure);
                auto max = std::max_element(z.begin(), z.end(), [](const auto& a, const auto& b) { return get<0>(a) < get<0>(b); });
                auto z_v   = get<0>(*max);
                return m_impl->computeEntropy(t, pressure, z_v) - entropy;
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
            return { flashPT(pressure, temp) };
        }

        // ===== If the fluid is not a compressed liquid nor a superheated vapor, the fluid is two-phase.
        auto vaporFraction = (s_l - entropy) / (s_l - s_v);
        return { m_impl->createEOSData(vaporFraction, temperature, pressure, z_v, phi_v), m_impl->createEOSData((1 - vaporFraction), temperature, pressure, z_l, phi_l) };
    }

    // ===== Saturation pressure at given temperature
    double EOSPengRobinson::saturationPressure(double temperature) const
    {
        return m_impl->computeSaturationPressure(temperature);
    }

    // ===== Saturation temperature at given pressure
    double EOSPengRobinson::saturationTemperature(double pressure) const
    {
        return m_impl->computeSaturationTemperature(pressure);
    }

}    // namespace PCProps::EquationOfState