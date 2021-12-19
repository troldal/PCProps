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

#include <cmath>
#include <limits>
#include <tuple>
#include <vector>

#include "PengRobinson.hpp"
#include <PhaseProperties.hpp>
#include <common/Globals.hpp>

#include <json/json.hpp>
#include <numerics.hpp>

using JSONString = std::string;

namespace PCProps::EquationOfState
{
    using PCProps::Globals::PI;
    using PCProps::Globals::R_CONST;
    using PCProps::Globals::STANDARD_P;
    using PCProps::Globals::STANDARD_T;
    using PCProps::Globals::MAX_ITER;
    using PCProps::Globals::TOLERANCE;

    class PengRobinson::impl
    {
    private:
        // ===== Basic fluid properties
        double m_criticalTemperature {};
        double m_criticalPressure {};
        double m_acentricFactor {};

        double m_normalFreezingPoint {};
        double m_normalBoilingPoint {};

        std::string m_name {};
        std::string m_CAS {};

        std::function<double(double)> m_idealGasCp {};
        std::function<double(double)> m_vaporPressure {};

        // ===== Calculated constants
        double m_ac {};
        double m_b {};
        double m_kappa {};

        mutable std::vector<PhaseProperties> m_phaseProps;

        /**
         * @brief Compute the 'a' coefficient for the Peng-robinson EOS.
         * @param temperature The temperature [K].
         * @return The 'a' coefficient.
         */
        inline double a(double temperature) const
        {
            auto result = m_ac * pow(1 + m_kappa * (1 - std::sqrt(temperature / m_criticalTemperature)), 2);
            if (std::isnan(result)) throw std::runtime_error("Numeric error: Coefficient 'a' could not be computed with T = " + std::to_string(temperature));
            return result;
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
            auto result = a(temperature) * pressure / pow(R_CONST, 2) / pow(temperature, 2);
            if (std::isnan(result))
                throw std::runtime_error("Numeric error: Coefficient 'A' could not be computed with T = " + std::to_string(temperature) + " and P = " + std::to_string(pressure));
            return result;
        }

        /**
         * @brief Compute the 'B' dimensionless coefficient for the Peng-Robinson EOS.
         * @param temperature The temperature [K].
         * @param pressure The pressure [Pa]
         * @return The 'B' coefficient.
         */
        inline double B(double temperature, double pressure) const
        {
            auto result = m_b * pressure / (R_CONST * temperature);
            if (std::isnan(result))
                throw std::runtime_error("Numeric error: Coefficient 'B' could not be computed with T = " + std::to_string(temperature) + " and P = " + std::to_string(pressure));
            return result;
        }

        /**
         * @brief Compute the 'alpha' coefficient for the Peng-Robinson EOS.
         * @param temperature The temperature [K].
         * @return The 'alpha' coefficient.
         */
        inline double alpha(double temperature) const
        {
            auto result = pow(1 + m_kappa * (1 - sqrt(temperature / m_criticalTemperature)), 2);
            if (std::isnan(result)) throw std::runtime_error("Numeric error: Coefficient 'alpha' could not be computed with T = " + std::to_string(temperature));
            return result;
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
            using std::cos;
            using std::sqrt;

            auto coeffA = A(temperature, pressure);
            auto coeffB = B(temperature, pressure);

            // ===== Compute the coefficients for solving Peng Robinson with respect to Z.
            auto a_0 = -(1.0 - coeffB);
            auto a_1 = (coeffA - 3.0 * pow(coeffB, 2.0) - 2.0 * coeffB);
            auto a_2 = -(coeffA * coeffB - pow(coeffB, 2.0) - pow(coeffB, 3.0));

            auto result = numeric::solve_cubic(a_0, a_1, a_2);
            result.erase(std::remove_if(result.begin(), result.end(), [](double root) { return root < 0.0; }), result.end());
            if (result.size() == 3) result.erase(result.begin() + 1);
            return result;
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

            auto result =
                exp(compressibility - 1 - log(compressibility - coeffB) -
                    coeffA / (coeffB * sqrt(8)) * log((compressibility + (1 + sqrt(2)) * coeffB) / (compressibility + (1 - sqrt(2)) * coeffB)));

            if (std::isnan(result))
                throw std::runtime_error(
                    "Numeric error: Fugacity coefficient could not be computed with T = " + std::to_string(temperature) + ", P = " + std::to_string(pressure) +
                    " and Z = " + std::to_string(compressibility));

            return result;
        }

        /**
         * @brief Compute the ideal gas enthalpy at the given T, relative the standard state.
         * @param temperature The temperature [K].
         * @return The ideal gas enthalpy [J/mol]
         */
        inline double idealGasEnthalpy(double temperature) const
        {
            using numeric::integrate;
            using PCProps::Globals::STANDARD_T;
            auto result = integrate([&](double t) { return m_idealGasCp(t); }, STANDARD_T, temperature);
            if (std::isnan(result)) throw std::runtime_error("Numeric error: Ideal gas enthalpy could not be computed with T = " + std::to_string(temperature));
            return result;
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
            using std::sqrt;

            auto coeffA = A(temperature, pressure);
            auto coeffB = B(temperature, pressure);

            auto result = (compressibility - 1.0 -
                           coeffA / (coeffB * sqrt(8)) * (1 + m_kappa * sqrt(temperature / m_criticalTemperature) / sqrt(alpha(temperature))) *
                               log((compressibility + (1 + sqrt(2)) * coeffB) / (compressibility + (1 - sqrt(2)) * coeffB))) *
                          R_CONST * temperature;

            if (std::isnan(result))
                throw std::runtime_error(
                    "Numeric error: Enthalpy departure could not be computed with T = " + std::to_string(temperature) + ", P = " + std::to_string(pressure) +
                    " and Z = " + std::to_string(compressibility));

            return result;
        }

        /**
         * @brief Compute the ideal gas entropy at the given T and P, relative the standard state.
         * @param temperature The temperature [K].
         * @param pressure The pressure [Pa].
         * @return The ideal gas entropy [J/mol-K]
         */
        inline double idealGasEntropy(double temperature, double pressure) const
        {
            using numeric::integrate;
            auto result = integrate([&](double temp) { return m_idealGasCp(temp) / temp; }, PCProps::Globals::STANDARD_T, temperature) - R_CONST * log(pressure / STANDARD_P);

            if (std::isnan(result))
                throw std::runtime_error("Numeric error: Ideal gas entropy could not be computed with T = " + std::to_string(temperature) + " and P = " + std::to_string(pressure));
            return result;
        }

        /**
         * @brief Compute the entropy departure relative to standard (ideal gas) state, at the given T and P.
         * @param temperature The temperature [K].
         * @param pressure The pressure [Pa]
         * @param compressibility The compressibility factor for the fluid [-]
         * @return The entropy departure (S-S^ig) for the fluid [J/mol-K]
         */
        inline double entropyDeparture(double temperature, double pressure, double compressibility) const
        {
            using std::log;
            using std::sqrt;

            auto coeffA = A(temperature, pressure);
            auto coeffB = B(temperature, pressure);

            auto result = (log(compressibility - coeffB) - coeffA / (coeffB * sqrt(8)) * m_kappa * sqrt(temperature / m_criticalTemperature) / sqrt(alpha(temperature)) *
                                                               log((compressibility + (1 + sqrt(2)) * coeffB) / (compressibility + (1 - sqrt(2)) * coeffB))) * R_CONST;

            if (std::isnan(result))
                throw std::runtime_error(
                    "Numeric error: Enthalpy departure could not be computed with T = " + std::to_string(temperature) + ", P = " + std::to_string(pressure) +
                    " and Z = " + std::to_string(compressibility));

            return result;
        }

        /**
         * @brief
         * @param temperature
         * @param molarVolume
         * @return
         */
        inline double calcPressure(double temperature, double molarVolume) const {
            auto result = R_CONST * temperature / (molarVolume - m_b) - (m_ac * alpha(temperature)) / (pow(molarVolume, 2) + 2 * m_b * molarVolume - pow(m_b, 2));

            if (std::isnan(result))
                throw std::runtime_error(
                    "Numeric error: Pressure could not be computed with T = " + std::to_string(temperature) + " and V = " + std::to_string(molarVolume));

            return result;
        }

    public:
        /**
         * @brief Constructor.
         * @param criticalTemperature The critical temperature [K].
         * @param criticalPressure The critical pressure [Pa]
         * @param acentricFactor The acentric factor [-]
         */
        explicit impl(const std::function<double(std::string)>& constants, const std::function<double(std::string, double)>& correlations)
            : m_criticalTemperature(constants("CriticalTemperature")),
              m_criticalPressure(constants("CriticalPressure")),
              m_acentricFactor(constants("AcentricFactor")),
              m_normalFreezingPoint(constants("NormalFreezingPoint")),
              m_normalBoilingPoint(constants("NormalBoilingPoint")),
              m_idealGasCp([=](double t) -> double { return correlations("IdealGasCp", t); }),
              m_vaporPressure([=](double t) -> double { return correlations("VaporPressure", t); }),
              m_ac(0.45723553 * pow(PCProps::Globals::R_CONST, 2) * pow(m_criticalTemperature, 2) / m_criticalPressure),
              m_b(0.07779607 * PCProps::Globals::R_CONST * m_criticalTemperature / m_criticalPressure),
              m_kappa(
                  m_acentricFactor <= 0.49 ? 0.37464 + 1.54226 * m_acentricFactor - 0.26992 * pow(m_acentricFactor, 2)
                                           : 0.379642 + 1.48503 * m_acentricFactor - 0.164423 * pow(m_acentricFactor, 2) + 0.016666 * pow(m_acentricFactor, 3))
        {}

        /**
         * @brief Compute the compressibilities and fugacity coefficients for the phases present at the given T and P.
         * @param temperature The temperature [K]
         * @param pressure The pressure [Pa].
         * @return A vector of std::pairs with the compressibility and fugacity coefficient for each phase.
         * @note A flash may yield Z and phi for two phases, even though only one phase is present. The individual
         * flash algorithm will discard the Z/Phi pair that is invalid.
         */
        inline std::vector<std::pair<double, double>> computeCompressibilityAndFugacity(double temperature, double pressure) const
        {
            // ===== Compute the compressibility factor(s)
            auto zs = computeCompressibilityFactors(temperature, pressure);

            // ===== Compute the fugacity coefficient for each Z and add the pair(s) to the results vector.
            std::vector<std::pair<double, double>> result;
            result.reserve(zs.size());
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
            if (temperature == m_criticalTemperature) return m_criticalPressure;

            using std::get;
            using std::min;

            // Define lambda function, for guessing the saturation pressure.
            // First, the saturation pressure is estimated using a correlation. If
            // the estimated pressure in the two-phase region of the van der Waals loop,
            // this pressure will be returned.
            // If not, a range of +/- 2% of the estimate will be defined. This range is sliced into
            // smaller and smaller segments, until a pressure in the two-phase region is found.
            // When found, the result will be returned.
            auto guessSaturationPressure = [&]() {
                // Guess vapor pressure using correlation.
                auto guess = m_vaporPressure(temperature);

                // If the guess leads to two compressibility factors, the guess is in the van der Waals loop,
                // and can be returned as the guess
                if (computeCompressibilityFactors(temperature, guess).size() == 2) return std::max(1.0, guess);

                // The findRange function takes the upper and lower bounds of a range,
                // and narrows in on the sub-range where the van der Waals loop, and
                // therefore the vapor pressure, can be found.
                // The input range must contain the van der Waals loop.
                auto findRange = [&](double lower, double upper) {
                    // ===== First, slice up the range in intervals
                    auto                sliceCount = 20;
                    std::vector<double> slices {};
                    for (int i = 0; i <= sliceCount; ++i) slices.emplace_back(lower + i * ((upper - lower) / sliceCount));

                    // ===== For each interval, compute the slope, dz/dp
                    std::vector<std::pair<double, double>> slopes;
                    for (auto iter = std::next(slices.begin()); iter != slices.end(); ++iter) {
                        auto z1 = computeCompressibilityFactors(temperature, *std::prev(iter));
                        auto z2 = computeCompressibilityFactors(temperature, *iter);
                        slopes.emplace_back(
                            std::make_pair((*std::prev(iter) + *iter) / 2, ((z1.front() + z1.back()) / 2 - (z2.front() + z2.back()) / 2) / (*std::prev(iter) - *iter)));
                    }

                    // ===== Select the four intervals with the steepest slope, and return the range.
                    std::sort(slopes.begin(), slopes.end(), [&](const auto& a, const auto& b) { return a.second < b.second; });
                    std::vector<double> result { slopes[0].first, slopes[1].first, slopes[2].first, slopes[3].first, slopes[4].first };
                    std::sort(result.begin(), result.end());
                    return std::vector<double> { result.front(), result.back() };
                };

                // ===== Assume that the true vapor pressure is within +/- 10% of the initial guess.
                std::vector<double> guesses { guess * 0.9, guess * 1.1 };

                // ===== Continually narrow the range until the van der Waals loop has been found,
                // ===== or the maximum number of iterations has been exceeded.
                int counter { 0 };
                while (true) {
                    if (counter >= MAX_ITER || computeCompressibilityFactors(temperature, (guesses.front() + guesses.back()) / 2.0).size() == 2)
                        return (guesses.front() + guesses.back()) / 2.0;

                    ++counter;
                    guesses = findRange(guesses.front(), guesses.back());
                }
            };

            // ===== If the temperature is less than the critical temperature, compute using the PR-EOS.
            if (temperature < m_criticalTemperature) {
                auto guess = guessSaturationPressure();
                int  counter { 0 };
                while (true) {
                    auto phi   = computeCompressibilityAndFugacity(temperature, guess);
                    auto phi_l = get<1>(phi.front());
                    auto phi_v = get<1>(phi.back());

                    // ===== Iterate using successive substitution. The minimum vapor pressure is 1.0 Pa, as lower
                    // ===== values may lead to numerical issues.
                    if (abs(phi_l / phi_v - 1) < TOLERANCE || counter >= MAX_ITER || guess < 1.0) return std::max(1.0, guess);
                    ++counter;
                    guess = guess * (phi_l / phi_v);
                }
            }

            // ===== Otherwise, compute the slope at the critical point, and extrapolate.
            auto slope = numeric::diff_backward([&](double t) { return m_vaporPressure(t); }, m_criticalTemperature);
            return m_criticalPressure + (temperature - m_criticalTemperature) * slope;
        }

        /**
         * @brief Compute the fluid saturation temperature at the given P.
         * @param pressure The pressure [Pa].
         * @return The saturation temperature [K].
         */
        inline double computeSaturationTemperature(double pressure) const
        {
            if (pressure == m_criticalPressure) return m_criticalTemperature;

            // ===== If the pressure is less than the critical pressure, compute using the PR-EOS.
            if (pressure < m_criticalPressure) {
                auto result = numeric::newton([&](double t) { return computeSaturationPressure(t) - pressure; }, m_criticalTemperature * 0.5);
                if (std::isnan(result)) throw std::runtime_error("Saturation temperature could not be calculated.");
                return result;
            }

            // ===== Otherwise, compute the slope at the critical point, and extrapolate.
            auto slope = numeric::diff_backward([&](double t) { return m_vaporPressure(t); }, m_criticalTemperature);
            return (pressure - m_criticalPressure) / slope + m_criticalTemperature;
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

        /**
         * @brief
         * @param temperature
         * @param pressure
         * @return
         */
        inline void computeThermodynamicProperties(double temperature, double pressure) const
        {
            // First, some housekeeping
            using PCProps::Globals::R_CONST;
            using std::get;
            using std::sqrt;
            auto eps = sqrt(std::numeric_limits<double>::epsilon());
            m_phaseProps.clear();

            // Calculate compressibility and fugacity coefficient for all phases at the given T and P.
            // For all phases, calculate basic thermodynamic properties.
            auto z_phi = computeCompressibilityAndFugacity(temperature, pressure);
            for (decltype(z_phi.size()) index = 0; index < z_phi.size(); ++index) {
                PhaseProperties data;

                auto determinePhaseType = [&]() {
                    if (z_phi.size() == 2 && index == 0) return PhaseType::Liquid;
                    if (z_phi.size() == 2 && index == 1) return PhaseType::Vapor;
                    if (z_phi.size() == 1 && pressure >= computeSaturationPressure(temperature)) return PhaseType::Liquid;
                    if (z_phi.size() == 1 && pressure < computeSaturationPressure(temperature)) return PhaseType::Vapor;
                    return PhaseType::Undefined;
                };

                data.Type                = determinePhaseType();
                data.Pressure            = pressure;
                data.Temperature         = temperature;
                data.Compressibility     = get<0>(z_phi[index]);
                data.FugacityCoefficient = get<1>(z_phi[index]);
                data.Enthalpy            = computeEnthalpy(temperature, pressure, data.Compressibility);
                data.Entropy             = computeEntropy(temperature, pressure, data.Compressibility);
                data.MolarVolume         = data.Compressibility * R_CONST * temperature / pressure;
                data.GibbsEnergy         = data.Enthalpy - temperature * data.Entropy;
                data.InternalEnergy      = data.Enthalpy - pressure * data.MolarVolume;
                data.HelmholzEnergy      = data.InternalEnergy - temperature * data.Entropy;
                data.VaporPressure       = computeSaturationPressure(temperature);
                data.CriticalPressure    = m_criticalPressure;
                data.CriticalTemperature = m_criticalTemperature;
                data.NormalFreezingPoint = m_normalFreezingPoint;
                data.NormalBoilingPoint  = m_normalBoilingPoint;

                auto dvdt = numeric::diff_central([&](double t) { return computeCompressibilityFactors(t, pressure)[index] * R_CONST * t / pressure; }, temperature);
                auto dvdp = numeric::diff_central([&](double p) { return computeCompressibilityFactors(temperature, p)[index] * R_CONST * temperature / p; }, pressure);
                auto dpdt = numeric::diff_central([&](double t) { return calcPressure(t, data.MolarVolume); }, temperature);
                auto dpdv = numeric::diff_central([&](double v) { return calcPressure(temperature, v); }, data.MolarVolume);

                data.Cp = numeric::diff_central([&](double t) { return computeEnthalpy(t, pressure, computeCompressibilityFactors(t, pressure)[index]); }, temperature);
                data.Cv = data.Cp + temperature * pow(dpdt, 2) / dpdv;
                data.ThermalExpansionCoefficient = (1.0 / data.MolarVolume) * dvdt;
                data.JouleThomsonCoefficient     = (temperature * dvdt - data.MolarVolume) / data.Cp;
                data.IsothermalCompressibility   = -(1.0 / data.MolarVolume) * dvdp;
                // TODO: This calculation does not seem to yield correct results!
                data.SpeedOfSound = sqrt(pow(data.MolarVolume, 2) * ((temperature / data.Cv) * pow(dpdt, 2) - dpdv));

                m_phaseProps.emplace_back(data);
            }
        }

        /**
         * @brief
         * @param pressure
         * @param temperature
         * @return
         */
        inline std::vector<PhaseProperties> flashPT(double pressure, double temperature) const
        {
            // ===== Compute compressibility factors and fugacity coefficients at given T and P.
            computeThermodynamicProperties(temperature, pressure);

            // ===== Set the molar flow to 1.0, as only one phase will result from a PT flash.
            for (auto& phase : m_phaseProps) phase.MolarFlow = 1.0;

            // ===== Return the phase with the lowest fugacity.
            return { *std::min_element(m_phaseProps.begin(), m_phaseProps.end(), [](const PhaseProperties& a, const PhaseProperties& b) {
                return a.FugacityCoefficient < b.FugacityCoefficient;
            }) };
        }

        /**
         * @brief
         * @param temperature
         * @param vaporFraction
         * @return
         */
        inline std::vector<PhaseProperties> flashTx(double temperature, double vaporFraction) const
        {
            // ===== If the temperature < Tc
            if (temperature < m_criticalTemperature) {
                // ===== First, calculate the saturation pressure at the specified pressure.
                auto pressure = computeSaturationPressure(temperature);
                computeThermodynamicProperties(temperature, pressure);

                // ===== If the specified vapor fraction is 1.0 (or higher), the fluid is a saturated vapor.
                if (vaporFraction >= 1.0) {
                    m_phaseProps.back().MolarFlow = 1.0;
                    return { m_phaseProps.back() };
                }

                // ===== If the specified vapor fraction is 0.0 (or lower), the fluid is a saturated liquid.
                if (vaporFraction <= 0.0) {
                    m_phaseProps.front().MolarFlow = 1.0;
                    return { m_phaseProps.front() };
                }

                // ===== If the vapor fraction is between 0.0 and 1.0, the fluid is two-phase.
                m_phaseProps[0].MolarFlow = 1 - vaporFraction;
                m_phaseProps[1].MolarFlow = vaporFraction;
                return { m_phaseProps[0], m_phaseProps[1] };
            }

            // ===== If the temperature >= Tc, calculate the hypothetical saturation conditions in the supercritical region.
            auto pressure = computeSaturationPressure(temperature);
            return flashPT(pressure, temperature);
        }

        /**
         * @brief
         * @param pressure
         * @param vaporFraction
         * @return
         */
        inline std::vector<PhaseProperties> flashPx(double pressure, double vaporFraction) const
        {
            auto temperature = computeSaturationTemperature(pressure);
            return flashTx(temperature, vaporFraction);
        }

        /**
         * @brief
         * @param pressure
         * @param enthalpy
         * @return
         */
        inline std::vector<PhaseProperties> flashPH(double pressure, double enthalpy) const
        {
            using std::get;

            // ===== Define objective function
            auto enthalpyObjFunction = [&](double t) {
                auto z_phi = computeCompressibilityAndFugacity(t, pressure);
                auto z     = get<0>(*std::min_element(z_phi.begin(), z_phi.end(), [](const auto& a, const auto& b) { return get<1>(a) < get<1>(b); }));
                return computeEnthalpy(t, pressure, z) - enthalpy;
            };

            // ===== If the fluid is supercritical, compute like so... (single phase)
            if (pressure >= m_criticalPressure) {
                // ===== Use the (hypothetical) saturation temperature as first guess, and compute until the enthalpy is found.
                return flashPT(pressure, numeric::newton(enthalpyObjFunction, computeSaturationTemperature(pressure)));
            }

            // ===== Otherwise, the fluid is sub-critical...
            // ===== First, calculate the saturation properties at the specified pressure.
            auto temperature  = computeSaturationTemperature(pressure);
            auto z_phi        = computeCompressibilityAndFugacity(temperature, pressure);
            auto [z_v, phi_v] = z_phi.back();
            auto [z_l, phi_l] = z_phi.front();
            auto h_v          = computeEnthalpy(temperature, pressure, z_v);
            auto h_l          = computeEnthalpy(temperature, pressure, z_l);

            // ===== If the specified enthalpy is lower than the saturated liquid enthalpy, the fluid is a compressed liquid.
            if (enthalpy < h_l) {
                return flashPT(pressure, numeric::newton(enthalpyObjFunction, temperature * 0.8));
            }

            // ===== If the specified enthalpy is higher than the saturated vapor entropy, the fluid is superheated vapor.
            if (enthalpy > h_v) {
                return flashPT(pressure, numeric::newton(enthalpyObjFunction, temperature * 1.2));
            }

            // ===== If the fluid is not a compressed liquid nor a superheated vapor, the fluid is two-phase.
            auto vaporFraction = (h_l - enthalpy) / (h_l - h_v);
            return flashPx(pressure, vaporFraction);
        }

        /**
         * @brief
         * @param pressure
         * @param entropy
         * @return
         */
        inline std::vector<PhaseProperties> flashPS(double pressure, double entropy) const
        {
            using std::get;

            // ===== Define objective function
            auto entropyObjFunction = [&](double t) {
                auto z_phi = computeCompressibilityAndFugacity(t, pressure);
                auto z     = get<0>(*std::min_element(z_phi.begin(), z_phi.end(), [](const auto& a, const auto& b) { return get<1>(a) < get<1>(b); }));
                return computeEntropy(t, pressure, z) - entropy;
            };

            if (pressure >= m_criticalPressure) {
                return flashPT(pressure, numeric::newton(entropyObjFunction, computeSaturationTemperature(pressure)));
            }

            // ===== First, calculate the saturation properties at the specified pressure.
            auto temperature  = computeSaturationTemperature(pressure);
            auto z_phi        = computeCompressibilityAndFugacity(temperature, pressure);
            auto [z_v, phi_v] = z_phi.back();
            auto [z_l, phi_l] = z_phi.front();
            auto s_v          = computeEntropy(temperature, pressure, z_v);
            auto s_l          = computeEntropy(temperature, pressure, z_l);

            // ===== If the specified entropy is lower than the saturated liquid entropy, the fluid is a compressed liquid.
            if (entropy < s_l) {
                return flashPT(pressure, numeric::newton(entropyObjFunction, temperature * 0.8));
            }

            // ===== If the specified entropy is higher than the saturated vapor entropy, the fluid is superheated vapor.
            if (entropy > s_v) {
                return flashPT(pressure, numeric::newton(entropyObjFunction, temperature * 1.2));
            }

            // ===== If the fluid is not a compressed liquid nor a superheated vapor, the fluid is two-phase.
            auto vaporFraction = (s_l - entropy) / (s_l - s_v);
            return flashPx(pressure, vaporFraction);
        }

        /**
         * @brief
         * @param temperature
         * @param volume
         * @return
         */
        inline std::vector<PhaseProperties> flashTV(double temperature, double volume) const
        {
            using std::pow;

            // ===== Fluid is supercritical
            if (temperature >= m_criticalTemperature) {
                return flashPT(calcPressure(temperature, volume), temperature);
            }

            flashTx(temperature, 0.5);

            // ===== Fluid is a compressed liquid or super-heated vapor
            if (volume < m_phaseProps.front().MolarVolume || volume > m_phaseProps.back().MolarVolume) {
                return flashPT(calcPressure(temperature, volume), temperature);
            }

            // ===== Fluid is multiphase
            auto vaporFraction = (volume - m_phaseProps.front().MolarVolume) / (m_phaseProps.back().MolarVolume - m_phaseProps.front().MolarVolume);
            return flashTx(temperature, vaporFraction);
        }

    };    // PengRobinson::impl

    // ===== Constructor, default
    PengRobinson::PengRobinson() : m_impl(nullptr) {};

    // ===== Constructor
    PengRobinson::PengRobinson(const std::function<double(std::string)>& constants, const std::function<double(std::string, double)>& correlations)
        : m_impl(std::make_unique<impl>(constants, correlations))
    {}

    // ===== Copy constructor
    PengRobinson::PengRobinson(const PengRobinson& other) : m_impl(other.m_impl ? std::make_unique<impl>(*other.m_impl) : nullptr) {};

    // ===== Move constructor
    PengRobinson::PengRobinson(PengRobinson&& other) noexcept = default;

    // ===== Destructor
    PengRobinson::~PengRobinson() = default;

    // ===== Copy assignment operator
    PengRobinson& PengRobinson::operator=(const PengRobinson& other)
    {
        PengRobinson copy = other;
        *this             = std::move(copy);
        return *this;
    };

    // ===== Move assignment operator
    PengRobinson& PengRobinson::operator=(PengRobinson&& other) noexcept = default;

    // ===== Initiates object with new pure component data
    void PengRobinson::init(const std::function<double(std::string)>& constants, const std::function<double(std::string, double)>& correlations)
    {
        m_impl = std::make_unique<impl>(constants, correlations);
    }

    // ===== Flash calculations at different specifications
    JSONString PengRobinson::flash(const std::string& specification, double var1, double var2) const
    {
        std::vector<nlohmann::json> result;

        if (specification == "PT")
            for (const auto& phase : m_impl->flashPT(var1, var2)) result.emplace_back(nlohmann::json::parse(phase.asJSON()));

        if (specification == "Tx")
            for (const auto& phase : m_impl->flashTx(var1, var2)) result.emplace_back(nlohmann::json::parse(phase.asJSON()));

        if (specification == "Px")
            for (const auto& phase : m_impl->flashPx(var1, var2)) result.emplace_back(nlohmann::json::parse(phase.asJSON()));

        if (specification == "PH")
            for (const auto& phase : m_impl->flashPH(var1, var2)) result.emplace_back(nlohmann::json::parse(phase.asJSON()));

        if (specification == "PS")
            for (const auto& phase : m_impl->flashPS(var1, var2)) result.emplace_back(nlohmann::json::parse(phase.asJSON()));

        if (specification == "TV")
            for (const auto& phase : m_impl->flashTV(var1, var2)) result.emplace_back(nlohmann::json::parse(phase.asJSON()));

        for (auto& phase : result) {
        }

        return nlohmann::json(result).dump();
    }

    // ===== Saturation pressure at given temperature
    double PengRobinson::saturationPressure(double temperature) const
    {
        return m_impl->computeSaturationPressure(temperature);
    }

    // ===== Saturation temperature at given pressure
    double PengRobinson::saturationTemperature(double pressure) const
    {
        return m_impl->computeSaturationTemperature(pressure);
    }

}    // namespace PCProps::EquationOfState