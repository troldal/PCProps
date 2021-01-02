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

#ifndef PCPROPS_EOSPENGROBINSON_HPP
#define PCPROPS_EOSPENGROBINSON_HPP

#include <functional>
#include <tuple>
#include <utility>

#include "EOSUtilities.hpp"

namespace PCProps::EquationOfState
{
    /**
     * @brief
     */
    class EOSPengRobinson
    {
        // ===== Basic fluid properties
        double m_criticalTemperature {};
        double m_criticalPressure {};
        double m_molecularWeight {};

        // ===== Calculated constants
        double m_kappa {};
        double m_b {};

        // ===== Temperature and/or pressure dependent coefficients
        std::function<double(double)>         m_alpha {};
        std::function<double(double)>         m_a {};
        std::function<double(double, double)> m_A {};
        std::function<double(double, double)> m_B {};

        // ===== User-supplied correlations
        std::function<double(double)>         m_vaporPressureFunction {};
        std::function<double(double)>         m_idealGasCpFunction {};
        std::function<double(double, double)> m_idealGasCpIntegralFunction {};
        std::function<double(double, double)> m_idealGasCpOverTemperatureIntegralFunction {};

        // ===== Private member functions
        PhaseData createEOSData(double moles, double temperature, double pressure, double z, double phi) const;

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
        EOSPengRobinson(
            double                               criticalTemperature,
            double                               criticalPressure,
            double                               acentricFactor,
            double                               molecularWeight,
            const std::function<double(double)>& vaporPressureFunction,
            const std::function<double(double)>& idealGasCpFunction);

        /**
         * @brief Copy constructor
         */
        EOSPengRobinson(const EOSPengRobinson& other);

        /**
         * @brief Move constructor
         */
        EOSPengRobinson(EOSPengRobinson&& other) noexcept;

        /**
         * @brief Destructor
         */
        ~EOSPengRobinson();

        /**
         * @brief Copy assignment operator
         */
        EOSPengRobinson& operator=(const EOSPengRobinson& other);

        /**
         * @brief Move assignment operator
         */
        EOSPengRobinson& operator=(EOSPengRobinson&& other) noexcept;

        /**
         * @brief Compute flash at specified temperature and pressure.
         * @param temperature The temperature [K]
         * @param pressure The pressure [Pa]
         * @param moles The number of moles in the fluid. Default is 1.0
         * @return The phase data for the phase resulting from the flash.
         */
        Phases flashPT(double pressure, double temperature, double moles = 1.0) const;

        /**
         * @brief Compute flash at specified temperature and vapor fraction.
         * @param temperature The temperature [K]
         * @param vaporFraction The vapor fraction [-]. Must be between 0.0 and 1.0.
         * @param moles The number of moles in the fluid. Default is 1.0
         * @return The phase data for the phase(s) resulting from the flash.
         */
        Phases flashTx(double temperature, double vaporFraction, double moles = 1.0) const;

        /**
         * @brief Compute flash at specified pressure and vapor fraction
         * @param pressure The pressure [Pa]
         * @param vaporFraction The vapor fraction [-]. Must be between 0.0 and 1.0.
         * @param moles The number of moles in the fluid. Default is 1.0
         * @return The phase data for the phase(s) resulting from the flash.
         */
        Phases flashPx(double pressure, double vaporFraction, double moles = 1.0) const;

        /**
         * @brief Compute flash at specified pressure and enthalpy.
         * @param pressure The pressure [Pa]
         * @param enthalpy The enthalpy [J/mol]
         * @param moles The number of moles in the fluid. Default is 1.0
         * @return The phase data for the phase(s) resulting from the flash.
         */
        Phases flashPH(double pressure, double enthalpy, double moles = 1.0) const;

        /**
         * @brief Compute flash at specified pressure and entropy.
         * @param pressure The pressure [Pa]
         * @param entropy The entropy [J/mol-K]
         * @param moles The number of moles in the fluid. Default is 1.0
         * @return The phase data for the phase(s) resulting from the flash.
         */
        Phases flashPS(double pressure, double entropy, double moles = 1.0) const;
    };

}    // namespace PCProps::EquationOfState
#endif    // PCPROPS_EOSPENGROBINSON_HPP
