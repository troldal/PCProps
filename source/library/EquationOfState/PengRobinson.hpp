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

#ifndef PCPROPS_PENGROBINSON_HPP
#define PCPROPS_PENGROBINSON_HPP

#include <functional>
#include <memory>

#include <library/PCPropsData.hpp>

namespace PCProps::EquationOfState
{
    /**
     * @brief
     */
    class PengRobinson
    {
    public:
        // =====================================================================
        // CONSTRUCTORS & ASSIGNMENT OPERATORS
        // =====================================================================

        /**
         * @brief Default constructor. All member variables set to default values.
         */
        PengRobinson();

        /**
         * @brief
         * @param criticalTemperature
         * @param criticalPressure
         * @param acentricFactor
         */
        PengRobinson(double criticalTemperature, double criticalPressure, double acentricFactor);

        /**
         * @brief Copy constructor
         */
        PengRobinson(const PengRobinson& other);

        /**
         * @brief Move constructor
         */
        PengRobinson(PengRobinson&& other) noexcept;

        /**
         * @brief Destructor
         */
        ~PengRobinson();

        /**
         * @brief Copy assignment operator
         */
        PengRobinson& operator=(const PengRobinson& other);

        /**
         * @brief Move assignment operator
         */
        PengRobinson& operator=(PengRobinson&& other) noexcept;

        // =====================================================================
        // MANIPULATORS
        // =====================================================================

        /**
         * @brief
         * @param criticalTemperature
         * @param criticalPressure
         * @param acentricFactor
         */
        void setProperties(double criticalTemperature, double criticalPressure, double acentricFactor);

        /**
         * @brief
         * @param vaporPressureFunction
         */
        void setVaporPressureFunction(const std::function<double(double)>& vaporPressureFunction);

        /**
         * @brief
         * @param idealGasCpFunction
         */
        void setIdealGasCpFunction(const std::function<double(double)>& idealGasCpFunction);

        /**
         * @brief
         * @param idealGasCpDerivativeFunction
         */
        void setIdealGasCpDerivativeFunction(const std::function<double(double)>& idealGasCpDerivativeFunction);

        /**
         * @brief
         * @param idealGasCpIntegralFunction
         */
        void setIdealGasCpIntegralFunction(const std::function<double(double)>& idealGasCpIntegralFunction);

        /**
         * @brief
         * @param idealGasOverTIntegralFunction
         */
        void setIdealGasCpOverTIntegralFunction(const std::function<double(double)>& idealGasOverTIntegralFunction);

        // =====================================================================
        // FLASH ALGORITHMS
        // =====================================================================

        /**
         * @brief Compute flash at specified temperature and pressure.
         * @param temperature The temperature [K]
         * @param pressure The pressure [Pa]
         * @return The phase data for the phase resulting from the flash.
         */
        PCPhases flashPT(double pressure, double temperature) const;

        /**
         * @brief Compute flash at specified temperature and vapor fraction.
         * @param temperature The temperature [K]
         * @param vaporFraction The vapor fraction [-]. Must be between 0.0 and 1.0.
         * @return The phase data for the phase(s) resulting from the flash.
         */
        PCPhases flashTx(double temperature, double vaporFraction) const;

        /**
         * @brief Compute flash at specified pressure and vapor fraction
         * @param pressure The pressure [Pa]
         * @param vaporFraction The vapor fraction [-]. Must be between 0.0 and 1.0.
         * @return The phase data for the phase(s) resulting from the flash.
         */
        PCPhases flashPx(double pressure, double vaporFraction) const;

        /**
         * @brief Compute flash at specified pressure and computeEnthalpy.
         * @param pressure The pressure [Pa]
         * @param enthalpy The computeEnthalpy [J/mol]
         * @return The phase data for the phase(s) resulting from the flash.
         */
        PCPhases flashPH(double pressure, double enthalpy) const;

        /**
         * @brief Compute flash at specified pressure and computeEntropy.
         * @param pressure The pressure [Pa]
         * @param entropy The computeEntropy [J/mol-K]
         * @return The phase data for the phase(s) resulting from the flash.
         */
        PCPhases flashPS(double pressure, double entropy) const;

        /**
         * @brief
         * @param temperature
         * @param volume
         * @return
         */
        PCPhases flashTV(double temperature, double volume) const;

        /**
         * @brief
         * @param temperature
         * @return
         */
        double saturationPressure(double temperature) const;

        /**
         * @brief
         * @param pressure
         * @return
         */
        double saturationTemperature(double pressure) const;

    private:
        class impl;
        std::unique_ptr<impl> m_impl;
    };

}    // namespace PCProps::EquationOfState
#endif    // PCPROPS_PENGROBINSON_HPP
