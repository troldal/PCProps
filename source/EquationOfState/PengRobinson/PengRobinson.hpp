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
#include <tuple>
#include <string>

namespace PCProps::EquationOfState
{
    /**
     * @brief
     */
    class PengRobinson
    {
        using JSONString = std::string;

    public:
        // =====================================================================
        // CONSTRUCTORS, ASSIGNMENT & INITIATION
        // =====================================================================

        /**
         * @brief Default constructor. The object only supports copying and moving. The init
         * function must be called before using the other member functions.
         */
        PengRobinson();

        /**
         * @brief
         * @param pureComponent
         */
        PengRobinson(const std::function<double(std::string)>& constants, const std::function<double(std::string, double)>& correlations);

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

        /**
         * @brief Initiates an existing object with a new pure component.
         * @param pureComponent An object with an interface compatible with IPureComponent
         */
        void init(const std::function<double(std::string)>& constants, const std::function<double(std::string, double)>& correlations);

        // =====================================================================
        // FLASH ALGORITHMS
        // =====================================================================

        /**
         * @brief Compute flash at specified temperature and pressure.
         * @param temperature The temperature [K]
         * @param pressure The pressure [Pa]
         * @return The phase data for the phase resulting from the flash.
         */
        JSONString flashPT(double pressure, double temperature) const;

        /**
         * @brief Compute flash at specified pressure and computeEnthalpy.
         * @param pressure The pressure [Pa]
         * @param enthalpy The computeEnthalpy [J/mol]
         * @return The phase data for the phase(s) resulting from the flash.
         */
        JSONString flashPH(double pressure, double enthalpy) const;

        /**
         * @brief Compute flash at specified pressure and computeEntropy.
         * @param pressure The pressure [Pa]
         * @param entropy The computeEntropy [J/mol-K]
         * @return The phase data for the phase(s) resulting from the flash.
         */
        JSONString flashPS(double pressure, double entropy) const;

        /**
         * @brief Compute flash at specified pressure and vapor fraction.
         * @details If the given pressure is higher than the saturation pressure, a hypothetical saturation temperature
         * is estimated by extrapolation of the vapor pressure curve.
         * This point will be located in the supercritical region, where there are no saturation conditions.
         * However, in order to ensure numerical stability, calculating the (unphysical) saturation conditions are preferred.
         * @param pressure The pressure [Pa]
         * @param vaporFraction The vapor fraction [-]. Must be between 0.0 and 1.0.
         * @return The phase data for the phase(s) resulting from the flash.
         * @note If the vaporFraction = 0.0 or 1.0, or if the temperature is higher than the critical temperature, only one phase is returned.
         */
        JSONString flashPx(double pressure, double vaporFraction) const;

        /**
         * @brief Compute flash at specified temperature and vapor fraction.
         * @details If the given temperature is higher than the saturation temperature, a hypothetical saturation pressure
         * is estimated by extrapolation of the vapor pressure curve.
         * This point will be located in the supercritical region, where there are no saturation conditions.
         * However, in order to ensure numerical stability, calculating the (unphysical) saturation conditions are preferred.
         * @param temperature The temperature [K]
         * @param vaporFraction The vapor fraction [-]. Must be between 0.0 and 1.0.
         * @return The phase data for the phase(s) resulting from the flash.
         * @note If the vaporFraction = 0.0 or 1.0, or if the pressure is higher than the critical pressure, only one phase is returned.
         */
        JSONString flashTx(double temperature, double vaporFraction) const;

        /**
         * @brief
         * @param temperature
         * @param volume
         * @return
         */
        JSONString flashTV(double temperature, double volume) const;

        /**
         * @brief Calculate the saturation pressure at the given temperature.
         * @param temperature The temperature [K]
         * @return The saturation pressure [Pa]
         * @warning Returns NaN if the temperature is higher than the critical temperature.
         */
        double saturationPressure(double temperature) const;

        /**
         * @brief Calculate the saturation temperature at the given pressure.
         * @param pressure The pressure [Pa]
         * @return The saturation temperature [K]
         * @warning Returns NaN if pressure is higher than the critical pressure.
         */
        double saturationTemperature(double pressure) const;

    private:
        class impl;
        std::unique_ptr<impl> m_impl;
    };

}    // namespace PCProps::EquationOfState
#endif    // PCPROPS_PENGROBINSON_HPP
