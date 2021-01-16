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

#ifndef PCPROPS_PCCOMPONENT_HPP
#define PCPROPS_PCCOMPONENT_HPP

#include <functional>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <type_traits>

#include <library/PCEquationOfState.hpp>
#include <library/PCPropsData.hpp>
#include <library/PCPropsException.hpp>
#include <library/Utilities/NamedTypes.hpp>

namespace PCProps
{
    /**
     * @brief The PCComponentData struct holds all the raw data and function objects that define a pure component.
     * An object of this type can be used to construct a PCComponent object.
     * @details A PCComponentData object can include all the data necessary to define a pure component, and as such
     * can be used directly. However, the intent is to use the object to construct a PCComponent object. Please see
     * documentation for PCComponent for further information.
     */
    struct PCComponentData
    {
        std::string name;    /**< The name of the component, e.g. Methane */
        std::string formula; /**< The formula of the component, e.g. CH4 */
        std::string casrn;   /** The CAS registration number for the component, e.g. 75-07-0 */
        std::string smiles;  /**< The SMILES string for the component, e.g. C */

        std::optional<double> molecularWeight {};         /**< The molecular weight [g/mol] of a component, e.g. 16.043 g/mol */
        std::optional<double> boilingTemperature {};      /**< The normal boiling temperature [K] (at 101325 Pa) of a component */
        std::optional<double> freezingTemperature {};     /**< The freezing point temperature [K] of a component */
        std::optional<double> criticalTemperature {};     /**< The critical temperature [K] of a component */
        std::optional<double> criticalPressure {};        /**< The critical pressure [Pa] of a component */
        std::optional<double> criticalVolume {};          /**< The critical molar volume [m3/mol] of a component */
        std::optional<double> criticalDensity {};         /**< The critical density [kg/m3] of a component */
        std::optional<double> criticalCompressibility {}; /**< The critical compressibility factor [-] of a component */
        std::optional<double> acentricFactor {};          /**< The acentric factor (omega) of a component */
        std::optional<double> dipoleMoment {};            /**< The dipole moment of the component */

        PCProps::PCEquationOfState            equationOfState {};
        std::function<double(double)>         idealGasCpCorrelation {};
        std::function<double(double)>         liquidCpCorrelation {};                  /**< The liquid Cp [HOLD] as a function of temperature [K] */
        std::function<double(double)>         vaporPressureCorrelation {};             /**< The vapor pressure [Pa] as a function of temperature [K] */
        std::function<double(double)>         surfaceTensionCorrelation {};            /**< The surface tension [HOLD] as a function of temperature [K] */
        std::function<double(double)>         heatOfVaporizationCorrelation {};        /**< The latent heat [HOLD] as a function of temperature [K] */
        std::function<double(double)>         vaporThermalConductivityCorrelation {};  /**< The vapor thermal conductivity [HOLD] as a function of temperature [K] */
        std::function<double(double)>         liquidThermalConductivityCorrelation {}; /**< The liquid thermal conductivity [HOLD] as a function of temperature [K] */
        std::function<double(double)>         saturatedVaporViscosityCorrelation {};   /**< The vapor viscosity [HOLD] as a function of temperature [K] */
        std::function<double(double)>         saturatedLiquidViscosityCorrelation {};  /**< The liquid viscosity [HOLD] as a function of temperature [K] */
        std::function<double(double)>         saturatedLiquidVolumeCorrelation {};     /**< The liquid density [HOLD] as a function of temperature [K] */
        std::function<double(double, double)> compressedLiquidVolumeCorrelation {};
    };

    /**
     * @brief A PCComponent object holds an PCComponentData object, as well as all required access functions.
     * @details In order for a PCComponent object to work properly, it should be constructed from a PCComponentData object.
     * The reason for dividing the definition of a pure component into two different classes is that the PCComponent does not
     * have any setter functions, an therefore cannot be modified after construction. The purpose is to avoid accidental
     * modification of pure component data after construction. If it is absolutely necessary to modify an existing PCComponent
     * object after construction, this should be done by copying the data object, modify the data and construct a new object
     * in the same place as the old object.
     */
    class PCComponent
    {
        // ===== Private Data Members
        PCComponentData m_data = {}; /**< The data object, holding raw data and correlations for a pure component */

        // ===== Public members
    public:
        /**
         * @brief Default constructor.
         */
        PCComponent();

        /**
         * @brief Constructor, taking a PCComponentData object as an argument.
         * @param data A PCComponentData object with the component data.
         */
        explicit PCComponent(const PCComponentData& data);

        /**
         * @brief Copy constructor.
         */
        PCComponent(const PCComponent& other);

        /**
         * @brief Move constructor.
         */
        PCComponent(PCComponent&& other) noexcept;

        /**
         * @brief Destructor.
         */
        ~PCComponent();

        /**
         * @brief Copy assignment operator.
         */
        PCComponent& operator=(const PCComponent& other);

        /**
         * @brief Move assignment operator.
         */
        PCComponent& operator=(PCComponent&& other) noexcept;

        /**
         * @brief Get a copy of the PCComponentData object contained in the current object.
         * @return A copy of the PCComponentData member.
         */
        PCComponentData& data();

        /**
         * @brief
         * @param pressure
         * @param temperature
         * @return
         */
        PCPhases flash(PCProps::Utilities::Pressure pressure, PCProps::Utilities::Temperature temperature) const;

        /**
         * @brief
         * @param pressure
         * @param vaporFraction
         * @return
         */
        PCPhases flash(PCProps::Utilities::Pressure pressure, PCProps::Utilities::VaporFraction vaporFraction) const;

        /**
         * @brief
         * @param temperature
         * @param vaporFraction
         * @return
         */
        PCPhases flash(PCProps::Utilities::Temperature temperature, PCProps::Utilities::VaporFraction vaporFraction) const;

        /**
         * @brief
         * @param pressure
         * @param enthalpy
         * @return
         */
        PCPhases flash(PCProps::Utilities::Pressure pressure, PCProps::Utilities::Enthalpy enthalpy) const;

        /**
         * @brief
         * @param pressure
         * @param entropy
         * @return
         */
        PCPhases flash(PCProps::Utilities::Pressure pressure, PCProps::Utilities::Entropy entropy) const;

        /**
         * @brief
         * @param temperature
         * @param volume
         * @return
         */
        PCPhases flash(PCProps::Utilities::Temperature temperature, PCProps::Utilities::Volume volume) const;

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

        /**
         * @brief Get the name of the component.
         * @return A std::string with the component name.
         */
        const std::string& name() const;

        /**
         * @brief Get the formula of the component.
         * @return A std::string with the component name.
         */
        const std::string& formula() const;

        /**
         * @brief Get the CAS number of the component.
         * @return A std::string with the CAS number.
         */
        const std::string& casrn() const;

        /**
         * @brief Get the SMILES string for the component.
         * @return A std::string with the SMILES string.
         */
        const std::string& smiles() const;

        void computeViscosity(PCPhases& phases) const;

    };
}    // namespace PCProps

#endif    // PCPROPS_PCCOMPONENT_HPP
