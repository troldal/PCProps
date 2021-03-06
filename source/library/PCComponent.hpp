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
#include <optional>
#include <string>

#include <library/PCPropsException.hpp>

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
        std::string casrn; /** The CAS registration number for the component, e.g. 75-07-0 */
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

        std::function<double(double)> vaporPressureFunction {};             /**< The vapor pressure [Pa] as a function of temperature [K] */
        std::function<double(double)> liquidDensityFunction {};             /**< The liquid density [HOLD] as a function of temperature [K] */
        std::function<double(double)> surfaceTensionFunction {};            /**< The surface tension [HOLD] as a function of temperature [K] */
        std::function<double(double)> heatOfVaporizationFunction {};        /**< The latent heat [HOLD] as a function of temperature [K] */
        std::function<double(double)> vaporThermalConductivityFunction {};  /**< The vapor thermal conductivity [HOLD] as a function of temperature [K] */
        std::function<double(double)> liquidThermalConductivityFunction {}; /**< The liquid thermal conductivity [HOLD] as a function of temperature [K] */
        std::function<double(double)> vaporViscosityFunction {};            /**< The vapor viscosity [HOLD] as a function of temperature [K] */
        std::function<double(double)> liquidViscosityFunction {};           /**< The liquid viscosity [HOLD] as a function of temperature [K] */
        std::function<double(double)> idealGasCpFunction {};                /**< The ideal gas Cp [HOLD] as a function of temperature [K] */
        std::function<double(double)> liquidCpFunction {};                  /**< The liquid Cp [HOLD] as a function of temperature [K] */
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
        PCComponentData data() const;

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

        /**
         * @brief Check if the molecular weight has been set.
         * @return true, if the value has been set; otherwise false.
         */
        bool hasMolecularWeight() const;

        /**
         * @brief Get the molecular weight [g/mol] of the component.
         * @return A double with the molecular weight.
         * @throws PCPropsException if the molecularWeight data member has not been set.
         */
        double molecularWeight() const;

        /**
         * @brief Check if the normal boiling point temperature has been set.
         * @return true, if the value has been set; otherwise false.
         */
        bool hasBoilingTemperature() const;

        /**
         * @brief Get the normal boiling point temperature [K] of the component.
         * @return A double with the normal boiling temperature.
         * @throws PCPropsException if the boilingTemperature data member has not been set.
         */
        double boilingTemperature() const;

        /**
         * @brief Check if the freezing point temperature has been set.
         * @return true, if the value has been set; otherwise false.
         */
        bool hasFreezingTemperature() const;

        /**
         * @brief Get the freezing temperature [K] of the component.
         * @return A double with the freezing temperature.
         * @throws PCPropsException if the freezingTemperature data member has not been set.
         */
        double freezingTemperature() const;

        /**
         * @brief Check if the critical temperature has been set.
         * @return true, if the value has been set; otherwise false.
         */
        bool hasCriticalTemperature() const;

        /**
         * @brief Get the critical temperature [K] of the component.
         * @return A double with the critical temperature.
         * @throws PCPropsException if the criticalTemperature data member has not been set.
         */
        double criticalTemperature() const;

        /**
         * @brief Check if the critical pressure has been set.
         * @return true, if the value has been set; otherwise false.
         */
        bool hasCriticalPressure() const;

        /**
         * @brief Get the critical pressure [Pa] of the component.
         * @return A double with the critical pressure.
         * @throws PCPropsException if the criticalPressure data member has not been set.
         */
        double criticalPressure() const;

        /**
         * @brief Check if the critical volume has been set.
         * @return true, if the value has been set; otherwise false.
         */
        bool hasCriticalVolume() const;

        /**
         * @brief Get the critical molar volume [m3/mol] of the component.
         * @return A double with the critical volume.
         * @throws PCPropsException if the criticalVolume data member has not been set.
         */
        double criticalVolume() const;

        /**
         * @brief Check if the critical density has been set.
         * @return true, if the value has been set; otherwise false.
         */
        bool hasCriticalDensity() const;

        /**
         * @brief Get the critical density [kg/m3] of the component.
         * @return A double with the critical density.
         * @throws PCPropsException if the criticalDensity data member has not been set.
         */
        double criticalDensity() const;

        /**
         * @brief Check if the critical compressibility has been set.
         * @return true, if the value has been set; otherwise false.
         */
        bool hasCriticalCompressibility() const;

        /**
         * @brief Get the critical compressibility factor [-] of the component.
         * @return A double with the critical compressibility factor.
         * @throws PCPropsException if the criticalCompressibility data member has not been set.
         */
        double criticalCompressibility() const;

        /**
         * @brief Check if the acentric factor has been set.
         * @return true, if the value has been set; otherwise false.
         */
        bool hasAcentricFactor() const;

        /**
         * @brief Get the acentric factor [-] of the component.
         * @return A double with the acentric factor.
         * @throws PCPropsException if the acentricFactor data member has not been set.
         */
        double acentricFactor() const;

        /**
         * @brief Check if the vapor pressure function has been set.
         * @return true, if the function has been set; otherwise false.
         */
        bool hasVaporPressureFunction() const;

        /**
         * @brief Compute the vapor pressure [Pa] of the component at the given temperature [K].
         * @param temperature The temperature at which to evaluate the vapor pressure.
         * @return The vapor pressure.
         * @throws PCPropsException if the vaporPressureFunction data member has not been set.
         */
        double vaporPressure(double temperature) const;

        /**
         * @brief Check if the liquid density function has been set.
         * @return true, if the function has been set; otherwise false.
         */
        bool hasLiquidDensityFunction() const;

        /**
         * @brief Compute the liquid density [HOLD] of the component at the given temperature [K].
         * @param temperature The temperature at which to evaluate the liquid density.
         * @return The liquid density.
         * @throws PCPropsException if the liquidDensityFunction data member has not been set.
         */
        double liquidDensity(double temperature) const;

        /**
         * @brief Check if the surface tension function has been set.
         * @return true, if the function has been set; otherwise false.
         */
        bool hasSurfaceTensionFunction() const;

        /**
         * @brief Compute the surface tension [HOLD] of the component at the given temperature [K].
         * @param temperature The temperature at which to evaluate the surface tension.
         * @return The surface tension.
         * @throws PCPropsException if the surfaceTensionFunction data member has not been set.
         */
        double surfaceTension(double temperature) const;

        /**
         * @brief Check if the heat of vaporization function has been set.
         * @return true, if the function has been set; otherwise false.
         */
        bool hasHeatOfVaporizationFunction() const;

        /**
         * @brief Compute the latent heat [HOLD] of the component at the given temperature [K].
         * @param temperature The temperature at which to evaluate the latent heat.
         * @return the latent heat (heat of vaporization).
         * @throws PCPropsException if the heatOfVaporizationFunction data member has not been set.
         */
        double heatOfVaporization(double temperature) const;

        /**
         * @brief Check if the vapor thermal conductivity function has been set.
         * @return true, if the function has been set; otherwise false.
         */
        bool hasVaporThermalConductivityFunction() const;

        /**
         * @brief Compute the vapor thermal conductivity [HOLD] of the component at the given temperature [K].
         * @param temperature The temperature at which to evaluate the vapor thermal conductivity.
         * @return The vapor thermal conductivity.
         * @throws PCPropsException if the vaporThermalConductivityFunction data member has not been set.
         */
        double vaporThermalConductivity(double temperature) const;

        /**
         * @brief Check if the liquid thermal conductivity function has been set.
         * @return true, if the function has been set; otherwise false.
         */
        bool hasLiquidThermalConductivityFunction() const;

        /**
         * @brief Compute the liquid thermal conductivity [HOLD] of the component at the given temperature [K].
         * @param temperature The temperature at which to evaluate the liquid thermal conductivity.
         * @return The liquid thermal conductivity.
         * @throws PCPropsException if the liquidThermalConductivityFunction data member has not been set.
         */
        double liquidThermalConductivity(double temperature) const;

        /**
         * @brief Check if the vapor viscosity function has been set.
         * @return true, if the function has been set; otherwise false.
         */
        bool hasVaporViscosityFunction() const;

        /**
         * @brief Compute the vapor viscosity [HOLD] of the component at the given temperature [K].
         * @param temperature The temperature at which to evaluate the vapor viscosity.
         * @return The vapor viscosity.
         * @throws PCPropsException if the vaporViscosityFunction data member has not been set.
         */
        double vaporViscosity(double temperature) const;

        /**
         * @brief Check if the liquid viscosity function has been set.
         * @return true, if the function has been set; otherwise false.
         */
        bool hasLiquidViscosityFunction() const;

        /**
         * @brief Compute the liquid viscosity [HOLD] of the component at the given temperature [K].
         * @param temperature The temperature at which to evaluate the liquid viscosity.
         * @return The liquid viscosity.
         * @throws PCPropsException if the liquidViscosityFunction data member has not been set.
         */
        double liquidViscosity(double temperature) const;

        /**
         * @brief Check if the ideal gas Cp function has been set.
         * @return true, if the function has been set; otherwise false.
         */
        bool hasIdealGasCpFunction() const;

        /**
         * @brief Compute the ideal gas Cp [HOLD] of the component at the given temperature [K].
         * @param temperature The temperature at which to evaluate the ideal gas Cp.
         * @return The ideal gas Cp.
         * @throws PCPropsException if the idealGasCpFunction data member has not been set.
         */
        double idealGasCp(double temperature) const;

        /**
         * @brief Check if the liquid Cp function has been set.
         * @return true, if the function has been set; otherwise false.
         */
        bool hasLiquidCpFunction() const;

        /**
         * @brief Compute the liquid Cp [HOLD] of the component at the given temperature [K].
         * @param temperature The temperature at which to evaluate the liquid Cp.
         * @return The liquid Cp.
         * @throws PCPropsException if the liquidCpFunction data member has not been set.
         */
        double liquidCp(double temperature) const;
    };
}    // namespace PCProps

#endif    // PCPROPS_PCCOMPONENT_HPP
