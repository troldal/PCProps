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

        std::optional<double> molarWeight {};             /**< The molecular weight [g/mol] of a component, e.g. 16.043 g/mol */
        std::optional<double> boilingTemperature {};      /**< The normal boiling temperature [K] (at 101325 Pa) of a component */
        std::optional<double> freezingTemperature {};     /**< The freezing point temperature [K] of a component */
        std::optional<double> criticalTemperature {};     /**< The critical temperature [K] of a component */
        std::optional<double> criticalPressure {};        /**< The critical pressure [Pa] of a component */
        std::optional<double> criticalVolume {};          /**< The critical molar volume [m3/mol] of a component */
        std::optional<double> criticalDensity {};         /**< The critical density [kg/m3] of a component */
        std::optional<double> criticalCompressibility {}; /**< The critical compressibility factor [-] of a component */
        std::optional<double> acentricFactor {};          /**< The acentric factor (omega) of a component */
        std::optional<double> dipoleMoment {};            /**< The dipole moment of the component */

        std::function<double(double)> saturatedLiquidVolumeCorrelation {}; /**< The liquid density [HOLD] as a function of temperature [K] */
        std::function<double(double)> idealGasCpCorrelation {};            /**<  */
        std::function<double(double)> liquidCpCorrelation {};              /**< The liquid Cp [HOLD] as a function of temperature [K] */
        std::function<double(double)> vaporPressureCorrelation {};         /**< The vapor pressure [Pa] as a function of temperature [K] */
        std::function<double(double)> surfaceTensionCorrelation {};        /**< The surface tension [HOLD] as a function of temperature [K] */
        std::function<double(double)> heatOfVaporizationCorrelation {};    /**< The latent heat [HOLD] as a function of temperature [K] */

        std::function<double(double)> saturatedVaporThermalConductivityCorrelation {};  /**< The vapor thermal conductivity [HOLD] as a function of temperature [K] */
        std::function<double(double)> saturatedLiquidThermalConductivityCorrelation {}; /**< The liquid thermal conductivity [HOLD] as a function of temperature [K] */
        std::function<double(double)> saturatedVaporViscosityCorrelation {};            /**< The vapor viscosity [HOLD] as a function of temperature [K] */
        std::function<double(double)> saturatedLiquidViscosityCorrelation {};           /**< The liquid viscosity [HOLD] as a function of temperature [K] */

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

    public:

        // ===== Constructors & Assignment Operators ===== //

        /**
         * @brief Default constructor.
         */
        PCComponent() = default;

        /**
         * @brief Constructor, taking a PCComponentData object as an argument.
         * @param data A PCComponentData object with the component data.
         */
        explicit PCComponent(const PCComponentData& data) : m_data(data) {}

        /**
         * @brief Copy constructor.
         */
        PCComponent(const PCComponent& other) = default;

        /**
         * @brief Move constructor.
         */
        PCComponent(PCComponent&& other) noexcept = default;

        /**
         * @brief Destructor.
         */
        ~PCComponent() = default;

        /**
         * @brief Copy assignment operator.
         */
        PCComponent& operator=(const PCComponent& other) = default;

        /**
         * @brief Move assignment operator.
         */
        PCComponent& operator=(PCComponent&& other) noexcept = default;


        // ===== Accessors (Constants) ===== //

        /**
         * @brief Get a copy of the PCComponentData object contained in the current object.
         * @return A copy of the PCComponentData member.
         */
        PCComponentData& data()
        {
            return m_data;
        }

        /**
         * @brief Get the name of the component.
         * @return A std::string with the component name.
         */
        const std::string& name() const
        {
            return m_data.name;
        }

        /**
         * @brief Get the formula of the component.
         * @return A std::string with the component name.
         */
        const std::string& formula() const
        {
            return m_data.formula;
        }

        /**
         * @brief Get the CAS number of the component.
         * @return A std::string with the CAS number.
         */
        const std::string& casrn() const
        {
            return m_data.casrn;
        }

        /**
         * @brief Get the SMILES string for the component.
         * @return A std::string with the SMILES string.
         */
        const std::string& smiles() const
        {
            return m_data.smiles;
        }



        bool molecularWeightIsValid() const {
            return m_data.molarWeight.has_value();
        }

        double molarWeight() const {
            return m_data.molarWeight.value();
        }

        bool boilingTemperatureIsValid() const {
            return m_data.boilingTemperature.has_value();
        }

        double boilingTemperature() const {
            return m_data.boilingTemperature.value();
        }

        bool freezingTemperatureIsValid() const {
            return m_data.freezingTemperature.has_value();
        }

        double freezingTemperature() const {
            return m_data.freezingTemperature.value();
        }

        bool criticalTemperatureIsValid() const {
            return m_data.criticalTemperature.has_value();
        }

        double criticalTemperature() const {
            return m_data.criticalTemperature.value();
        }

        bool criticalPressureIsValid() const {
            return m_data.criticalPressure.has_value();
        }

        double criticalPressure() const {
            return m_data.criticalPressure.value();
        }

        bool criticalVolumeIsValid() const {
            return m_data.criticalVolume.has_value();
        }

        double criticalVolume() const {
            return m_data.criticalVolume.value();
        }

        bool criticalDensityIsValid() const {
            return m_data.criticalDensity.has_value();
        }

        double criticalDensity() const {
            return m_data.criticalDensity.value();
        }

        bool criticalCompressibilityIsValid() const {
            return m_data.criticalCompressibility.has_value();
        }

        double criticalCompressibility() const {
            return m_data.criticalCompressibility.value();
        }

        bool acentricFactorIsValid() const {
            return m_data.acentricFactor.has_value();
        }

        double acentricFactor() const {
            return m_data.acentricFactor.value();
        }

        bool dipoleMomentIsValid() const {
            return m_data.dipoleMoment.has_value();
        }

        double dipoleMoment() const {
            return m_data.dipoleMoment.value();
        }

        // ===== Accessors (Temperature Dependent) ===== //

        bool satLiquidVolumeIsValid() const {
            return static_cast<bool>(m_data.saturatedLiquidVolumeCorrelation);
        }

        double satLiquidVolume(double temperature) const {
            return m_data.saturatedLiquidVolumeCorrelation(temperature);
        }

        bool idealGasCpIsValid() const {
            return static_cast<bool>(m_data.idealGasCpCorrelation);
        }

        double idealGasCp(double temperature) const {
            return m_data.idealGasCpCorrelation(temperature);
        }

        bool satLiquidCpIsValid() const {
            return static_cast<bool>(m_data.liquidCpCorrelation);
        }

        double satLiquidCp(double temperature) const {
            return m_data.liquidCpCorrelation(temperature);
        }

        bool vaporPressureIsValid() const {
            return static_cast<bool>(m_data.vaporPressureCorrelation);
        }

        double vaporPressure(double temperature) const {
            return m_data.vaporPressureCorrelation(temperature);
        }

        bool surfaceTensionIsValid() const{
            return static_cast<bool>(m_data.surfaceTensionCorrelation);
        }

        double surfaceTension(double temperature) const {
            return m_data.surfaceTensionCorrelation(temperature);
        }

        bool heatOfVaporizationIsValid() const {
            return static_cast<bool>(m_data.heatOfVaporizationCorrelation);
        }

        double heatOfVaporization(double temperature) const {
            return m_data.heatOfVaporizationCorrelation(temperature);
        }

        bool satVaporThermalConductivityIsValid() const {
            return static_cast<bool>(m_data.saturatedVaporThermalConductivityCorrelation);
        }

        double satVaporThermalConductivity(double temperature) const {
            return m_data.saturatedVaporThermalConductivityCorrelation(temperature);
        }

        bool satLiquidThermalConductivityIsValid() const {
            return static_cast<bool>(m_data.saturatedLiquidThermalConductivityCorrelation);
        }

        double satLiquidThermalConductivity(double temperature) const {
            return m_data.saturatedLiquidThermalConductivityCorrelation(temperature);
        }

        bool satVaporViscosityIsValid() const {
            return static_cast<bool>(m_data.saturatedVaporViscosityCorrelation);
        }

        double satVaporViscosity(double temperature) const {
            return m_data.saturatedVaporViscosityCorrelation(temperature);
        }

        bool satLiquidViscosityIsValid() const {
            return static_cast<bool>(m_data.saturatedLiquidViscosityCorrelation);
        }

        double satLiquidViscosity(double temperature) const {
            return m_data.saturatedLiquidViscosityCorrelation(temperature);
        }



//        void computeViscosity(PCPhases& phases) const;

    };
}    // namespace PCProps

#endif    // PCPROPS_PCCOMPONENT_HPP
