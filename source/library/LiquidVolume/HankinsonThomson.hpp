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

#ifndef PCPROPS_HANKINSONTHOMSON_HPP
#define PCPROPS_HANKINSONTHOMSON_HPP

namespace PCProps::LiquidVolume
{
    /**
     * @brief The SLVHankinsonThomson class encapsulates the Hankinson & Thomson method for estimating
     * the saturated molar volume for liquids.
     * @details The SLVHankinsonThomson class encapsulates the Hankinson & Thomson method for estimating
     * the saturated molar volume for liquids:
     * \f[ V_{s} = V^* \cdot V^{(0)} \cdot \left [    1 - \omega_{SRK} \cdot V^{(\delta)}    \right] \f]
     * \f[ V^{(0)} = 1 +
     *  A \cdot \left(1-\frac{T}{T_{c}}\right)^{\frac{1}{3}} +
     *  B \cdot \left(1-\frac{T}{T_{c}}\right)^{\frac{2}{3}} +
     *  C \cdot \left(1-\frac{T}{T_{c}}\right) +
     *  D \cdot \left(1-\frac{T}{T_{c}}\right)^{\frac{4}{3}} \f]
     * \f[ V^{(\delta)} = \frac{E +
     * F \cdot \left(\frac{T}{T_{c}}\right) +
     * G \cdot \left(    \frac{T}{T_{c}}   \right)^2 +
     * H \cdot \left(    \frac{T}{T_{c}}   \right)^3}{\frac{T}{T_{c}} - 1.00001} \f]
     *
     * The Hankinson-Thomson method is also known as the COSTALD method. It is a popular method, due to it's accuracy
     * and ease of use. It depends on the critical temperature, a parameter called the characteristic volume and the
     * acentric factor adjusted to give the SRK equation of state the best fit to existing vapor pressure data. However,
     * several authors report that if those parameters are not available, the critical volume and the true acentric factor
     * may be used with little loss of accuracy. Alternatively, the characteristic volume can be estimated using the
     * following correlation:
     * \f[ V^* = \frac{R \cdot T_{c}}{P_{c}}\left(a + b \cdot \omega_{SRK} + c \cdot \omega_{SRK}^2 \right) \f]
     *
     * The coefficient a, b and c are constants specific for the component type and can be found in Reid et. al (4th edition)
     * Again, if the \f$ omega_{SRK} \f$ parameter is not available, the true acentric factor can be used.
     *
     * It is obvious that the COSTALD method is based on the Yen-Woods equation. The COSTALD equation can be rearranged
     * to resemble the Yen-Woods equation more closely. However, this would essentially make the Yen-Woods coefficients
     * temperature dependent. For this reason, the COSTALD method has to be implemented separately from the SLVYenWoods class.
     *
     * Note that the coefficients can only be set in the constructor or using the factory functions;
     * if coefficients needs to be changed after construction, a new object has to be created.
     *
     */
    class HankinsonThomson final
    {
        double m_criticalTemperature  = 0.0;
        double m_characteristicVolume = 0.0;
        double m_acentricFactor       = 0.0;


    public:
        /**
         * @brief Constructor, taking critical temperature [K], characteristic volume [m3/mol] and SRK acentric factor [-]
         * @details This constructor is private. To create an SLVHankinsonThomson object, use one of the static factory functions.
         * @param criticalTemperature The critical temperature [K]
         * @param characteristicVolume The characteristic volume [m3/mol]
         * @param acentricFactor The SRK acentric factor [-]
         * @note If the characteristic volume and/or the SRK acentric factor is not available, the critical volume and
         * the true acentric factor can be used with little loss of accuracy.
         */
        HankinsonThomson(double criticalTemperature, double characteristicVolume, double acentricFactor);

        /**
         * @brief An enum class enumerating the different fluid types available when estimating the characteristic volume.
         */
        enum class FluidType {
            Paraffin,
            Olefin,
            Cycloparaffin,
            Aromatic,
            Hydrocarbon,
            SulfurCompound,
            FluoroCarbon,
            CryogenicLiquid,
            CondensableGas
        };

        /**
         * @brief Copy constructor.
         */
        HankinsonThomson(const HankinsonThomson& other);

        /**
         * @brief Move constructor
         */
        HankinsonThomson(HankinsonThomson&& other) noexcept;

        /**
         * @brief Destructor
         */
        ~HankinsonThomson();

        /**
         * @brief Copy assignment operator.
         */
        HankinsonThomson& operator=(const HankinsonThomson& other);

        /**
         * @brief Move assignment operator.
         */
        HankinsonThomson& operator=(HankinsonThomson&& other) noexcept;

        /**
         * @brief Function call operator, taking temperature [K] as an argument and returns the liquid molar volume [m3/mol]
         * @param temperature The temperature [K]
         * @return The saturated liquid molar volume [m3/mol]
         */
        double operator()(double temperature) const;

        /**
         * @brief Static factory function, creating an SLVHankinsonThomson object from the critical temperature [K],
         * characteristic volume [m3/mol] and SRK acentric factor [-]
         * @param criticalTemperature The critical temperature [K]
         * @param characteristicVolume The characteristic volume [m3/mol]
         * @param acentricFactor The SRK acentric factor [-]
         * @note If the characteristic volume and/or the SRK acentric factor is not available, the critical volume and
         * the true acentric factor can be used with little loss of accuracy.
         * @return An SLVHankinsonThomson created from the input parameters.
         */
        static HankinsonThomson createFromCharacteristicVolume(double criticalTemperature, double characteristicVolume, double acentricFactor);

        /**
         * @brief Static factory function, creating an SLVHankinsonThomson object from the critical temperature [K],
         * critical pressure [Pa] and SRK acentric factor [-]
         * @details This method uses the critical properties to estimate the characteristic volume of fluid with the given type.
         * @param criticalTemperature The critical temperature [K]
         * @param criticalPressure The critical pressure [Pa]
         * @param acentricFactor The SRK acentric factor [-]
         * @param type The type of fluid. The default fluid type is an average hydrocarbon.
         * @note If the SRK acentric factor is not available, the true acentric factor can be used with little loss of accuracy.
         * @return An SLVHankinsonThomson created from critical properties and en estimated characteristic volume.
         */
        static HankinsonThomson
            createFromEstimatedProperties(double criticalTemperature, double criticalPressure, double acentricFactor, HankinsonThomson::FluidType type = FluidType::Hydrocarbon);
    };

}    // namespace PCProps::LiquidVolume
#endif    // PCPROPS_HANKINSONTHOMSON_HPP
