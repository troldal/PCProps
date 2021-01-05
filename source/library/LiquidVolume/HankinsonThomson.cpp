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

#include "HankinsonThomson.hpp"

namespace PCProps::LiquidVolume
{
    // ===== Constructor, taking critical properties as arguments
    HankinsonThomson::HankinsonThomson(double criticalTemperature, double characteristicVolume, double acentricFactor)
        : m_criticalTemperature(criticalTemperature),
          m_characteristicVolume(characteristicVolume),
          m_acentricFactor(acentricFactor)

    {}

    // =====  Copy constructor
    HankinsonThomson::HankinsonThomson(const HankinsonThomson& other) = default;

    // ===== Move constructor
    HankinsonThomson::HankinsonThomson(HankinsonThomson&& other) noexcept = default;

    // ===== Destructor
    HankinsonThomson::~HankinsonThomson() = default;

    // ===== Copy assignment operator
    HankinsonThomson& HankinsonThomson::operator=(const HankinsonThomson& other) = default;

    // ===== Move assignment operator
    HankinsonThomson& HankinsonThomson::operator=(HankinsonThomson&& other) noexcept = default;

    // ===== Function call operator, taking temperature [K] as argument and returns saturated liquid volume [m3/mol]
    double HankinsonThomson::operator()(double temperature) const
    {
        using std::pow;
        double tr = temperature / m_criticalTemperature;

        double V_0 = 1 - 1.528160 * pow(1 - tr, 1.0 / 3.0) + 1.439070 * pow(1 - tr, 2.0 / 3.0) - 0.814460 * (1 - tr) + 0.190454 * pow(1 - tr, 4.0 / 3.0);

        double V_Delta = (-0.296123 + 0.386914 * tr - 0.0427258 * pow(tr, 2) - 0.0480645 * pow(tr, 3)) / (tr - 1.00001);

        return m_characteristicVolume * V_0 * (1 - m_acentricFactor * V_Delta);
    }

    // ===== Static factory function,creating an SLVHankinsonThomson object from the characteristic volume.
    HankinsonThomson HankinsonThomson::createFromCharacteristicVolume(double criticalTemperature, double characteristicVolume, double acentricFactor)
    {
        return HankinsonThomson(criticalTemperature, characteristicVolume, acentricFactor);
    }

    // ===== Static factory function,creating an SLVHankinsonThomson object from estimate of the characteristic volume.
    HankinsonThomson HankinsonThomson::createFromEstimatedProperties(double criticalTemperature, double criticalPressure, double acentricFactor, HankinsonThomson::FluidType type)
    {
        using std::pow;
        double a = 0.2851686;
        double b = -0.06379110;
        double c = 0.01379173;

        switch (type) {
            case FluidType::Paraffin:
                a = 0.2905331;
                b = -0.0857958;
                c = 0.02276965;
                break;

            case FluidType::Olefin:
                a = 0.3070619;
                b = -0.2368581;
                c = 0.2834693;
                break;

            case FluidType::Cycloparaffin:
                a = 0.6564296;
                b = -3.391715;
                c = 7.442388;
                break;

            case FluidType::Aromatic:
                a = 0.2717636;
                b = -0.05759377;
                c = 0.05527757;
                break;

            case FluidType::SulfurCompound:
                a = 0.3053426;
                b = -0.1703247;
                c = 0.1753972;
                break;

            case FluidType::FluoroCarbon:
                a = 0.5218098;
                b = -2.349616;
                c = 5.407302;
                break;

            case FluidType::CryogenicLiquid:
                a = 0.2960998;
                b = -0.05468500;
                c = -0.1901563;
                break;

            case FluidType::CondensableGas:
                a = 0.2828447;
                b = -0.1183987;
                c = 0.1050570;
                break;

            default:    // ===== FluidType::Hydrocarbons
                break;
        }

        double V_char = (8.31446261815324 * criticalTemperature / criticalPressure) * (a + b * acentricFactor + c * pow(acentricFactor, 2));

        return HankinsonThomson(criticalTemperature, V_char, acentricFactor);
    }

}    // namespace PCProps::LiquidVolume
