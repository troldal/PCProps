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
#include <iostream>

#include "SLVRackett.hpp"

namespace PCProps::LiquidVolume
{
    // ===== Constructor, default
    SLVRackett::SLVRackett() = default;

    // ===== Constructor, taking the four Rackett coefficients as arguments
    SLVRackett::SLVRackett(double coeffA, double coeffB, double coeffC, double coeffD)
        : m_A { coeffA },
          m_B { coeffB },
          m_C { coeffC },
          m_D { coeffD }
    {}

    // ===== Copy constructor
    SLVRackett::SLVRackett(const SLVRackett& other) = default;

    // ===== Move constructor
    SLVRackett::SLVRackett(SLVRackett&& other) noexcept = default;

    // ===== Destructor
    SLVRackett::~SLVRackett() = default;

    // ===== Copy assignment operator
    SLVRackett& SLVRackett::operator=(const SLVRackett& other) = default;

    // ===== Move assignment operator
    SLVRackett& SLVRackett::operator=(SLVRackett&& other) noexcept = default;

    // ===== Function call operator, taking temperature [K] as argument, returning saturated liquid density [mol/m3]
    double SLVRackett::operator()(double temperature) const
    {
        using std::pow;
        return m_A * pow(m_B, 1 + pow(1 - temperature / m_C, m_D));
    }

    // ===== Static factory, constructing an object from the four Rackett coefficients
    SLVRackett SLVRackett::createFromCoefficients(double coeffA, double coeffB, double coeffC, double coeffD)
    {
        return SLVRackett(coeffA, coeffB, coeffC, coeffD);
    }

    // ===== Static factory, constructing an object from critical prpoerties only
    SLVRackett SLVRackett::createFromCriticalProperties(double criticalTemperature, double criticalPressure, double criticalCompressibility)
    {
        return SLVRackett(
            (8.31446261815324 * criticalTemperature) / criticalPressure,
            criticalCompressibility,
            criticalTemperature,
            2.0 / 7.0);
    }

    // ===== Static factory, constructing object using the Yamada-Gunn relation
    SLVRackett SLVRackett::createFromAcentricFactor(double criticalTemperature, double criticalPressure, double acentricFactor)
    {
        return SLVRackett(
            (8.31446261815324 * criticalTemperature) / criticalPressure,
            (0.29056 - 0.08775 * acentricFactor),
            criticalTemperature,
            2.0 / 7.0);
    }

    // ===== Static factory, constructing an object from a known reference point and the acentric factor.
    SLVRackett SLVRackett::createFromReferencePointA(
        double criticalTemperature,
        double experimentalTemperature,
        double experimentalVolume,
        double acentricFactor)
    {
        using std::pow;
        auto k = pow(1.0 - experimentalTemperature / criticalTemperature, 2.0 / 7.0);
        auto z = (0.29056 - 0.08775 * acentricFactor);

        return SLVRackett(experimentalVolume / pow(z, 1 + k), z, criticalTemperature, 2.0 / 7.0);
    }

    // ===== Static factory, constructing an object from a known reference point and the critical compressibility.
    SLVRackett SLVRackett::createFromReferencePointB(
        double criticalTemperature,
        double experimentalTemperature,
        double experimentalVolume,
        double criticalCompressibility)
    {
        using std::pow;
        auto k = pow(1.0 - experimentalTemperature / criticalTemperature, 2.0 / 7.0);

        return SLVRackett(
            experimentalVolume / pow(criticalCompressibility, 1 + k),
            criticalCompressibility,
            criticalTemperature,
            2.0 / 7.0);
    }

}    // namespace PCProps::LiquidVolume