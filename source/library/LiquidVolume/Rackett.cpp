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

#include "Rackett.hpp"

namespace PCProps::LiquidVolume
{
    // ===== Constructor, default
    Rackett::Rackett() = default;

    // ===== Constructor, taking the four Rackett coefficients as arguments
    Rackett::Rackett(double coeffA, double coeffB, double coeffC, double coeffD) : m_A { coeffA }, m_B { coeffB }, m_C { coeffC }, m_D { coeffD } {}

    // ===== Constructor, taking DIPPR coefficients as arguments
    Rackett::Rackett(const Rackett::CreateFromDIPPR& c) : Rackett(1 / (1000 * c.A), c.B, c.C, c.D) {}

    Rackett::Rackett(const Rackett::CreateFromYaws& c) : Rackett(c.molecularWeight / (c.A * c.B * 1E6), c.B, c.Tc, c.n) {}

    // ===== Copy constructor
    Rackett::Rackett(const Rackett& other) = default;

    // ===== Move constructor
    Rackett::Rackett(Rackett&& other) noexcept = default;

    // ===== Destructor
    Rackett::~Rackett() = default;

    // ===== Copy assignment operator
    Rackett& Rackett::operator=(const Rackett& other) = default;

    // ===== Move assignment operator
    Rackett& Rackett::operator=(Rackett&& other) noexcept = default;

    // ===== Function call operator, taking temperature [K] as argument, returning saturated liquid density [mol/m3]
    double Rackett::operator()(double temperature) const
    {
        using std::pow;
        return m_A * pow(m_B, 1 + pow(1 - temperature / m_C, m_D));
    }

    // ===== Static factory, constructing an object from the four Rackett coefficients
    Rackett Rackett::createFromCoefficients(double coeffA, double coeffB, double coeffC, double coeffD)
    {
        return Rackett(coeffA, coeffB, coeffC, coeffD);
    }

    // ===== Static factory, constructing an object from critical prpoerties only
    Rackett Rackett::createFromCriticalProperties(double criticalTemperature, double criticalPressure, double criticalCompressibility)
    {
        return Rackett((8.31446261815324 * criticalTemperature) / criticalPressure, criticalCompressibility, criticalTemperature, 2.0 / 7.0);
    }

    // ===== Static factory, constructing object using the Yamada-Gunn relation
    Rackett Rackett::createFromAcentricFactor(double criticalTemperature, double criticalPressure, double acentricFactor)
    {
        return Rackett((8.31446261815324 * criticalTemperature) / criticalPressure, (0.29056 - 0.08775 * acentricFactor), criticalTemperature, 2.0 / 7.0);
    }

    // ===== Static factory, constructing an object from a known reference point and the acentric factor.
    Rackett Rackett::createFromReferencePointA(double criticalTemperature, double experimentalTemperature, double experimentalVolume, double acentricFactor)
    {
        using std::pow;
        auto k = pow(1.0 - experimentalTemperature / criticalTemperature, 2.0 / 7.0);
        auto z = (0.29056 - 0.08775 * acentricFactor);

        return Rackett(experimentalVolume / pow(z, 1 + k), z, criticalTemperature, 2.0 / 7.0);
    }

    // ===== Static factory, constructing an object from a known reference point and the critical compressibility.
    Rackett Rackett::createFromReferencePointB(double criticalTemperature, double experimentalTemperature, double experimentalVolume, double criticalCompressibility)
    {
        using std::pow;
        auto k = pow(1.0 - experimentalTemperature / criticalTemperature, 2.0 / 7.0);

        return Rackett(experimentalVolume / pow(criticalCompressibility, 1 + k), criticalCompressibility, criticalTemperature, 2.0 / 7.0);
    }

}    // namespace PCProps::LiquidVolume