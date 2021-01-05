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

#include "AntoineExtended.hpp"

namespace PCProps::VaporPressure
{
    // ===== Constructor, default
    AntoineExtended::AntoineExtended() = default;

    // ===== Constructor, taking coefficienta A-G as arguments
    AntoineExtended::AntoineExtended(double A, double B, double C, double D, double E, double F, double G)
        : m_coeffA { A },
          m_coeffB { B },
          m_coeffC { C },
          m_coeffD { D },
          m_coeffE { E },
          m_coeffF { F },
          m_coeffG { G }
    {}

    // ===== Constructor, creating an object from DIPPR coefficients
    AntoineExtended::AntoineExtended(const AntoineExtended::CreateFromDIPPR& coefficients)
        : AntoineExtended { coefficients.A, coefficients.B, 0.0, 0.0, coefficients.C, coefficients.D, coefficients.E }
    {}

    // ===== Constructor, creating an object from Yaws (1999) coefficients
    AntoineExtended::AntoineExtended(const AntoineExtended::CreateFromYaws& c)
        : AntoineExtended(log(133.322368) + c.A * log(10), c.B * log(10), 0.0, c.D * log(10), c.C, c.E * log(10), 2)
    {}

    // ===== Copy constructor
    AntoineExtended::AntoineExtended(const AntoineExtended& other) = default;

    // ===== Move constructor
    AntoineExtended::AntoineExtended(AntoineExtended&& other) noexcept = default;

    // ===== Destructor
    AntoineExtended::~AntoineExtended() = default;

    // ===== Copy assignment operator
    AntoineExtended& AntoineExtended::operator=(const AntoineExtended& other) = default;

    // ===== Move assignment operator
    AntoineExtended& AntoineExtended::operator=(AntoineExtended&& other) noexcept = default;

    // ===== Function call operator
    double AntoineExtended::operator()(double temperature) const
    {
        using std::exp;
        using std::log;
        using std::pow;
        return exp(m_coeffA + m_coeffB / (temperature + m_coeffC) + m_coeffD * temperature + m_coeffE * log(temperature) + m_coeffF * pow(temperature, static_cast<int>(m_coeffG)));
    }

}    // namespace PCProps::VaporPressure
