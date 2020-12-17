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

#include "VPAntoineExt.hpp"

namespace PCProps::VaporPressure
{
    // ===== Constructor, default
    VPAntoineExt::VPAntoineExt() = default;

    // ===== Constructor, taking coefficienta A-G as arguments
    VPAntoineExt::VPAntoineExt(double A, double B, double C, double D, double E, double F, double G) : m_coefficients { A, B, C, D, E, F, G } {}

    // ===== Copy constructor
    VPAntoineExt::VPAntoineExt(const VPAntoineExt& other) = default;

    // ===== Move constructor
    VPAntoineExt::VPAntoineExt(VPAntoineExt&& other) noexcept = default;

    // ===== Destructor
    VPAntoineExt::~VPAntoineExt() = default;

    // ===== Copy assignment operator
    VPAntoineExt& VPAntoineExt::operator=(const VPAntoineExt& other) = default;

    // ===== Move assignment operator
    VPAntoineExt& VPAntoineExt::operator=(VPAntoineExt&& other) noexcept = default;

    // ===== Function call operator
    double VPAntoineExt::operator()(double temperature) const
    {
        using std::log;
        using std::pow;
        using std::exp;
        return exp(m_coefficients[0] +
                   m_coefficients[1] / (temperature + m_coefficients[2]) +
                   m_coefficients[3] * temperature +
                   m_coefficients[4] * log(temperature) +
                   m_coefficients[5] * pow(temperature, static_cast<int>(m_coefficients[6])));
    }

    // ===== Return array of coefficients
    std::array<double, 7> VPAntoineExt::coefficients() const
    {
        return m_coefficients;
    }

}    // namespace PCProps::VaporPressure
