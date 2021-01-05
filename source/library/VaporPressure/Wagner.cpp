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

#include "Wagner.hpp"

namespace PCProps::VaporPressure {
    // ===== Constructor, default
    Wagner::Wagner() = default;

    // ===== Constructor, taking critical properties and coefficients as arguments
    Wagner::Wagner(double criticalTemperature, double criticalPressure, double A, double B, double C, double D, VPWagnerForm form)
        : m_criticalTemperature { criticalTemperature },
          m_criticalPressure { criticalPressure },
          m_coefficients { A, B, C, D },
          m_expC { form == VPWagnerForm::Form25 ? 2.5 : 3.0 },
          m_expD { form == VPWagnerForm::Form25 ? 5.0 : 6.0 }
    {}

    // ===== Copy constructor
    Wagner::Wagner(const Wagner& other) = default;

    // ===== Move constructor
    Wagner::Wagner(Wagner&& other) noexcept = default;

    // ===== Destructor
    Wagner::~Wagner() = default;

    // ===== Copy assignment operator
    Wagner& Wagner::operator=(const Wagner& other) = default;

    // ===== Move assignment operator
    Wagner& Wagner::operator=(Wagner&& other) noexcept = default;

    // ===== Function call operator
    double Wagner::operator()(double temperature) const
    {
        using std::exp;
        using std::log;
        using std::pow;
        auto tr = temperature / m_criticalTemperature;

        return exp((1.0 / tr) *
                   (m_coefficients[0] * (1 - tr) + m_coefficients[1] * pow(1 - tr, 1.5) + m_coefficients[2] * pow(1 - tr, m_expC) + m_coefficients[3] * pow(1 - tr, m_expD))) * m_criticalPressure;
    }

    // ===== Getter for the critical temperature [K]
    double Wagner::criticalTemperature() const
    {
        return m_criticalTemperature;
    }

    // ===== Getter for the critical pressure [Pa]
    double Wagner::criticalPressure() const
    {
        return m_criticalPressure;
    }

    // ===== Getter for the equation coefficients
    std::array<double, 4> Wagner::coefficients() const
    {
        return m_coefficients;
    }

    // ===== Getter for the form of the Wagner equation
    VPWagnerForm Wagner::form() const
    {
        if (m_expC == 3 and m_expD == 6) return VPWagnerForm::Form36;

        return VPWagnerForm::Form25;
    }

} // namespace PCProps::VaporPressure