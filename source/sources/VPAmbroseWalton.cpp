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

#include "VPAmbroseWalton.hpp"

namespace PCProps::VaporPressure {

    // ===== Constructor, default
    VPAmbroseWalton::VPAmbroseWalton() = default;

    // ===== Constructor, taking critical properties and acentric factor as arguments
    VPAmbroseWalton::VPAmbroseWalton(double criticalTemperature, double criticalPressure, double acentricFactor)
        : m_criticalTemperature {criticalTemperature},
          m_criticalPressure {criticalPressure},
          m_acentricFactor {acentricFactor}
    {}

    // ===== Copy constructor
    VPAmbroseWalton::VPAmbroseWalton(const VPAmbroseWalton& other) = default;

    // ===== Move constructor
    VPAmbroseWalton::VPAmbroseWalton(VPAmbroseWalton&& other) noexcept = default;

    // ===== Destructor
    VPAmbroseWalton::~VPAmbroseWalton() = default;

    // ===== Copy assignment operator
    VPAmbroseWalton& VPAmbroseWalton::operator=(const VPAmbroseWalton& other) = default;

    // ===== Move assignment operator
    VPAmbroseWalton& VPAmbroseWalton::operator=(VPAmbroseWalton&& other) noexcept = default;

    // ===== Function call operator
    double VPAmbroseWalton::operator()(double temperature)
    {
        using std::pow;
        using std::exp;
        auto tau = 1 - (temperature / m_criticalTemperature);
        auto f0 = (-5.97616 * tau + 1.29874 * pow(tau, 1.5) - 0.60394 * pow(tau, 2.5) - 1.06841 * pow(tau, 5)) / (1 - tau);
        auto f1 = (-5.03365 * tau + 1.11505 * pow(tau, 1.5) - 5.41217 * pow(tau, 2.5) - 7.46628 * pow(tau, 5)) / (1 - tau);
        auto f2 = (-0.64771 * tau + 2.41539 * pow(tau, 1.5) - 4.26979 * pow(tau, 2.5) + 3.25259 * pow(tau, 5)) / (1 - tau);

        return exp(f0 + m_acentricFactor * f1 + pow(m_acentricFactor, 2) * f2) * m_criticalPressure;
    }

    // ===== Get the critical temperature used for the vapor pressure estimation
    double VPAmbroseWalton::criticalTemperature() const
    {
        return m_criticalTemperature;
    }

    // ===== Get the critical pressure used for the vapor pressure estimation
    double VPAmbroseWalton::criticalPressure() const
    {
        return m_criticalPressure;
    }

    // ===== Get the acentric factor used for the vapor pressure estimation
    double VPAmbroseWalton::acentricFactor() const
    {
        return m_acentricFactor;
    }

} // namespace PCProps::VaporPressure