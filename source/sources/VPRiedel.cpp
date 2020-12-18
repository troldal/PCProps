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

#include "VPRiedel.hpp"

using namespace std;

namespace PCProps::VaporPressure {

    // ===== Constructor, default
    VPRiedel::VPRiedel() = default;

    // ===== Constructor
    VPRiedel::VPRiedel(double boilingTemperature, double criticalTemperature, double criticalPressure, VPRiedelType type)
        : m_criticalTemperature {criticalTemperature},
          m_criticalPressure {criticalPressure}
    {
        double tbr = boilingTemperature / criticalTemperature;
        double h = tbr * std::log(criticalPressure/101325.0)/(1.0 - tbr);
        double K = [&]() {
            switch (type) {
                case VPRiedelType::Organic:
                    return 0.0838;

                case VPRiedelType::Acid:
                    return -0.12 + 0.025 * h;

                case VPRiedelType::Alcohol:
                    return 0.373 - 0.030 * h;
            }
        }();

        double psi = -35.0 + 36.0/tbr + 42.0 * std::log(tbr) - std::pow(tbr, 6);
        double alpha_c = (3.758 * K * psi + std::log(criticalPressure/101325.0)) / (K * psi - std::log(tbr));


        m_coefficients[3] = K * (alpha_c - 3.758);
        m_coefficients[2] = alpha_c - 42.0 * m_coefficients[3];
        m_coefficients[1] = -36.0 * m_coefficients[3];
        m_coefficients[0] = 35.0 * m_coefficients[3];
    }

    // ===== Copy constructor
    VPRiedel::VPRiedel(const VPRiedel& other) = default;

    // ===== Move constructor
    VPRiedel::VPRiedel(VPRiedel&& other) noexcept = default;

    // ===== Destructor
    VPRiedel::~VPRiedel() = default;

    // ===== Copy assignment operator
    VPRiedel& VPRiedel::operator=(const VPRiedel& other) = default;

    // ===== Move assignment operator
    VPRiedel& VPRiedel::operator=(VPRiedel&& other) noexcept = default;

    // ===== Function Call Operator
    double VPRiedel::operator()(double temperature) const
    {
        using std::log;
        using std::pow;
        auto tr = temperature / m_criticalTemperature;
        return exp(m_coefficients[0] +
                   m_coefficients[1] / tr +
                   m_coefficients[2] * log(tr) +
                   m_coefficients[3] * pow(tr, 6)) * m_criticalPressure;
    }

    // ===== Get the critical temperature used in the vapor pressure estimation.
    double VPRiedel::criticalTemperature() const
    {
        return m_criticalTemperature;
    }

    // ===== Get the critical pressure used in the vapor pressure estimation.
    double VPRiedel::criticalPressure() const
    {
        return m_criticalPressure;
    }

    // ===== Coefficients
    std::array<double, 4> VPRiedel::coefficients() const
    {
        return m_coefficients;
    }

} // namespace PCProps::VaporPressure