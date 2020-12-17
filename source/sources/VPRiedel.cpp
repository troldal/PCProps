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
#include <tuple>
#include <nlohmann/json.hpp>

#include "headers/VPRiedel.hpp"

namespace {
    std::tuple<double, double, double, double> initialize(double tboil, double tcrit, double pcrit) {
        double tbr = tboil / tcrit;

        double psi = -35.0 + 36.0/tbr + 42.0 * std::log(tbr) - std::pow(tbr, 6);
        double alpha_c = (3.758 * 0.0838 * psi + std::log(pcrit/101325.0)) / (0.0838 * psi - std::log(tbr));
        double h = tbr * std::log(pcrit/101325.0)/(1.0 - tbr);

        double D = 0.0838 * (alpha_c - 3.758);
        double C = alpha_c - 42.0 * D;
        double B = -36.0 * D;
        double A = 35.0 * D;

        A = A - C * std::log(tcrit) + std::log(pcrit);
        B = B * tcrit;
        D = D * std::pow(1.0/tcrit, 6);

        return std::make_tuple(A, B, C, D);
    }
} // namespace

using namespace std;
using json = nlohmann::json;

namespace PCProps::VaporPressure {

    // ===== Constructor
    VPRiedel::VPRiedel(double tboil, double tcrit, double pcrit) {
        auto coefficients = initialize(tboil, tcrit, pcrit);

        m_A = std::get<0>(coefficients);
        m_B = std::get<1>(coefficients);
        m_C = std::get<2>(coefficients);
        m_D = std::get<3>(coefficients);
    }

    // ===== Constructor
    VPRiedel::VPRiedel(const string& jsonstring) {
        auto js = json::parse(jsonstring);
        auto coefficients = initialize(js["tboil"], js["tcrit"], js["pcrit"]);

        m_A = std::get<0>(coefficients);
        m_B = std::get<1>(coefficients);
        m_C = std::get<2>(coefficients);
        m_D = std::get<3>(coefficients);

    }

    // ===== Function Call Operator
    double VPRiedel::operator()(double temperature) const
    {
        return std::exp(m_A + m_B/ temperature + m_C * std::log(temperature) + m_D * std::pow(temperature, 6));
    }

    // ===== Coefficients
    std::string VPRiedel::coefficients() const
    {
        json coeff;
        coeff["A"] = m_A;
        coeff["B"] = m_B;
        coeff["C"] = 0.0;
        coeff["D"] = 0.0;
        coeff["E"] = m_C;
        coeff["F"] = m_D;
        coeff["G"] = 6;

        return coeff.dump(4);
    }

} // namespace PCProps::VaporPressure