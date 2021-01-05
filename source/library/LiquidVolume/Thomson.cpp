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

#include "Thomson.hpp"

namespace PCProps::LiquidVolume
{
    // ===== Constructor
    Thomson::Thomson(
        double                               criticalTemperature,
        double                               criticalPressure,
        double                               acentricFactor,
        const std::function<double(double)>& satVolumeFunction,
        const std::function<double(double)>& vaporPressureFunction)
        : m_saturatedVolumeFunction(satVolumeFunction),
          m_vaporPressureFunction(vaporPressureFunction),
          m_criticalTemperature(criticalTemperature),
          m_criticalPressure(criticalPressure),
          m_acentricFactor(acentricFactor)
    {}

    // ===== Copy constructor
    Thomson::Thomson(const Thomson& other) = default;

    // ===== Move constructor
    Thomson::Thomson(Thomson&& other) noexcept = default;

    // ===== Destructor
    Thomson::~Thomson() = default;

    // ===== Copy assignment operator
    Thomson& Thomson::operator=(const Thomson& other) = default;

    // ===== Move assignment operator
    Thomson& Thomson::operator=(Thomson&& other) noexcept = default;

    // ===== Function call operator
    double Thomson::operator()(double temperature, double pressure)
    {
        using std::exp;
        using std::log;
        using std::pow;

        double tr = temperature / m_criticalTemperature;
        double C  = 0.0861488 + 0.0344483 * m_acentricFactor;
        double B  = m_criticalPressure * (-1 - 9.070217 * pow(1 - tr, 1.0 / 3.0) + 62.45326 * pow(1 - tr, 2.0 / 3.0) - 135.1102 * (1 - tr) +
                                         exp(4.79594 + 0.250047 * m_acentricFactor + 1.14188 * pow(m_acentricFactor, 2)) * pow(1 - tr, 4.0 / 3.0));

        double psat = m_vaporPressureFunction(temperature);
        double vsat = m_saturatedVolumeFunction(temperature);

        return vsat * (1 - C * log((B + pressure) / (B + psat)));
    }

    // ===== Static factory function
    Thomson Thomson::create(
        double                               criticalTemperature,
        double                               criticalPressure,
        double                               acentricFactor,
        const std::function<double(double)>& satVolumeFunction,
        const std::function<double(double)>& vaporPressureFunction)
    {
        return Thomson(criticalTemperature, criticalPressure, acentricFactor, satVolumeFunction, vaporPressureFunction);
    }

}    // namespace PCProps::LiquidVolume