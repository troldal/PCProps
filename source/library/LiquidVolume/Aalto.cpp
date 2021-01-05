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

#include "Aalto.hpp"

namespace PCProps::LiquidVolume
{
    // ===== Constructor
    Aalto::Aalto(
        double                               criticalTemperature,
        double                               criticalPressure,
        double                               acentricFactor,
        const std::function<double(double)>& satVolumeFunction,
        const std::function<double(double)>& vaporPressureFunction)
        : m_saturatedVolumeFunction(satVolumeFunction),
          m_vaporPressureFunction(vaporPressureFunction),
          m_criticalTemperature(criticalTemperature),
          m_criticalPressure(criticalPressure),
          B { 0.164813 - 0.0914427 * acentricFactor }
    {}

    // ===== Copy constructor
    Aalto::Aalto(const Aalto& other) = default;

    // ===== Move constructor
    Aalto::Aalto(Aalto&& other) noexcept = default;

    // ===== Destructor
    Aalto::~Aalto() = default;

    // ===== Copy assignment operator
    Aalto& Aalto::operator=(const Aalto& other) = default;

    // ===== Move assignment operator
    Aalto& Aalto::operator=(Aalto&& other) noexcept = default;

    // ===== Function call operator
    double Aalto::operator()(double temperature, double pressure) const
    {
        using std::pow;
        double tr = temperature / m_criticalTemperature;

        double A = -170.335 - 28.578 * tr + 124.809 * pow(tr, 3) - 55.5393 * pow(tr, 6) + 130.01 / tr;

        return m_saturatedVolumeFunction(temperature) * ((A * m_criticalPressure + pow(C, pow(D - tr, B)) * (pressure - m_vaporPressureFunction(temperature))) /
                                                         (A * m_criticalPressure + C * (pressure - m_vaporPressureFunction(temperature))));
    }

    // ===== Static factory function
    Aalto Aalto::create(
        double                               criticalTemperature,
        double                               criticalPressure,
        double                               acentricFactor,
        const std::function<double(double)>& satVolumeFunction,
        const std::function<double(double)>& vaporPressureFunction)
    {
        return Aalto(criticalTemperature, criticalPressure, acentricFactor, satVolumeFunction, vaporPressureFunction);
    }

}    // namespace PCProps::LiquidVolume