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

#include "SLVYenWoods.hpp"

namespace PCProps::LiquidVolume
{
    // ===== Constructor, taking critical properties and Yen-Woods coefficients as arguments
    SLVYenWoods::SLVYenWoods(
        double            criticalTemperature,
        double            criticalVolume,
        double            coeffA,
        double            coeffB,
        double            coeffC,
        double            coeffD,
        SLVYenWoods::Form type)
        : m_type(type),
          m_criticalTemperature(criticalTemperature),
          m_criticalVolume(criticalVolume),
          m_A(coeffA),
          m_B(coeffB),
          m_C(coeffC),
          m_D(coeffD)
    {}

    // ===== Copy constructor
    SLVYenWoods::SLVYenWoods(const SLVYenWoods& other) = default;

    // ===== Move constructor
    SLVYenWoods::SLVYenWoods(SLVYenWoods&& other) noexcept = default;

    // ===== Destructor
    SLVYenWoods::~SLVYenWoods() = default;

    // ===== Copy assignment operator
    SLVYenWoods& SLVYenWoods::operator=(const SLVYenWoods& other) = default;

    // ===== Move assignment operator
    SLVYenWoods& SLVYenWoods::operator=(SLVYenWoods&& other) noexcept = default;

    // ===== Function call operator, taking temperature [K] as and argument and returns the liquid molar volume [m3/mol]
    double SLVYenWoods::operator()(double temperature)
    {
        using std::pow;
        auto tau = 1 - (temperature / m_criticalTemperature);

        if (m_type == Form::Original)
            return m_criticalVolume / (1 + m_A * pow(tau, 1.0 / 3.0) + m_B * pow(tau, 2.0 / 3.0) + m_C * tau + m_D * pow(tau, 4.0 / 3.0));

        return m_criticalVolume / (1 + m_A * pow(tau, 0.35) + m_B * pow(tau, 2.0 / 3.0) + m_C * tau + m_D * pow(tau, 4.0 / 3.0));
    }

    // ===== Static factory function for creating an SLVYenWoods object from coefficients corresponding to the original form.
    SLVYenWoods SLVYenWoods::createFromOriginalYenWoodsCoefficients(
        double criticalTemperature,
        double criticalVolume,
        double coeffA,
        double coeffB,
        double coeffC,
        double coeffD)
    {
        return SLVYenWoods(criticalTemperature, criticalVolume, coeffA, coeffB, coeffC, coeffD, Form::Original);
    }

    // ===== Static factory function for creating an SLVYenWoods object from coefficients corresponding to the modified form.
    SLVYenWoods SLVYenWoods::createFromModifiedYenWoodsCoefficients(
        double criticalTemperature,
        double criticalVolume,
        double coeffA,
        double coeffB,
        double coeffC,
        double coeffD)
    {
        return SLVYenWoods(criticalTemperature, criticalVolume, coeffA, coeffB, coeffC, coeffD, Form::Modified);
    }

    // ===== Static factory function for creating an SLVYenWoods object from PPDS coefficients.
    SLVYenWoods SLVYenWoods::createFromPPDSCoefficients(
        double criticalTemperature,
        double criticalVolume,
        double molecularWeight,
        double coeffA,
        double coeffB,
        double coeffC,
        double coeffD)
    {
        double A = coeffA * 1000 * criticalVolume / molecularWeight;
        double B = coeffB * 1000 * criticalVolume / molecularWeight;
        double C = coeffC * 1000 * criticalVolume / molecularWeight;
        double D = coeffD * 1000 * criticalVolume / molecularWeight;

        return SLVYenWoods(criticalTemperature, criticalVolume, A, B, C, D, Form::Modified);
    }

    // ===== Static factory function for creating an SLVYenWoods object from DIPPR 116 coefficients.
    SLVYenWoods SLVYenWoods::createFromDIPPR116Coefficients(
        double criticalTemperature,
        double criticalVolume,
        double coeffA,
        double coeffB,
        double coeffC,
        double coeffD)
    {
        double A = coeffA * 1000 * criticalVolume;
        double B = coeffB * 1000 * criticalVolume;
        double C = coeffC * 1000 * criticalVolume;
        double D = coeffD * 1000 * criticalVolume;

        return SLVYenWoods(criticalTemperature, criticalVolume, A, B, C, D, Form::Modified);
    }

    // ===== Static factory function for creating an SLVYenWoods object from the Yen & Woods estimation procedure.
    SLVYenWoods SLVYenWoods::createFromYenWoodsEstimation(double criticalTemperature, double criticalVolume, double criticalCompressibility)
    {
        using std::pow;

        double A = 17.4425 - 214.578 * criticalCompressibility + 989.625 * pow(criticalCompressibility, 2) -
                   1522.06 * pow(criticalCompressibility, 3);

        double B = [&]() {
            if (criticalCompressibility <= 0.26)
                return -3.28257 + 13.6377 * criticalCompressibility + 107.4844 * pow(criticalCompressibility, 2) -
                       384.211 * pow(criticalCompressibility, 3);

            return 60.2091 - 402.063 * criticalCompressibility + 501.0 * pow(criticalCompressibility, 2) +
                   641.0 * pow(criticalCompressibility, 3);
        }();

        double C = 0.0;

        double D = 0.93 - B;

        return SLVYenWoods(criticalTemperature, criticalVolume, A, B, C, D, Form::Original);
    }

}    // namespace PCProps::LiquidVolume