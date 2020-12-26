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

#include <array>
#include <cmath>
#include <numeric>
#include <tuple>

#include "SLVElbro.hpp"

namespace
{
    using SLVElbroGroupDef = std::tuple<double, double, double>;
    constexpr std::array<SLVElbroGroupDef, 36> ElbroGroups {
        SLVElbroGroupDef { 18.960, 45.58, 0.0 },    SLVElbroGroupDef { 12.520, 12.94, 0.0 },    SLVElbroGroupDef { 6.297, -21.92, 0.0 },
        SLVElbroGroupDef { 1.296, -59.66, 0.0 },    SLVElbroGroupDef { 10.090, 17.37, 0.0 },    SLVElbroGroupDef { 23.580, 24.43, 0.0 },
        SLVElbroGroupDef { 18.16, -8.589, 0.0 },    SLVElbroGroupDef { 8.925, -31.86, 0.0 },    SLVElbroGroupDef { 7.369, -83.60, 0.0 },
        SLVElbroGroupDef { 20.63, 31.43, 0.0 },     SLVElbroGroupDef { 6.761, 23.97, 0.0 },     SLVElbroGroupDef { -0.3971, -14.10, 0.0 },
        SLVElbroGroupDef { 39.46, -110.60, 23.31 }, SLVElbroGroupDef { 40.92, -193.20, 32.21 }, SLVElbroGroupDef { 41.20, -164.20, 22.78 },
        SLVElbroGroupDef { 42.18, -67.17, 22.58 },  SLVElbroGroupDef { 48.56, -170.40, 32.15 }, SLVElbroGroupDef { 25.17, -185.60, 28.59 },
        SLVElbroGroupDef { 12.090, 45.25, 0.0 },    SLVElbroGroupDef { 42.82, -20.50, 16.42 },  SLVElbroGroupDef { 49.73, -154.10, 33.19 },
        SLVElbroGroupDef { 43.28, -168.70, 33.25 }, SLVElbroGroupDef { 14.23, 11.93, 0.0 },     SLVElbroGroupDef { 43.06, -147.20, 20.93 },
        SLVElbroGroupDef { 16.66, 74.31, 0.0 },     SLVElbroGroupDef { 14.41, 28.54, 0.0 },     SLVElbroGroupDef { 35.07, -199.7, 40.93 },
        SLVElbroGroupDef { 30.12, -247.3, 40.69 },  SLVElbroGroupDef { 25.29, 49.11, 0.0 },     SLVElbroGroupDef { 17.40, 27.24, 0.0 },
        SLVElbroGroupDef { 37.62, -179.1, 32.47 },  SLVElbroGroupDef { 36.45, 54.31, 0.0 },     SLVElbroGroupDef { 48.74, 65.53, 0.0 },
        SLVElbroGroupDef { 23.51, 9.303, 0.0 },     SLVElbroGroupDef { 86.71, -555.5, 97.9 },   SLVElbroGroupDef { 17.41, -22.18, 0.0 }
    };

}    // namespace

namespace PCProps::LiquidVolume
{
    SLVElbro::SLVElbro(const std::vector<std::function<double(double)>>& groups) : m_groups { groups } {}

    SLVElbro::SLVElbro(const SLVElbro& other) = default;

    SLVElbro::SLVElbro(SLVElbro&& other) noexcept = default;

    SLVElbro::~SLVElbro() = default;

    SLVElbro& SLVElbro::operator=(const SLVElbro& other) = default;

    SLVElbro& SLVElbro::operator=(SLVElbro&& other) noexcept = default;

    double SLVElbro::operator()(double temperature) const
    {
        double result = 0.0;
        for (const auto& item : m_groups) result += item(temperature);
        return result / 1000000;
    }

    SLVElbro SLVElbro::create(const std::vector<SLVElbroGroup>& groups)
    {
        using std::get;
        using std::pow;

        std::vector<std::function<double(double)>> lambdas {};

        for (const auto& item : groups) {
            auto coeffs = ElbroGroups.at(static_cast<uint64_t>(item.first - 1));
            auto func   = [=](double temperature) -> double {
                return get<0>(coeffs) * item.second + get<1>(coeffs) * item.second / 1000 * temperature +
                       get<2>(coeffs) * item.second / 100000 * pow(temperature, 2);
            };

            lambdas.emplace_back(func);
        }

        return SLVElbro(lambdas);
    }

}    // namespace PCProps::LiquidVolume