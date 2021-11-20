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

#ifndef PCPROPS_ELBRO_HPP
#define PCPROPS_ELBRO_HPP

#include <array>
#include <cmath>
#include <functional>
#include <tuple>
#include <utility>
#include <vector>

namespace PCProps::LiquidVolume::detail
{
    /**
     * @brief
     */
    using SLVElbroGroupDef = std::tuple<double, double, double>;

    /**
     * @brief
     */
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

}    // namespace PCProps::LiquidVolume::detail

namespace PCProps::LiquidVolume
{
    /**
     * @brief The SLVElbroGroup is an alias for a std::pair<int, int>, used to define a single Elbro group
     * for a component. The first element is the Elbro group index, and the second element is the count. For
     * example, Hexadecane consists of two -CH3 groups and 14 -CH2 groups. -CH3 has an index of 1 and -CH2
     * has an index 2. Hence the SLVElbroGroup definitions are <1, 2> and <2, 14>
     */
    using SLVElbroGroup = std::pair<int, int>;

    /**
     * @brief The SLVElbro class implements the Elbro group contribution method for estimating saturated liquid volumes.
     * @details The Elbro method can be used to estimate saturated liquid molar volumes.
     * Documentation of how to use can be found in Reid et. al (5th edition).
     *
     * The group definitions are tabulated below:
     * Index  | Group
     * :----- | :-----
     * 1  | CH3
     * 2  | CH2
     * 3  | CH
     * 4  | C
     * 5  | ACH
     * 6  | ACCH3
     * 7  | ACCH2
     * 8  | ACCH
     * 9  | ACC
     * 10 | CH2=
     * 11 | CH=
     * 12 | C=
     * 13 | CH2OH
     * 14 | CHOH
     * 15 | ACOH
     * 16 | CH3CO
     * 17 | CH2CO
     * 18 | CHCO
     * 19 | CHOH
     * 20 | CH3COO
     * 21 | CH2COO
     * 22 | CHCOO
     * 23 | COO
     * 24 | ACCOO
     * 25 | CH3O
     * 26 | CH2O
     * 27 | CHOH
     * 28 | COO
     * 29 | CH2Cl
     * 30 | CHCl
     * 31 | CCl
     * 32 | CHCl2
     * 33 | CCl3
     * 34 | ACCl
     * 35 | Si
     * 36 | SiO
     */
    class Elbro
    {
        std::vector<std::function<double(double)>> m_groups; /**< Collection of polynomials used to calculate group contributions */

        /**
         * @brief Constructor (private), taking a std::vector of function objects as parameter.
         * @param groups A std::vector of function objects for calculating individual group contributions.
         * @note This is a private constructor. To create a SLVElbro object, use the SLVElbro::create static factory function.
         */
        explicit Elbro(const std::vector<std::function<double(double)>>& groups) : m_groups { groups } {}

    public:

        /**
         * @brief
         */
        Elbro() = default;

        /**
         * @brief Constructor, creating an SLVElbro object, taking a std::vector of group IDs and counts as argument.
         * @param groups A std::vector of SLVElbroGroup's (aka std::pair<int, int>) with group IDs and counts.
         * @return An SLVElbro object.
         */
        explicit Elbro(const std::vector<SLVElbroGroup>& groups) {

            using std::get;
            using std::pow;

            std::vector<std::function<double(double)>> lambdas {};

            for (const auto& item : groups) {
                auto coeffs = detail::ElbroGroups.at(static_cast<uint64_t>(item.first - 1));
                auto func   = [=](double temperature) -> double {
                       return get<0>(coeffs) * item.second + get<1>(coeffs) * item.second / 1000 * temperature + get<2>(coeffs) * item.second / 100000 * pow(temperature, 2);
                };

                lambdas.emplace_back(func);
            }

            *this = std::move(Elbro(lambdas));
        }

        /**
         * @brief Constructor creating an SLVElbro object, taking a collection of group IDs and counts
         * as argument. The collection can be any container supporting iterators.
         * @details If using a map (e.g. std::map or std::unordered_map), the Elbro group ID is the key, and the count is the value.
         * @tparam Container The type of the container used
         * @param groups A collection of SLVElbroGroup's (aka std::pair<int, int>) with group IDs and counts.
         * @return An SLVElbro object.
         */
        template<typename Container>
        explicit Elbro(const Container& groups)
        {
            std::vector<SLVElbroGroup> groupvec;
            for (const auto& item : groups) groupvec.template emplace_back(item);

            *this = std::move(Elbro(groupvec));
        }

        /**
         * @brief Copy constructor
         */
        Elbro(const Elbro& other) = default;

        /**
         * @brief Move constructor
         */
        Elbro(Elbro&& other) noexcept = default;

        /**
         * @brief Destructor
         */
        ~Elbro() = default;

        /**
         * @brief Copy assignment operator
         */
        Elbro& operator=(const Elbro& other) = default;

        /**
         * @brief Move assignment operator
         */
        Elbro& operator=(Elbro&& other) noexcept = default;

        /**
         * @brief Function call operator, taking temperature [K] as argument and returns saturated liquid molar volume [m3/mol]
         * @param temperature The temperature [K]
         * @return Saturated liquid molar volume [m3/mol]
         */
        double operator()(double temperature) const
        {
            double result = 0.0;
            for (const auto& item : m_groups) result += item(temperature);
            return result / 1000000;
        }
    };
}    // namespace PCProps::LiquidVolume
#endif    // PCPROPS_ELBRO_HPP
