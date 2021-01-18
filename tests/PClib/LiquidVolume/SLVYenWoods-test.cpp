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

#include <LiquidVolume/YenWoods.hpp>
#include <catch/catch.hpp>

using PCProps::LiquidVolume::YenWoods;

TEST_CASE("SLVYenWoods produces correct saturated liquid volume calculations")
{
    SECTION("Example from Caleb Bell GitHub repo")
    {
        auto water = YenWoods(YenWoods::CreateFromYenWoodsEstimation{647.14, 55.45E-6, 0.245});

        REQUIRE(water(300.0) == Approx(1.7695330765295693e-05));
    }

    SECTION("PPDS Coefficients: Test data from VDI Heat Map")
    {
        auto water = YenWoods(YenWoods::CreateFromPPDSCoefficients{647.14, 55.9472E-6, 18.02, 1094.0233, -1813.2295, 3863.9557, -2479.813});

        REQUIRE(water(300.0) == Approx(1.80811e-05));
    }

    SECTION("DIPPR 116 Coefficients: Test results calculated in Excel")
    {
        auto water = YenWoods(YenWoods::CreateFromDIPPR116Coefficients{647.14, 55.9472E-6, 58.606, -95.396, 213.89, -141.26});

        REQUIRE(water(300.0) == Approx(1.81202E-05));
    }

    SECTION("Copy constructor")
    {
        auto water = YenWoods(YenWoods::CreateFromPPDSCoefficients{647.14, 55.9472E-6, 18.02, 1094.0233, -1813.2295, 3863.9557, -2479.813});
        auto copy(water);

        REQUIRE(copy(300) == Approx(1.80811e-05));
    }

    SECTION("Move constructor")
    {
        auto water = YenWoods(YenWoods::CreateFromPPDSCoefficients{647.14, 55.9472E-6, 18.02, 1094.0233, -1813.2295, 3863.9557, -2479.813});
        auto copy(std::move(water));

        REQUIRE(copy(300) == Approx(1.80811e-05));
    }

    SECTION("Coefficients from Perry's using copy assignment")
    {
        auto water = YenWoods(YenWoods::CreateFromPPDSCoefficients{647.14, 55.9472E-6, 18.02, 1094.0233, -1813.2295, 3863.9557, -2479.813});
        auto copy  = YenWoods(YenWoods::CreateFromModifiedYenWoodsCoefficients{0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
        copy       = water;

        REQUIRE(copy(300) == Approx(1.80811e-05));
    }

    SECTION("Coefficients from Perry's using move assignment")
    {
        auto water = YenWoods(YenWoods::CreateFromPPDSCoefficients{647.14, 55.9472E-6, 18.02, 1094.0233, -1813.2295, 3863.9557, -2479.813});
        auto copy  = YenWoods(YenWoods::CreateFromOriginalYenWoodsCoefficients{0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
        copy       = std::move(water);

        REQUIRE(copy(300) == Approx(1.80811e-05));
    }
}