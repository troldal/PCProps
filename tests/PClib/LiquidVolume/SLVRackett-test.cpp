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

#include <LiquidVolume/Rackett.hpp>
#include <catch/catch.hpp>

using PCProps::LiquidVolume::Rackett;

TEST_CASE("SLVRackett produces correct saturated liquid volume calculations")
{
    SECTION("Default constructed SLVRackett object")
    {
        auto rackett = Rackett();

        REQUIRE(rackett(300.0) == Approx(0.0));
    }

    SECTION("Example 4-5a from Poling et. al")
    {
        auto trifluoroethane = Rackett(Rackett::CreateFromCriticalProperties{346.3, 37.92E5, 0.255});

        // ===== The example in Poling et. al errorneously gives 89.726E-6 as the result
        REQUIRE(trifluoroethane(300.0) == Approx(89.7395E-6));
    }

    SECTION("Example 4-5b from Poling et. al")
    {
        auto trifluoroethane = Rackett(Rackett::CreateFromReferencePointA{346.3, 245, 75.38E-6, 0.259});

        // ===== The example in Poling et. al errorneously gives 90.59E-6 as the result
        REQUIRE(trifluoroethane(300.0) == Approx(90.7765E-6));
    }

    SECTION("Example 4-5 from Poling et. al, using Yamada-Gunn")
    {
        auto trifluoroethane = Rackett(Rackett::CreateFromYamadaGunn {346.3, 37.92E5, 0.259});

        REQUIRE(trifluoroethane(300.0) == Approx(96.8963E-6));
    }

    SECTION("Example 4-5 from Poling et. al, using a reference point and critical compressibility")
    {
        auto trifluoroethane = Rackett(Rackett::CreateFromReferencePointB{346.3, 245, 75.38E-6, 0.255});

        REQUIRE(trifluoroethane(300.0) == Approx(91.4074E-6));
    }

    SECTION("Example 1 from Perry's")
    {
        auto acetonitrile = Rackett(Rackett::CreateFromCriticalProperties{545.5, 4.83E6, 0.184});

        REQUIRE(acetonitrile(376.69) == Approx(1.0 / 19.42E3).epsilon(0.001));
    }

    SECTION("Example 2 from Perry's")
    {
        auto acetonitrile = Rackett(Rackett::CreateFromCriticalProperties{545.5, 4.83E6, 0.202});

        // ===== This example uses Z_RA instead of the critical compressibility
        REQUIRE(acetonitrile(376.69) == Approx(1.0 / 16.577E3).epsilon(0.01));
    }

    SECTION("Coefficients from Perry's")
    {
        auto acetonitrile = PCProps::LiquidVolume::Rackett(1.0 / 1.3064E3, 0.22597, 545.5, 0.28678);

        REQUIRE(acetonitrile(229.32) == Approx(1.0 / 20.628E3).epsilon(0.001));
        REQUIRE(acetonitrile(545.5) == Approx(1.0 / 5.7813E3).epsilon(0.001));
    }

    SECTION("Coefficients from Perry's using copy constructor")
    {
        auto acetonitrile = PCProps::LiquidVolume::Rackett(1.0 / 1.3064E3, 0.22597, 545.5, 0.28678);
        auto copy(acetonitrile);

        REQUIRE(copy(229.32) == Approx(1.0 / 20.628E3).epsilon(0.001));
        REQUIRE(copy(545.5) == Approx(1.0 / 5.7813E3).epsilon(0.001));
    }

    SECTION("Coefficients from Perry's using move constructor")
    {
        auto acetonitrile = PCProps::LiquidVolume::Rackett(1.0 / 1.3064E3, 0.22597, 545.5, 0.28678);
        auto copy(std::move(acetonitrile));

        REQUIRE(copy(229.32) == Approx(1.0 / 20.628E3).epsilon(0.001));
        REQUIRE(copy(545.5) == Approx(1.0 / 5.7813E3).epsilon(0.001));
    }

    SECTION("Coefficients from Perry's using copy assignment")
    {
        auto acetonitrile = PCProps::LiquidVolume::Rackett(1.0 / 1.3064E3, 0.22597, 545.5, 0.28678);
        auto copy         = Rackett();
        copy              = acetonitrile;

        REQUIRE(copy(229.32) == Approx(1.0 / 20.628E3).epsilon(0.001));
        REQUIRE(copy(545.5) == Approx(1.0 / 5.7813E3).epsilon(0.001));
    }

    SECTION("Coefficients from Perry's using move assignment")
    {
        auto acetonitrile = PCProps::LiquidVolume::Rackett(1.0 / 1.3064E3, 0.22597, 545.5, 0.28678);
        auto copy         = Rackett();
        copy              = std::move(acetonitrile);

        REQUIRE(copy(229.32) == Approx(1.0 / 20.628E3).epsilon(0.001));
        REQUIRE(copy(545.5) == Approx(1.0 / 5.7813E3).epsilon(0.001));
    }
}