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

#include <LiquidVolume/HankinsonThomson.hpp>
#include <catch.hpp>

using PCProps::LiquidVolume::HankinsonThomson;

TEST_CASE("SLVHankinsonThomson produces correct saturated liquid volume calculations")
{
    SECTION("Example from Reid et. al 4th Edition")
    {
        auto isobutane = PCProps::LiquidVolume::HankinsonThomson(408.04, 256.8E-6, 0.1825);

        REQUIRE(isobutane(310.93) == Approx(108.9E-6).epsilon(0.001));
    }

    SECTION("Example using estimated properties")
    {
        auto isobutane2 = HankinsonThomson::createFromEstimatedProperties(408.04, 3640000, 0.1825);

        REQUIRE(isobutane2(310.93) == Approx(108.255E-6).epsilon(0.001));
    }

    SECTION("Example from Gmehling et. al 2nd Edition")
    {
        auto n_hexane = PCProps::LiquidVolume::HankinsonThomson(507.8, 386.8E-6, 0.3002);

        REQUIRE(n_hexane(293.15) == Approx(136.85E-6).epsilon(0.001));
    }

    SECTION("Copy constructor")
    {
        auto n_hexane = PCProps::LiquidVolume::HankinsonThomson(507.8, 386.8E-6, 0.3002);
        auto copy(n_hexane);

        REQUIRE(copy(293.15) == Approx(136.85E-6).epsilon(0.001));
    }

    SECTION("Move constructor")
    {
        auto n_hexane = PCProps::LiquidVolume::HankinsonThomson(507.8, 386.8E-6, 0.3002);
        auto copy(std::move(n_hexane));

        REQUIRE(copy(293.15) == Approx(136.85E-6).epsilon(0.001));
    }

    SECTION("Coefficients from Perry's using copy assignment")
    {
        auto n_hexane = PCProps::LiquidVolume::HankinsonThomson(507.8, 386.8E-6, 0.3002);
        auto copy     = PCProps::LiquidVolume::HankinsonThomson(0.0, 0.0, 0.0);
        copy          = n_hexane;

        REQUIRE(copy(293.15) == Approx(136.85E-6).epsilon(0.001));
    }

    SECTION("Coefficients from Perry's using move assignment")
    {
        auto n_hexane = PCProps::LiquidVolume::HankinsonThomson(507.8, 386.8E-6, 0.3002);
        auto copy     = PCProps::LiquidVolume::HankinsonThomson(0.0, 0.0, 0.0);
        copy          = std::move(n_hexane);

        REQUIRE(copy(293.15) == Approx(136.85E-6).epsilon(0.001));
    }
}