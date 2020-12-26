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

#include <LiquidVolume/SLVElbro.hpp>
#include <catch.hpp>
#include <list>
#include <utility>

using PCProps::LiquidVolume::SLVElbro;

TEST_CASE("SLVElbro produces correct saturated liquid volume calculations")
{
    SECTION("Example 1 from Reid et. al 5th Edition, with std::vector")
    {
        auto hexadecane = SLVElbro::create({ { 1, 2 }, { 2, 14 } });

        REQUIRE(hexadecane(298.15) == Approx(294.39E-6).epsilon(0.001));
    }

    SECTION("Example 1 from Reid et. al 5th Edition, with std::list")
    {
        auto hexadecane = SLVElbro::create(std::list { std::make_pair(1, 2), std::make_pair(2, 14) });

        REQUIRE(hexadecane(298.15) == Approx(294.39E-6).epsilon(0.001));
    }

    SECTION("Example 2 from Reid et. al 5th Edition")
    {
        auto polymethylacrylate = SLVElbro::create({ { 1, 1 }, { 2, 1 }, { 22, 1 } });

        REQUIRE(polymethylacrylate(298.15) == Approx(71.4E-6).epsilon(0.001));
    }

    SECTION("Copy constructor")
    {
        auto hexadecane = SLVElbro::create({ { 1, 2 }, { 2, 14 } });
        auto copy(hexadecane);

        REQUIRE(copy(298.15) == Approx(294.39E-6).epsilon(0.001));
    }

    SECTION("Move constructor")
    {
        auto hexadecane = SLVElbro::create({ { 1, 2 }, { 2, 14 } });
        auto copy(std::move(hexadecane));

        REQUIRE(copy(298.15) == Approx(294.39E-6).epsilon(0.001));
    }

    SECTION("Coefficients from Perry's using copy assignment")
    {
        auto hexadecane = SLVElbro::create({ { 1, 2 }, { 2, 14 } });
        auto copy       = SLVElbro::create({});
        copy            = hexadecane;

        REQUIRE(copy(298.15) == Approx(294.39E-6).epsilon(0.001));
    }

    SECTION("Coefficients from Perry's using move assignment")
    {
        auto hexadecane = SLVElbro::create({ { 1, 2 }, { 2, 14 } });
        auto copy       = SLVElbro::create({});
        copy            = std::move(hexadecane);

        REQUIRE(copy(298.15) == Approx(294.39E-6).epsilon(0.001));
    }
}