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

#include <LiquidVolume/CLVAalto.hpp>
#include <catch.hpp>

using PCProps::LiquidVolume::CLVAalto;

TEST_CASE("CLVAalto produces correct saturated liquid volume calculations")
{
    SECTION("Example 1 from Reid et. al 5th Edition @ 300K")
    {
        auto ammonia_psat = [](double _) { return 10.61E5; };
        auto ammonia_vsat = [](double _) { return 28.38E-6; };
        auto ammonia      = CLVAalto::create(405.4, 113.53E5, 0.256, ammonia_vsat, ammonia_psat);

        REQUIRE(ammonia(300.0, 400E5) == Approx(27.19E-6).epsilon(0.001));
    }

    SECTION("Example 1 from Reid et. al 5th Edition @ 400K")
    {
        auto ammonia_psat = [](double _) { return 102.97E5; };
        auto ammonia_vsat = [](double _) { return 49.15E-6; };
        auto ammonia      = CLVAalto::create(405.4, 113.53E5, 0.256, ammonia_vsat, ammonia_psat);

        REQUIRE(ammonia(400.0, 400E5) == Approx(35.60E-6).epsilon(0.001));
    }

    SECTION("Example 2 from Reid et. al 5th Edition")
    {
        auto ammonia_psat = [](double _) { return 1E5; };
        auto ammonia_vsat = [](double _) { return 127.31E-6; };
        auto ammonia      = CLVAalto::create(705.7, 45.6E5, 0.452, ammonia_vsat, ammonia_psat);

        REQUIRE(ammonia(503.15, 3000E5) == Approx(112.97E-6).epsilon(0.001));
    }

    SECTION("Copy constructor")
    {
        auto ammonia_psat = [](double _) { return 10.61E5; };
        auto ammonia_vsat = [](double _) { return 28.38E-6; };
        auto ammonia      = CLVAalto::create(405.4, 113.53E5, 0.256, ammonia_vsat, ammonia_psat);

        auto copy(ammonia);

        REQUIRE(copy(300.0, 400E5) == Approx(27.19E-6).epsilon(0.001));
    }

    SECTION("Move constructor")
    {
        auto ammonia_psat = [](double _) { return 10.61E5; };
        auto ammonia_vsat = [](double _) { return 28.38E-6; };
        auto ammonia      = CLVAalto::create(405.4, 113.53E5, 0.256, ammonia_vsat, ammonia_psat);

        auto copy(std::move(ammonia));

        REQUIRE(copy(300.0, 400E5) == Approx(27.19E-6).epsilon(0.001));
    }

    SECTION("Coefficients from Perry's using copy assignment")
    {
        auto ammonia_psat = [](double _) { return 10.61E5; };
        auto ammonia_vsat = [](double _) { return 28.38E-6; };
        auto ammonia      = CLVAalto::create(405.4, 113.53E5, 0.256, ammonia_vsat, ammonia_psat);
        auto copy         = CLVAalto::create(0.0, 0.0, 0.0, {}, {});
        copy              = ammonia;

        REQUIRE(copy(300.0, 400E5) == Approx(27.19E-6).epsilon(0.001));
    }

    SECTION("Coefficients from Perry's using move assignment")
    {
        auto ammonia_psat = [](double _) { return 10.61E5; };
        auto ammonia_vsat = [](double _) { return 28.38E-6; };
        auto ammonia      = CLVAalto::create(405.4, 113.53E5, 0.256, ammonia_vsat, ammonia_psat);
        auto copy         = CLVAalto::create(0.0, 0.0, 0.0, {}, {});
        copy              = std::move(ammonia);

        REQUIRE(copy(300.0, 400E5) == Approx(27.19E-6).epsilon(0.001));
    }
}