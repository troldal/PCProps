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

#include <VaporPressure/AmbroseWalton.hpp>
#include <catch.hpp>

TEST_CASE("VPAmbroseWalton Test")
{
    auto psat = PCProps::VaporPressure::AmbroseWalton { 632.35, 45.1911E5, 0.249857 };

    SECTION("Default Construction")
    {
        psat = PCProps::VaporPressure::AmbroseWalton {};

        REQUIRE(psat.criticalTemperature() == Approx(0.0));
        REQUIRE(psat.criticalPressure() == Approx(0.0));
        REQUIRE(psat.acentricFactor() == Approx(0.0));

        REQUIRE(psat(300.0) != psat(300.0));
    }

    SECTION("Object Construction")
    {
        REQUIRE(psat.criticalTemperature() == Approx(632.35));
        REQUIRE(psat.criticalPressure() == Approx(45.1911E5));
        REQUIRE(psat.acentricFactor() == Approx(0.249857));

        REQUIRE(psat(300.00) == Approx(0.0178E5).epsilon(0.001));
        REQUIRE(psat(350.00) == Approx(0.174E5).epsilon(0.001));
        REQUIRE(psat(400.00) == Approx(0.885E5).epsilon(0.001));
        REQUIRE(psat(450.00) == Approx(2.98E5).epsilon(0.001));
        REQUIRE(psat(500.00) == Approx(7.67E5).epsilon(0.001));
        REQUIRE(psat(550.00) == Approx(16.43E5).epsilon(0.001));
        REQUIRE(psat(600.00) == Approx(31.14E5).epsilon(0.001));

    }

    SECTION("Object Copy Construction and Assignment")
    {
        auto psat2 {psat};

        REQUIRE(psat2(300.00) == Approx(0.0178E5).epsilon(0.001));
        REQUIRE(psat2(350.00) == Approx(0.174E5).epsilon(0.001));
        REQUIRE(psat2(400.00) == Approx(0.885E5).epsilon(0.001));
        REQUIRE(psat2(450.00) == Approx(2.98E5).epsilon(0.001));
        REQUIRE(psat2(500.00) == Approx(7.67E5).epsilon(0.001));
        REQUIRE(psat2(550.00) == Approx(16.43E5).epsilon(0.001));
        REQUIRE(psat2(600.00) == Approx(31.14E5).epsilon(0.001));

        auto psat3 = PCProps::VaporPressure::AmbroseWalton {};
        psat3 = psat2;

        REQUIRE(psat3(300.00) == Approx(0.0178E5).epsilon(0.001));
        REQUIRE(psat3(350.00) == Approx(0.174E5).epsilon(0.001));
        REQUIRE(psat3(400.00) == Approx(0.885E5).epsilon(0.001));
        REQUIRE(psat3(450.00) == Approx(2.98E5).epsilon(0.001));
        REQUIRE(psat3(500.00) == Approx(7.67E5).epsilon(0.001));
        REQUIRE(psat3(550.00) == Approx(16.43E5).epsilon(0.001));
        REQUIRE(psat3(600.00) == Approx(31.14E5).epsilon(0.001));
    }

    SECTION("Object Move Construction and Assignment")
    {
        auto psat2 {std::move(psat)};

        REQUIRE(psat2(300.00) == Approx(0.0178E5).epsilon(0.001));
        REQUIRE(psat2(350.00) == Approx(0.174E5).epsilon(0.001));
        REQUIRE(psat2(400.00) == Approx(0.885E5).epsilon(0.001));
        REQUIRE(psat2(450.00) == Approx(2.98E5).epsilon(0.001));
        REQUIRE(psat2(500.00) == Approx(7.67E5).epsilon(0.001));
        REQUIRE(psat2(550.00) == Approx(16.43E5).epsilon(0.001));
        REQUIRE(psat2(600.00) == Approx(31.14E5).epsilon(0.001));

        auto psat3 = PCProps::VaporPressure::AmbroseWalton {};
        psat3 = std::move(psat2);

        REQUIRE(psat3(300.00) == Approx(0.0178E5).epsilon(0.001));
        REQUIRE(psat3(350.00) == Approx(0.174E5).epsilon(0.001));
        REQUIRE(psat3(400.00) == Approx(0.885E5).epsilon(0.001));
        REQUIRE(psat3(450.00) == Approx(2.98E5).epsilon(0.001));
        REQUIRE(psat3(500.00) == Approx(7.67E5).epsilon(0.001));
        REQUIRE(psat3(550.00) == Approx(16.43E5).epsilon(0.001));
        REQUIRE(psat3(600.00) == Approx(31.14E5).epsilon(0.001));
    }

}