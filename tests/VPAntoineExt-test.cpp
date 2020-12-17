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

#include <VPAntoineExt.hpp>
#include <catch.hpp>

TEST_CASE("VPAntoineExt")
{
    auto psat = PCProps::VaporPressure::VPAntoineExt {1.9369E+02, -8.0367E+03, 0.00E+00, 0.00E+00, -2.9502E+01, 4.3678E-02, 1.0};

    SECTION("Default Construction")
    {
        psat = PCProps::VaporPressure::VPAntoineExt {};
        auto coeffs = psat.coefficients();

        REQUIRE(coeffs[0] == Approx(0.0));
        REQUIRE(coeffs[1] == Approx(0.0));
        REQUIRE(coeffs[2] == Approx(0.0));
        REQUIRE(coeffs[3] == Approx(0.0));
        REQUIRE(coeffs[4] == Approx(0.0));
        REQUIRE(coeffs[5] == Approx(0.0));
        REQUIRE(coeffs[6] == Approx(0.0));

        REQUIRE(psat(300.0) == Approx(1.0));
    }

    SECTION("Object Construction")
    {
        auto coeffs = psat.coefficients();

        REQUIRE(coeffs[0] == Approx(1.9369E+02));
        REQUIRE(coeffs[1] == Approx(-8.0367E+03));
        REQUIRE(coeffs[2] == Approx(0.0));
        REQUIRE(coeffs[3] == Approx(0.0));
        REQUIRE(coeffs[4] == Approx(-2.9502E+01));
        REQUIRE(coeffs[5] == Approx(4.3678E-02));
        REQUIRE(coeffs[6] == Approx(1.0));

        REQUIRE(psat(150.15) == Approx(3.23E-1).epsilon(0.001));
        REQUIRE(psat(466.00) == Approx(5.565E6).epsilon(0.001));
    }

    SECTION("Object Copy Construction and Assignment")
    {
        auto psat2 {psat};

        REQUIRE(psat2(150.15) == Approx(3.23E-1).epsilon(0.001));
        REQUIRE(psat2(466.00) == Approx(5.565E6).epsilon(0.001));

        auto psat3 = PCProps::VaporPressure::VPAntoineExt {};
        psat3 = psat2;

        REQUIRE(psat3(150.15) == Approx(3.23E-1).epsilon(0.001));
        REQUIRE(psat3(466.00) == Approx(5.565E6).epsilon(0.001));
    }

    SECTION("Object Move Construction and Assignment")
    {
        auto psat2 {std::move(psat)};

        REQUIRE(psat2(150.15) == Approx(3.23E-1).epsilon(0.001));
        REQUIRE(psat2(466.00) == Approx(5.565E6).epsilon(0.001));

        auto psat3 = PCProps::VaporPressure::VPAntoineExt {};
        psat3 = std::move(psat2);

        REQUIRE(psat3(150.15) == Approx(3.23E-1).epsilon(0.001));
        REQUIRE(psat3(466.00) == Approx(5.565E6).epsilon(0.001));
    }

}