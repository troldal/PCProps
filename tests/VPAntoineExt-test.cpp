//
// Created by Kenneth Balslev on 17/12/2020.
//

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