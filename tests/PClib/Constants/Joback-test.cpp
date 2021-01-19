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

#include <ConstantData/Joback.hpp>
#include <catch/catch.hpp>

using PCProps::ConstantData::JobackGroup;
using PCProps::ConstantData::Joback;

TEST_CASE("CDJoback Test")
{
    SECTION("Property Estimation - Wikipedia Example")
    {
        Joback acetone(std::vector<JobackGroup> { { 1, 2 }, { 24, 1 } }, 58.08, 10, 0.0);

        REQUIRE(acetone.boilingTemperatureIsValid());
        REQUIRE(acetone.meltingTemperatureIsValid());
        REQUIRE(acetone.criticalTemperatureIsValid());
        REQUIRE(acetone.criticalPressureIsValid());
        REQUIRE(acetone.criticalVolumeIsValid());
        REQUIRE(acetone.enthalpyOfFormationIsValid());
        REQUIRE(acetone.gibbsEnergyOfFormationIsValid());
        REQUIRE(acetone.enthalpyOfFusionIsValid());
        REQUIRE(acetone.enthalpyOfVaporizationIsValid());
        REQUIRE(acetone.idealGasCpIsValid());
        REQUIRE(acetone.liquidViscosityIsValid());

        REQUIRE(acetone.boilingTemperature() == Approx(322.11).epsilon(0.001));
        REQUIRE(acetone.meltingTemperature() == Approx(173.5).epsilon(0.001));
        REQUIRE(acetone.criticalTemperature() == Approx(500.5590).epsilon(0.001));
        REQUIRE(acetone.criticalPressure() == Approx(48.0250E5).epsilon(0.001));
        REQUIRE(acetone.criticalVolume() == Approx(209.5000E-6).epsilon(0.001));
        REQUIRE(acetone.enthalpyOfFormation() == Approx(-217.8300E3).epsilon(0.001));
        REQUIRE(acetone.gibbsEnergyOfFormation() == Approx(-154.5400E3).epsilon(0.001));
        REQUIRE(acetone.enthalpyOfFusion() == Approx(5.1250E3).epsilon(0.001));
        REQUIRE(acetone.enthalpyOfVaporization() == Approx(29.0180E3).epsilon(0.001));
        REQUIRE(acetone.idealGasCp(300.0) == Approx(75.3264).epsilon(0.001));
        REQUIRE(acetone.liquidViscosity(300.0) == Approx(0.0002942).epsilon(0.001));
    }

    SECTION("Property Estimation - Perry Example")
    {
        Joback o_xylene(std::vector<JobackGroup> { { 14, 4 }, { 15, 2 }, { 1, 2 } }, 122.16, 18, 417.58);

        REQUIRE(o_xylene.criticalTemperatureIsValid());
        REQUIRE(o_xylene.criticalPressureIsValid());
        REQUIRE(o_xylene.criticalVolumeIsValid());

        REQUIRE(o_xylene.criticalTemperature() == Approx(630.37).epsilon(0.001));
        REQUIRE(o_xylene.criticalPressure() == Approx(35.86E5).epsilon(0.001));
        REQUIRE(o_xylene.criticalVolume() == Approx(375.5E-6).epsilon(0.001));
    }

    SECTION("Property Estimation - Poling & Prausnitz Example")
    {
        Joback ethylphenol_a(std::vector<JobackGroup> { { 1, 1 }, { 2, 1 }, { 14, 4 }, { 15, 2 }, { 21, 1 } }, 106.16, 19, 477.67);

        REQUIRE(ethylphenol_a.criticalTemperatureIsValid());
        REQUIRE(ethylphenol_a.criticalPressureIsValid());
        REQUIRE(ethylphenol_a.criticalVolumeIsValid());

        REQUIRE(ethylphenol_a.criticalTemperature() == Approx(698.1).epsilon(0.001));
        REQUIRE(ethylphenol_a.criticalPressure() == Approx(44.09E5).epsilon(0.001));
        REQUIRE(ethylphenol_a.criticalVolume() == Approx(341.5E-6).epsilon(0.001));

        Joback ethylphenol_b(std::vector<JobackGroup> { { 1, 1 }, { 2, 1 }, { 14, 4 }, { 15, 2 }, { 21, 1 } }, 106.16, 19);

        REQUIRE(ethylphenol_b.criticalTemperatureIsValid());
        REQUIRE(ethylphenol_b.criticalTemperature() == Approx(715.5).epsilon(0.001));
    }
}