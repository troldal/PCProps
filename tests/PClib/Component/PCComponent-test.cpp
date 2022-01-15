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

#include "../../../source/PureComponent/PureComponent.hpp"
#include <HeatCapacity/AlyLee.hpp>
#include <catch/catch.hpp>

using namespace PCProps::HeatCapacity;

TEST_CASE("PCComponent Test")
{
    auto pcd = PCProps::PCComponentData {};

    pcd.name    = "ACETALDEHYDE";
    pcd.formula = "C2H4O";
    pcd.casrn   = "75-07-0";
    pcd.smiles  = "CC=O";

    pcd.molarWeight             = 44.053;
    pcd.boilingTemperature      = 293.300;
    pcd.freezingTemperature     = 0.0;
    pcd.criticalTemperature     = 466.00;
    pcd.criticalPressure        = 5550000.000;
    pcd.criticalVolume          = 0.000154;
    pcd.criticalDensity         = 286.05844;
    pcd.criticalCompressibility = 0.2206;
    pcd.acentricFactor          = 0.2907;

//    pcd.equationOfState = PCProps::EquationOfState::PengRobinson {};
    pcd.idealGasCpCorrelation = AlyLee(AlyLee::CreateFromDIPPR { 0.5192E5, 1.9245E5, 1.6265E3, 1.168E5, 723.6 });
    pcd.vaporPressureCorrelation          = [](double temperature) { return temperature; };
    pcd.satLiquidVolumeCorrelation                    = [](double temperature) { return temperature; };
    pcd.surfaceTensionCorrelation         = [](double temperature) { return temperature; };
    pcd.heatOfVaporizationCorrelation     = [](double temperature) { return temperature; };
    pcd.satVaporThermalConductivityCorrelation        = [](double temperature) { return temperature; };
    pcd.satLiquidThermalConductivityCorrelation       = [](double temperature) { return temperature; };
    pcd.satVaporViscosityCorrelation                  = [](double temperature) { return temperature; };
    pcd.satLiquidViscosityCorrelation                 = [](double temperature) { return temperature; };
    pcd.liquidCpCorrelation                  = [](double temperature) { return temperature; };

    auto pc = PCProps::PureComponent(pcd);

    /**
     * @brief Test that accessing values and functions of a default constructed PCComponent object will trow an exception,
     * and that the 'hasXXX' member functions will all yield false.
     */
    SECTION("Default Construction")
    {
        pcd = PCProps::PCComponentData {};
        pc  = PCProps::PureComponent(pcd);

        REQUIRE(pc.name().empty());
        REQUIRE(pc.formula().empty());
        REQUIRE(pc.casrn().empty());
        REQUIRE(pc.smiles().empty());

//        REQUIRE(pc.hasMolecularWeight() == false);
//        REQUIRE(pc.hasBoilingTemperature() == false);
//        REQUIRE(pc.hasFreezingTemperature() == false);
//        REQUIRE(pc.hasCriticalTemperature() == false);
//        REQUIRE(pc.hasCriticalPressure() == false);
//        REQUIRE(pc.hasCriticalVolume() == false);
//        REQUIRE(pc.hasCriticalDensity() == false);
//        REQUIRE(pc.hasCriticalCompressibility() == false);
//        REQUIRE(pc.hasAcentricFactor() == false);
//
//        REQUIRE_THROWS(pc.molecularWeight());
//        REQUIRE_THROWS(pc.boilingTemperature());
//        REQUIRE_THROWS(pc.criticalTemperature());
//        REQUIRE_THROWS(pc.criticalPressure());
//        REQUIRE_THROWS(pc.criticalVolume());
//        REQUIRE_THROWS(pc.criticalDensity());
//        REQUIRE_THROWS(pc.criticalCompressibility());
//        REQUIRE_THROWS(pc.acentricFactor());
//
//        REQUIRE(pc.hasVaporPressureFunction() == false);
//        REQUIRE(pc.hasLiquidDensityFunction() == false);
//        REQUIRE(pc.hasSurfaceTensionFunction() == false);
//        REQUIRE(pc.hasHeatOfVaporizationFunction() == false);
//        REQUIRE(pc.hasVaporThermalConductivityFunction() == false);
//        REQUIRE(pc.hasLiquidThermalConductivityFunction() == false);
//        REQUIRE(pc.hasVaporViscosityFunction() == false);
//        REQUIRE(pc.hasLiquidViscosityFunction() == false);
//        REQUIRE(pc.hasIdealGasCpFunction() == false);
//        REQUIRE(pc.hasLiquidCpFunction() == false);
//
//        REQUIRE_THROWS(pc.vaporPressure(0.0));
//        REQUIRE_THROWS(pc.liquidDensity(0.0));
//        REQUIRE_THROWS(pc.surfaceTension(0.0));
//        REQUIRE_THROWS(pc.heatOfVaporization(0.0));
//        REQUIRE_THROWS(pc.vaporThermalConductivity(0.0));
//        REQUIRE_THROWS(pc.liquidThermalConductivity(0.0));
//        REQUIRE_THROWS(pc.vaporViscosity(0.0, 0));
//        REQUIRE_THROWS(pc.liquidViscosity(0.0, 0));
//        REQUIRE_THROWS(pc.idealGasCp(0.0));
//        REQUIRE_THROWS(pc.liquidCp(0.0));
    }

    /**
     * @brief Test that an object initiated with metadata return the correct results.
     */
    SECTION("Object Construction - Metadata")
    {
        REQUIRE(pc.name() == "ACETALDEHYDE");
        REQUIRE(pc.formula() == "C2H4O");
        REQUIRE(pc.casrn() == "75-07-0");
        REQUIRE(pc.smiles() == "CC=O");
    }

    /**
     * @brief Test that a correctly initiated object returns the correct pure component constants.
     */
    SECTION("Object Construction - Constants")
    {
//        REQUIRE(pc.hasMolecularWeight());
//        REQUIRE(pc.hasBoilingTemperature());
//        REQUIRE(pc.hasFreezingTemperature());
//        REQUIRE(pc.hasCriticalTemperature());
//        REQUIRE(pc.hasCriticalPressure());
//        REQUIRE(pc.hasCriticalVolume());
//        REQUIRE(pc.hasCriticalDensity());
//        REQUIRE(pc.hasCriticalCompressibility());
//        REQUIRE(pc.hasAcentricFactor());
//
//        REQUIRE(pc.molecularWeight() == Approx(44.053));
//        REQUIRE(pc.boilingTemperature() == Approx(293.300));
//        REQUIRE(pc.freezingTemperature() == Approx(0.0));
//        REQUIRE(pc.criticalTemperature() == Approx(466.00));
//        REQUIRE(pc.criticalPressure() == Approx(5550000.00));
//        REQUIRE(pc.criticalVolume() == Approx(0.000154));
//        REQUIRE(pc.criticalDensity() == Approx(286.05844));
//        REQUIRE(pc.criticalCompressibility() == Approx(0.2206));
//        REQUIRE(pc.acentricFactor() == Approx(0.2907));
    }

    /**
     * @brief Test that a correctly initiated object returns the correct values, when the temperature
     * dependent correlation functions are called.
     */
    SECTION("Object Construction - Functions")
    {
//        REQUIRE(pc.hasVaporPressureFunction());
//        REQUIRE(pc.hasLiquidDensityFunction());
//        REQUIRE(pc.hasSurfaceTensionFunction());
//        REQUIRE(pc.hasHeatOfVaporizationFunction());
//        REQUIRE(pc.hasVaporThermalConductivityFunction());
//        REQUIRE(pc.hasLiquidThermalConductivityFunction());
//        REQUIRE(pc.hasVaporViscosityFunction());
//        REQUIRE(pc.hasLiquidViscosityFunction());
//        REQUIRE(pc.hasIdealGasCpFunction());
//        REQUIRE(pc.hasLiquidCpFunction());
//
//        REQUIRE(pc.vaporPressure(1.1) == Approx(1.1));
//        REQUIRE(pc.liquidDensity(2.2) == Approx(2.2));
//        REQUIRE(pc.surfaceTension(3.3) == Approx(3.3));
//        REQUIRE(pc.heatOfVaporization(4.4) == Approx(4.4));
//        REQUIRE(pc.vaporThermalConductivity(5.5) == Approx(5.5));
//        REQUIRE(pc.liquidThermalConductivity(6.6) == Approx(6.6));
//        REQUIRE(pc.vaporViscosity(7.7, 0) == Approx(7.7));
//        REQUIRE(pc.liquidViscosity(8.8, 0) == Approx(8.8));
//        REQUIRE(pc.idealGasCp(9.9) == Approx(9.9));
//        REQUIRE(pc.liquidCp(10.01) == Approx(10.01));
    }
}