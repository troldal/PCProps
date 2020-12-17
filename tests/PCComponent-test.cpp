//
// Created by Kenneth Balslev on 17/12/2020.
//

#include <PCComponent.hpp>
#include <catch.hpp>

TEST_CASE("PCComponent")
{
    auto pcd = PCProps::PCComponentData {};

    pcd.name    = "ACETALDEHYDE";
    pcd.formula = "C2H4O";
    pcd.casrn   = "75-07-0";
    pcd.smiles  = "CC=O";

    pcd.molecularWeight         = 44.053;
    pcd.boilingTemperature      = 293.300;
    pcd.freezingTemperature     = 0.0;
    pcd.criticalTemperature     = 466.00;
    pcd.criticalPressure        = 5550000.000;
    pcd.criticalVolume          = 0.000154;
    pcd.criticalDensity         = 286.05844;
    pcd.criticalCompressibility = 0.2206;
    pcd.acentricFactor          = 0.2907;

    pcd.vaporPressureFunction             = [](double temperature) { return temperature; };
    pcd.liquidDensityFunction             = [](double temperature) { return temperature; };
    pcd.surfaceTensionFunction            = [](double temperature) { return temperature; };
    pcd.heatOfVaporizationFunction        = [](double temperature) { return temperature; };
    pcd.vaporThermalConductivityFunction  = [](double temperature) { return temperature; };
    pcd.liquidThermalConductivityFunction = [](double temperature) { return temperature; };
    pcd.vaporViscosityFunction            = [](double temperature) { return temperature; };
    pcd.liquidViscosityFunction           = [](double temperature) { return temperature; };
    pcd.idealGasCpFunction                = [](double temperature) { return temperature; };
    pcd.liquidCpFunction                  = [](double temperature) { return temperature; };

    auto pc = PCProps::PCComponent(pcd);

    /**
     * @brief Test that accessing values and functions of a default constructed PCComponent object will trow an exception,
     * and that the 'hasXXX' member functions will all yield false.
     */
    SECTION("Default Construction")
    {
        pcd = PCProps::PCComponentData {};
        pc  = PCProps::PCComponent(pcd);

        REQUIRE(pc.name().empty());
        REQUIRE(pc.formula().empty());
        REQUIRE(pc.casrn().empty());
        REQUIRE(pc.smiles().empty());

        REQUIRE(pc.hasMolecularWeight() == false);
        REQUIRE(pc.hasBoilingTemperature() == false);
        REQUIRE(pc.hasFreezingTemperature() == false);
        REQUIRE(pc.hasCriticalTemperature() == false);
        REQUIRE(pc.hasCriticalPressure() == false);
        REQUIRE(pc.hasCriticalVolume() == false);
        REQUIRE(pc.hasCriticalDensity() == false);
        REQUIRE(pc.hasCriticalCompressibility() == false);
        REQUIRE(pc.hasAcentricFactor() == false);

        REQUIRE_THROWS(pc.molecularWeight());
        REQUIRE_THROWS(pc.boilingTemperature());
        REQUIRE_THROWS(pc.criticalTemperature());
        REQUIRE_THROWS(pc.criticalPressure());
        REQUIRE_THROWS(pc.criticalVolume());
        REQUIRE_THROWS(pc.criticalDensity());
        REQUIRE_THROWS(pc.criticalCompressibility());
        REQUIRE_THROWS(pc.acentricFactor());

        REQUIRE(pc.hasVaporPressureFunction() == false);
        REQUIRE(pc.hasLiquidDensityFunction() == false);
        REQUIRE(pc.hasSurfaceTensionFunction() == false);
        REQUIRE(pc.hasHeatOfVaporizationFunction() == false);
        REQUIRE(pc.hasVaporThermalConductivityFunction() == false);
        REQUIRE(pc.hasLiquidThermalConductivityFunction() == false);
        REQUIRE(pc.hasVaporViscosityFunction() == false);
        REQUIRE(pc.hasLiquidViscosityFunction() == false);
        REQUIRE(pc.hasIdealGasCpFunction() == false);
        REQUIRE(pc.hasLiquidCpFunction() == false);

        REQUIRE_THROWS(pc.vaporPressure(0.0));
        REQUIRE_THROWS(pc.liquidDensity(0.0));
        REQUIRE_THROWS(pc.surfaceTension(0.0));
        REQUIRE_THROWS(pc.heatOfVaporization(0.0));
        REQUIRE_THROWS(pc.vaporThermalConductivity(0.0));
        REQUIRE_THROWS(pc.liquidThermalConductivity(0.0));
        REQUIRE_THROWS(pc.vaporViscosity(0.0));
        REQUIRE_THROWS(pc.liquidViscosity(0.0));
        REQUIRE_THROWS(pc.idealGasCp(0.0));
        REQUIRE_THROWS(pc.liquidCp(0.0));
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
        REQUIRE(pc.hasMolecularWeight());
        REQUIRE(pc.hasBoilingTemperature());
        REQUIRE(pc.hasFreezingTemperature());
        REQUIRE(pc.hasCriticalTemperature());
        REQUIRE(pc.hasCriticalPressure());
        REQUIRE(pc.hasCriticalVolume());
        REQUIRE(pc.hasCriticalDensity());
        REQUIRE(pc.hasCriticalCompressibility());
        REQUIRE(pc.hasAcentricFactor());

        REQUIRE(pc.molecularWeight() == Approx(44.053));
        REQUIRE(pc.boilingTemperature() == Approx(293.300));
        REQUIRE(pc.freezingTemperature() == Approx(0.0));
        REQUIRE(pc.criticalTemperature() == Approx(466.00));
        REQUIRE(pc.criticalPressure() == Approx(5550000.00));
        REQUIRE(pc.criticalVolume() == Approx(0.000154));
        REQUIRE(pc.criticalDensity() == Approx(286.05844));
        REQUIRE(pc.criticalCompressibility() == Approx(0.2206));
        REQUIRE(pc.acentricFactor() == Approx(0.2907));
    }

    /**
     * @brief Test that a correctly initiated object returns the correct values, when the temperature
     * dependent correlation functions are called.
     */
    SECTION("Object Construction - Functions")
    {
        REQUIRE(pc.hasVaporPressureFunction());
        REQUIRE(pc.hasLiquidDensityFunction());
        REQUIRE(pc.hasSurfaceTensionFunction());
        REQUIRE(pc.hasHeatOfVaporizationFunction());
        REQUIRE(pc.hasVaporThermalConductivityFunction());
        REQUIRE(pc.hasLiquidThermalConductivityFunction());
        REQUIRE(pc.hasVaporViscosityFunction());
        REQUIRE(pc.hasLiquidViscosityFunction());
        REQUIRE(pc.hasIdealGasCpFunction());
        REQUIRE(pc.hasLiquidCpFunction());

        REQUIRE(pc.vaporPressure(1.1) == Approx(1.1));
        REQUIRE(pc.liquidDensity(2.2) == Approx(2.2));
        REQUIRE(pc.surfaceTension(3.3) == Approx(3.3));
        REQUIRE(pc.heatOfVaporization(4.4) == Approx(4.4));
        REQUIRE(pc.vaporThermalConductivity(5.5) == Approx(5.5));
        REQUIRE(pc.liquidThermalConductivity(6.6) == Approx(6.6));
        REQUIRE(pc.vaporViscosity(7.7) == Approx(7.7));
        REQUIRE(pc.liquidViscosity(8.8) == Approx(8.8));
        REQUIRE(pc.idealGasCp(9.9) == Approx(9.9));
        REQUIRE(pc.liquidCp(10.01) == Approx(10.01));
    }
}