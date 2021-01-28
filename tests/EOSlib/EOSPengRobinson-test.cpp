//
// Created by Kenneth Balslev on 03/01/2021.
//

#include <HeatCapacity/AlyLee.hpp>
#include <PengRobinson/PengRobinson.hpp>
#include <catch/catch.hpp>
#include <common/PropertyData.hpp>

using PCProps::EquationOfState::PengRobinson;
using PCProps::PCPhase;
using PCProps::EquationOfState::PengRobinson;
using PCProps::PCPhases;

using namespace PCProps;
using namespace PCProps::EquationOfState;
using namespace PCProps::HeatCapacity;

TEST_CASE("PengRobinson Test")
{
    auto tc    = 369.83;
    auto pc    = 4.248E6;
    auto omega = 0.1523;
    auto propane = PengRobinson(tc, pc, omega);

    auto cp = AlyLee(AlyLee::CreateFromDIPPR { 0.5192E5, 1.9245E5, 1.6265E3, 1.168E5, 723.6 });
    propane.setIdealGasCpFunction(cp);

    SECTION("PT Flash of propane @ 273.15 K and 1 bar")
    {
        auto result = propane.flashPT(100000.0, 273.15)[0];

        REQUIRE(result[PCTemperature] == Approx(273.15).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.978939).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.97931788809).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.022233).epsilon(0.001));
        REQUIRE(result[PCEnthalpy] == Approx(-1899.261857).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(-6.496219).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-4122.525469).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(-124.819694).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(-2348.083305).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 223.15 K and 1 bar - Compressed Liquid")
    {
        auto result = propane.flashPT(100000.0, 223.15)[0];

        REQUIRE(result[PCTemperature] == Approx(223.15).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.003766).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.69350105382).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.000070).epsilon(0.01));
        REQUIRE(result[PCEnthalpy] == Approx(-24084.251701).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(-101.724616).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-24091.238463).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(-1384.403549).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(-1391.390312).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 323.15 K and 1 bar - Superheated Vapor")
    {
        auto result = propane.flashPT(100000.0, 323.15)[0];

        REQUIRE(result[PCTemperature] == Approx(323.15).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.986760).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.98690029547).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.026512).epsilon(0.001));
        REQUIRE(result[PCEnthalpy] == Approx(1803.116466).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(5.926803).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-848.128978).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(-112.130043).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(-2763.375487).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 469.8 K and 1 bar - T > Tc")
    {
        auto result = propane.flashPT(100000.0, 469.8)[0];

        REQUIRE(result[PCTemperature] == Approx(469.8).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.995756).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.99576428359).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.038896).epsilon(0.001));
        REQUIRE(result[PCEnthalpy] == Approx(15577.940459).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(40.740519).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(11688.385381).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(-3561.955560).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(-7451.510638).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 323.15 K and 50 bar - P > Pc")
    {
        auto result = propane.flashPT(5000000.0, 323.15)[0];

        REQUIRE(result[PCTemperature] == Approx(323.15).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(5000000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.171348).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.3014710365).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.000092).epsilon(0.001));
        REQUIRE(result[PCEnthalpy] == Approx(-13089.953760).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(-62.826659).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-13550.335573).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(7212.481187).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(6752.099374).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 469.8 K and 50 bar - Supercritical")
    {
        auto result = propane.flashPT(5000000.0, 469.8)[0];

        REQUIRE(result[PCTemperature] == Approx(469.8).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(5000000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.794804).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.8106014669).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.000621).epsilon(0.001));
        REQUIRE(result[PCEnthalpy] == Approx(12596.727773).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(3.579012).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(9492.115417).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(10915.307766).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(7810.695410).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 10 K and 10 Pa - Low P & T")
    {
        auto result = propane.flashPT(10, 10)[0];

        REQUIRE(result[PCTemperature] == Approx(10.0).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(10.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.000007).epsilon(0.1));
//        REQUIRE(result[PCFugacityCoefficient] == Approx(0.0).epsilon(0.01));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.000056).epsilon(0.01));
        REQUIRE(result[PCEnthalpy] == Approx(-43434.5901840318).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(-417.1895727118).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-43434.5907494721).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(-39262.6944569137).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(-39262.695022354).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ Tc and Pc - @ Critical Point")
    {
        auto result = propane.flashPT(pc, tc)[0];

        REQUIRE(result[PCTemperature] == Approx(369.83).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(4248000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.305823).epsilon(0.1));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.6426352828).epsilon(0.01));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.000221).epsilon(0.01));
        REQUIRE(result[PCEnthalpy] == Approx(-2372.295636).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(-32.215245).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-3312.681919).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(9541.868409).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(8601.482126).epsilon(0.001));
    }

    SECTION("Px Flash of propane @ 1 bar and x = 0.0 - Saturated Liquid")
    {
        auto result = propane.flashPx(100000.0, 0.0)[0];

        REQUIRE(result[PCTemperature] == Approx(230.657526).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.003697).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.96789922051).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.000071).epsilon(0.01));
        REQUIRE(result[PCEnthalpy] == Approx(-23379.777150).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(-98.619779).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-23386.868147).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(-632.382962).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(-639.473959).epsilon(0.001));
    }

    SECTION("Px Flash of propane @ 1 bar and x = 1.0 - Saturated Vapor")
    {
        auto result = propane.flashPx(100000.0, 1.0)[0];

        REQUIRE(result[PCTemperature] == Approx(230.657526).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.966927).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.96789922052).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.018544).epsilon(0.001));
        REQUIRE(result[PCEnthalpy] == Approx(-4667.418728).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(-17.493623).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-6521.785497).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(-632.382962).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(-2486.749730).epsilon(0.001));
    }

    SECTION("Px Flash of propane @ 1 bar and x = 0.5 - Two-Phase")
    {
        auto result = propane.flashPx(100000.0, 0.5);

        REQUIRE(result[1][PCTemperature] == Approx(230.657526).epsilon(0.001));
        REQUIRE(result[1][PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[1][PCCompressibility] == Approx(0.966927).epsilon(0.001));
        REQUIRE(result[1][PCFugacityCoefficient] == Approx(0.96789922052).epsilon(0.001));
        REQUIRE(result[1][PCMolarFlow] == Approx(0.5).epsilon(0.001));
        REQUIRE(result[1][PCMolarVolume] == Approx(0.018544).epsilon(0.001));
        REQUIRE(result[1][PCEnthalpy] == Approx(-4667.418728).epsilon(0.001));
        REQUIRE(result[1][PCEntropy] == Approx(-17.493623).epsilon(0.001));
        REQUIRE(result[1][PCInternalEnergy] == Approx(-6521.785497).epsilon(0.001));
        REQUIRE(result[1][PCGibbsEnergy] == Approx(-632.382962).epsilon(0.001));
        REQUIRE(result[1][PCHelmholzEnergy] == Approx(-2486.749730).epsilon(0.001));

        REQUIRE(result[0][PCTemperature] == Approx(230.657526).epsilon(0.001));
        REQUIRE(result[0][PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[0][PCCompressibility] == Approx(0.003697).epsilon(0.001));
        REQUIRE(result[0][PCFugacityCoefficient] == Approx(0.96789922051).epsilon(0.001));
        REQUIRE(result[0][PCMolarFlow] == Approx(0.5).epsilon(0.001));
        REQUIRE(result[0][PCMolarVolume] == Approx(0.000071).epsilon(0.01));
        REQUIRE(result[0][PCEnthalpy] == Approx(-23379.777150).epsilon(0.001));
        REQUIRE(result[0][PCEntropy] == Approx(-98.619779).epsilon(0.001));
        REQUIRE(result[0][PCInternalEnergy] == Approx(-23386.868147).epsilon(0.001));
        REQUIRE(result[0][PCGibbsEnergy] == Approx(-632.382962).epsilon(0.001));
        REQUIRE(result[0][PCHelmholzEnergy] == Approx(-639.473959).epsilon(0.001));
    }

    SECTION("Tx Flash of propane @ 223.15 K and x = 0.0 - Saturated Liquid")
    {
        auto result = propane.flashTx(223.15, 0.0)[0];

        REQUIRE(result[PCTemperature] == Approx(223.15).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(71038.119363).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.002675).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.9751737591).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.000070).epsilon(0.01));
        REQUIRE(result[PCEnthalpy] == Approx(-24085.412742).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(-101.720751).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-24090.376323).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(-1386.427112).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(-1391.390693).epsilon(0.001));
    }

    SECTION("Tx Flash of propane @ 223.15 K and x = 1.0 - Saturated Vapor")
    {
        auto result = propane.flashTx(223.15, 1.0)[0];

        REQUIRE(result[PCTemperature] == Approx(223.15).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(71038.119363).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.974596).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.9751737591).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.025454).epsilon(0.001));
        REQUIRE(result[PCEnthalpy] == Approx(-5075.219460).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(-16.530551).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-6883.457001).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(-1386.427112).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(-3194.664652).epsilon(0.001));
    }

    SECTION("Tx Flash of propane @ 223.15 K and x = 0.5 - Two-Phase")
    {
        auto result = propane.flashTx(223.15, 0.5);

        REQUIRE(result[1][PCTemperature] == Approx(223.15).epsilon(0.001));
        REQUIRE(result[1][PCPressure] == Approx(71038.119363).epsilon(0.001));
        REQUIRE(result[1][PCCompressibility] == Approx(0.974596).epsilon(0.001));
        REQUIRE(result[1][PCFugacityCoefficient] == Approx(0.9751737591).epsilon(0.001));
        REQUIRE(result[1][PCMolarFlow] == Approx(0.5).epsilon(0.001));
        REQUIRE(result[1][PCMolarVolume] == Approx(0.025454).epsilon(0.001));
        REQUIRE(result[1][PCEnthalpy] == Approx(-5075.219460).epsilon(0.001));
        REQUIRE(result[1][PCEntropy] == Approx(-16.530551).epsilon(0.001));
        REQUIRE(result[1][PCInternalEnergy] == Approx(-6883.457001).epsilon(0.001));
        REQUIRE(result[1][PCGibbsEnergy] == Approx(-1386.427112).epsilon(0.001));
        REQUIRE(result[1][PCHelmholzEnergy] == Approx(-3194.664652).epsilon(0.001));

        REQUIRE(result[0][PCTemperature] == Approx(223.15).epsilon(0.001));
        REQUIRE(result[0][PCPressure] == Approx(71038.119363).epsilon(0.001));
        REQUIRE(result[0][PCCompressibility] == Approx(0.002675).epsilon(0.001));
        REQUIRE(result[0][PCFugacityCoefficient] == Approx(0.9751737591).epsilon(0.001));
        REQUIRE(result[0][PCMolarFlow] == Approx(0.5).epsilon(0.001));
        REQUIRE(result[0][PCMolarVolume] == Approx(0.000070).epsilon(0.01));
        REQUIRE(result[0][PCEnthalpy] == Approx(-24085.412742).epsilon(0.001));
        REQUIRE(result[0][PCEntropy] == Approx(-101.720751).epsilon(0.001));
        REQUIRE(result[0][PCInternalEnergy] == Approx(-24090.376323).epsilon(0.001));
        REQUIRE(result[0][PCGibbsEnergy] == Approx(-1386.427112).epsilon(0.001));
        REQUIRE(result[0][PCHelmholzEnergy] == Approx(-1391.390693).epsilon(0.001));
    }

    SECTION("PH Flash of propane @ 1 bar and H = -14000 - Two-Phase")
    {
        auto result = propane.flashPH(100000, -14000);

        REQUIRE(result[1][PCTemperature] == Approx(230.657526).epsilon(0.001));
        REQUIRE(result[1][PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[1][PCCompressibility] == Approx(0.966927).epsilon(0.001));
        REQUIRE(result[1][PCFugacityCoefficient] == Approx(0.96789922052).epsilon(0.001));
        REQUIRE(result[1][PCMolarFlow] == Approx(0.501261).epsilon(0.001));
        REQUIRE(result[1][PCMolarVolume] == Approx(0.018544).epsilon(0.001));
        REQUIRE(result[1][PCEnthalpy] == Approx(-4667.418728).epsilon(0.001));
        REQUIRE(result[1][PCEntropy] == Approx(-17.493623).epsilon(0.001));
        REQUIRE(result[1][PCInternalEnergy] == Approx(-6521.785497).epsilon(0.001));
        REQUIRE(result[1][PCGibbsEnergy] == Approx(-632.382962).epsilon(0.001));
        REQUIRE(result[1][PCHelmholzEnergy] == Approx(-2486.749730).epsilon(0.001));

        REQUIRE(result[0][PCTemperature] == Approx(230.657526).epsilon(0.001));
        REQUIRE(result[0][PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[0][PCCompressibility] == Approx(0.003697).epsilon(0.001));
        REQUIRE(result[0][PCFugacityCoefficient] == Approx(0.96789922051).epsilon(0.001));
        REQUIRE(result[0][PCMolarFlow] == Approx(0.498739).epsilon(0.001));
        REQUIRE(result[0][PCMolarVolume] == Approx(0.000071).epsilon(0.01));
        REQUIRE(result[0][PCEnthalpy] == Approx(-23379.777150).epsilon(0.001));
        REQUIRE(result[0][PCEntropy] == Approx(-98.619779).epsilon(0.001));
        REQUIRE(result[0][PCInternalEnergy] == Approx(-23386.868147).epsilon(0.001));
        REQUIRE(result[0][PCGibbsEnergy] == Approx(-632.382962).epsilon(0.001));
        REQUIRE(result[0][PCHelmholzEnergy] == Approx(-639.473959).epsilon(0.001));
    }

    SECTION("PH Flash of propane @ 1 bar and H = 2000 - Superheated Vapor")
    {
        auto result = propane.flashPH(100000, 2000)[0];

        REQUIRE(result[PCTemperature] == Approx(325.627654).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.987043).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.98717657406).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.026723).epsilon(0.001));
        REQUIRE(result[PCEnthalpy] == Approx(2000).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(6.533741).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-672.338547).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(-127.566650).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(-2799.905197).epsilon(0.001));
    }

    SECTION("PH Flash of propane @ 1 bar and H = 16000 - T > Tc")
    {
        auto result = propane.flashPH(100000, 16000)[0];

        REQUIRE(result[PCTemperature] == Approx(473.704635).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.995872).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.99587960365).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.039223).epsilon(0.001));
        REQUIRE(result[PCEnthalpy] == Approx(16000).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(41.635185).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(12077.660810).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(-3722.780350).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(-7645.119541).epsilon(0.001));
    }

    SECTION("PH Flash of propane @ 1 bar and H = -24000 - Compressed Liquid")
    {
        auto result = propane.flashPH(100000, -24000)[0];

        REQUIRE(result[PCTemperature] == Approx(224.055958).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.003757).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.72299863177).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.000070).epsilon(0.01));
        REQUIRE(result[PCEnthalpy] == Approx(-24000).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(-101.347824).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-24006.998906).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(-1292.416115).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(-1299.415022).epsilon(0.001));
    }

    SECTION("PH Flash of propane @ 50 bar and H = -20000 - P > Pc")
    {
        auto result = propane.flashPH(5000000, -20000)[0];

        REQUIRE(result[PCTemperature] == Approx(263.308763).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(5000000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.171527).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.0749889902).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.000075).epsilon(0.01));
        REQUIRE(result[PCEnthalpy] == Approx(-20000).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(-86.373130).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-20375.518173).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(2742.802087).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(2367.283914).epsilon(0.001));
    }

    SECTION("PH Flash of propane @ 50 bar and H = 12000 - Supercritical")
    {
        auto result = propane.flashPH(5000000, 12000)[0];

        REQUIRE(result[PCTemperature] == Approx(464.976124).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(5000000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.785893).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.8040055372).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.000608).epsilon(0.01));
        REQUIRE(result[PCEnthalpy] == Approx(12000).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(2.302270).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(8961.717250).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(10929.499302).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(7891.216551).epsilon(0.001));
    }

    SECTION("PS Flash of propane @ 1 bar and S = -60 - Two-Phase")
    {
        auto result = propane.flashPS(100000, -60);

        REQUIRE(result[1][PCTemperature] == Approx(230.657526).epsilon(0.001));
        REQUIRE(result[1][PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[1][PCCompressibility] == Approx(0.966927).epsilon(0.001));
        REQUIRE(result[1][PCFugacityCoefficient] == Approx(0.96789922052).epsilon(0.001));
        REQUIRE(result[1][PCMolarFlow] == Approx(0.476046).epsilon(0.001));
        REQUIRE(result[1][PCMolarVolume] == Approx(0.018544).epsilon(0.001));
        REQUIRE(result[1][PCEnthalpy] == Approx(-4667.418728).epsilon(0.001));
        REQUIRE(result[1][PCEntropy] == Approx(-17.493623).epsilon(0.001));
        REQUIRE(result[1][PCInternalEnergy] == Approx(-6521.785497).epsilon(0.001));
        REQUIRE(result[1][PCGibbsEnergy] == Approx(-632.382962).epsilon(0.001));
        REQUIRE(result[1][PCHelmholzEnergy] == Approx(-2486.749730).epsilon(0.001));

        REQUIRE(result[0][PCTemperature] == Approx(230.657526).epsilon(0.001));
        REQUIRE(result[0][PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[0][PCCompressibility] == Approx(0.003697).epsilon(0.001));
        REQUIRE(result[0][PCFugacityCoefficient] == Approx(0.96789922051).epsilon(0.001));
        REQUIRE(result[0][PCMolarFlow] == Approx(0.523954).epsilon(0.001));
        REQUIRE(result[0][PCMolarVolume] == Approx(0.000071).epsilon(0.01));
        REQUIRE(result[0][PCEnthalpy] == Approx(-23379.777150).epsilon(0.001));
        REQUIRE(result[0][PCEntropy] == Approx(-98.619779).epsilon(0.001));
        REQUIRE(result[0][PCInternalEnergy] == Approx(-23386.868147).epsilon(0.001));
        REQUIRE(result[0][PCGibbsEnergy] == Approx(-632.382962).epsilon(0.001));
        REQUIRE(result[0][PCHelmholzEnergy] == Approx(-639.473959).epsilon(0.001));
    }

    SECTION("PS Flash of propane @ 1 bar and S = 7 - Superheated Vapor")
    {
        auto result = propane.flashPS(100000, 7)[0];

        REQUIRE(result[PCTemperature] == Approx(327.532918).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.987255).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.98738413477).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.026886).epsilon(0.001));
        REQUIRE(result[PCEnthalpy] == Approx(2152.271027).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(7.0).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-536.281691).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(-140.459399).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(-2829.012117).epsilon(0.001));
    }

    SECTION("PS Flash of propane @ 1 bar and S = 40 - T > Tc")
    {
        auto result = propane.flashPS(100000, 40)[0];

        REQUIRE(result[PCTemperature] == Approx(466.573938).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.995657).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.99566632742).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.038625).epsilon(0.001));
        REQUIRE(result[PCEnthalpy] == Approx(15231.239828).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(40.0).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(11368.776266).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(-3431.717744).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(-7294.181306).epsilon(0.001));
    }

    SECTION("PS Flash of propane @ 1 bar and S = -100 - Compressed Liquid")
    {
        auto result = propane.flashPS(100000, -100)[0];

        REQUIRE(result[PCTemperature] == Approx(227.308593).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.003727).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.83688539684).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.000070).epsilon(0.01));
        REQUIRE(result[PCEnthalpy] == Approx(-23695.822221).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(-100.0).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-23702.865692).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(-964.962920).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(-972.006390).epsilon(0.001));
    }

    SECTION("PS Flash of propane @ 50 bar and S = -85 - P > Pc")
    {
        auto result = propane.flashPS(5000000, -85)[0];

        REQUIRE(result[PCTemperature] == Approx(266.835388).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(5000000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.170734).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.0833455104).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.000076).epsilon(0.01));
        REQUIRE(result[PCEnthalpy] == Approx(-19636.022175).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(-85.0).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-20014.811861).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(3044.985805).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(2666.196118).epsilon(0.001));
    }

    SECTION("PS Flash of propane @ 50 bar and S = 3 - Supercritical")
    {
        auto result = propane.flashPS(5000000, 3.0)[0];

        REQUIRE(result[PCTemperature] == Approx(467.605973).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(5000000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.790814).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.8076329386).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.000615).epsilon(0.01));
        REQUIRE(result[PCEnthalpy] == Approx(12325.343803).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(3.0).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(9250.744824).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(10922.525896).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(7847.926917).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 273.15 K and 1 bar - Copy Constructor")
    {
        auto temp   = propane;
        auto result = temp.flashPT(100000.0, 273.15)[0];

        REQUIRE(result[PCTemperature] == Approx(273.15).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.978939).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.97931788809).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.022233).epsilon(0.001));
        REQUIRE(result[PCEnthalpy] == Approx(-1899.261857).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(-6.496219).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-4122.525469).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(-124.819694).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(-2348.083305).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 273.15 K and 1 bar - Move Constructor")
    {
        auto temp   = std::move(propane);
        auto result = temp.flashPT(100000.0, 273.15)[0];

        REQUIRE(result[PCTemperature] == Approx(273.15).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.978939).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.97931788809).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.022233).epsilon(0.001));
        REQUIRE(result[PCEnthalpy] == Approx(-1899.261857).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(-6.496219).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-4122.525469).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(-124.819694).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(-2348.083305).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 273.15 K and 1 bar - Copy Assignment")
    {
        PengRobinson temp {};
        temp        = propane;
        auto result = temp.flashPT(100000.0, 273.15)[0];

        REQUIRE(result[PCTemperature] == Approx(273.15).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.978939).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.97931788809).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.022233).epsilon(0.001));
        REQUIRE(result[PCEnthalpy] == Approx(-1899.261857).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(-6.496219).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-4122.525469).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(-124.819694).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(-2348.083305).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 273.15 K and 1 bar - Move Assignment")
    {
        PengRobinson temp {};
        temp        = std::move(propane);
        auto result = temp.flashPT(100000.0, 273.15)[0];

        REQUIRE(result[PCTemperature] == Approx(273.15).epsilon(0.001));
        REQUIRE(result[PCPressure] == Approx(100000.0).epsilon(0.001));
        REQUIRE(result[PCCompressibility] == Approx(0.978939).epsilon(0.001));
        REQUIRE(result[PCFugacityCoefficient] == Approx(0.97931788809).epsilon(0.001));
        REQUIRE(result[PCMolarFlow] == Approx(1.0).epsilon(0.001));
        REQUIRE(result[PCMolarVolume] == Approx(0.022233).epsilon(0.001));
        REQUIRE(result[PCEnthalpy] == Approx(-1899.261857).epsilon(0.001));
        REQUIRE(result[PCEntropy] == Approx(-6.496219).epsilon(0.001));
        REQUIRE(result[PCInternalEnergy] == Approx(-4122.525469).epsilon(0.001));
        REQUIRE(result[PCGibbsEnergy] == Approx(-124.819694).epsilon(0.001));
        REQUIRE(result[PCHelmholzEnergy] == Approx(-2348.083305).epsilon(0.001));
    }
}