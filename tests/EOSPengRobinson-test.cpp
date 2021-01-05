//
// Created by Kenneth Balslev on 03/01/2021.
//

#include <catch.hpp>
#include <library/EquationOfState/EOSPengRobinson.hpp>
#include <library/EquationOfState/EOSUtilities.hpp>
#include <library/HeatCapacity/AlyLee.hpp>
#include <library/VaporPressure/AmbroseWalton.hpp>

using PCProps::EquationOfState::EOSPengRobinson;
using PCProps::EquationOfState::PhaseDataElement;
using PCProps::HeatCapacity::AlyLee;
using PCProps::VaporPressure::AmbroseWalton;

using namespace PCProps::EquationOfState;

TEST_CASE("EOSPengRobinson Test")
{
    auto tc    = 369.83;
    auto pc    = 4.248E6;
    auto omega = 0.1523;
    auto mw    = 44.096;

    auto            PSat = AmbroseWalton(tc, pc, omega);
    auto            igCp = AlyLee(0.5192E5, 1.9245E5, 1.6265E3, 1.168E5, 723.6);
    EOSPengRobinson propane(tc, pc, omega, mw, PSat, igCp);

    SECTION("PT Flash of propane @ 273.15 K and 1 bar")
    {
        auto result = propane.flashPT(100000.0, 273.15);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;
        REQUIRE(get<Temperature>(result[0]) == Approx(273.15).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.978939).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(97931.788809).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.022233).epsilon(0.001));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-1899.261857).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-6.496219).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-4122.525469).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-124.819694).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-2348.083305).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 223.15 K and 1 bar - Compressed Liquid")
    {
        auto result = propane.flashPT(100000.0, 223.15);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;
        REQUIRE(get<Temperature>(result[0]) == Approx(223.15).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.003766).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(69350.105382).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.000070).epsilon(0.01));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-24084.251701).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-101.724616).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-24091.238463).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-1384.403549).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-1391.390312).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 323.15 K and 1 bar - Superheated Vapor")
    {
        auto result = propane.flashPT(100000.0, 323.15);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;
        REQUIRE(get<Temperature>(result[0]) == Approx(323.15).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.986760).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(98690.029547).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.026512).epsilon(0.001));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(1803.116466).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(5.926803).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-848.128978).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-112.130043).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-2763.375487).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 469.8 K and 1 bar - T > Tc")
    {
        auto result = propane.flashPT(100000.0, 469.8);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;
        REQUIRE(get<Temperature>(result[0]) == Approx(469.8).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.995756).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(99576.428359).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.038896).epsilon(0.001));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(15577.940459).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(40.740519).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(11688.385381).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-3561.955560).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-7451.510638).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 323.15 K and 50 bar - P > Pc")
    {
        auto result = propane.flashPT(5000000.0, 323.15);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;
        REQUIRE(get<Temperature>(result[0]) == Approx(323.15).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(5000000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.171348).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(1507355.182464).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.000092).epsilon(0.001));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-13089.953760).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-62.826659).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-13550.335573).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(7212.481187).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(6752.099374).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 469.8 K and 50 bar - Supercritical")
    {
        auto result = propane.flashPT(5000000.0, 469.8);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;
        REQUIRE(get<Temperature>(result[0]) == Approx(469.8).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(5000000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.794804).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(4053007.334728).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.000621).epsilon(0.001));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(12596.727773).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(3.579012).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(9492.115417).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(10915.307766).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(7810.695410).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 1 K and 1 Pa - Low P & T")
    {
        auto result = propane.flashPT(1, 1);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;
        REQUIRE(get<Temperature>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.000007).epsilon(0.1));
        REQUIRE(get<Fugacity>(result[0]) == Approx(0.0).epsilon(0.01));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.000056).epsilon(0.01));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-45132.566226).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-952.525226).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-45132.566282).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-44180.041000).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-44180.041056).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ Tc and Pc - @ Critical Point")
    {
        auto result = propane.flashPT(pc, tc);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;
        REQUIRE(get<Temperature>(result[0]) == Approx(369.83).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(4248000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.305823).epsilon(0.1));
        REQUIRE(get<Fugacity>(result[0]) == Approx(2729914.681402).epsilon(0.01));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.000221).epsilon(0.01));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-2372.295636).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-32.215245).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-3312.681919).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(9541.868409).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(8601.482126).epsilon(0.001));
    }

    SECTION("Px Flash of propane @ 1 bar and x = 0.0 - Saturated Liquid")
    {
        auto result = propane.flashPx(100000.0, 0.0);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;

        REQUIRE(get<Temperature>(result[0]) == Approx(230.657526).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.003697).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(96789.922051).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.000071).epsilon(0.01));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-23379.777150).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-98.619779).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-23386.868147).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-632.382962).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-639.473959).epsilon(0.001));
    }

    SECTION("Px Flash of propane @ 1 bar and x = 1.0 - Saturated Vapor")
    {
        auto result = propane.flashPx(100000.0, 1.0);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;

        REQUIRE(get<Temperature>(result[0]) == Approx(230.657526).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.966927).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(96789.922052).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.018544).epsilon(0.001));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-4667.418728).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-17.493623).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-6521.785497).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-632.382962).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-2486.749730).epsilon(0.001));
    }

    SECTION("Px Flash of propane @ 1 bar and x = 0.5 - Two-Phase")
    {
        auto result = propane.flashPx(100000.0, 0.5);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;

        REQUIRE(get<Temperature>(result[0]) == Approx(230.657526).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.966927).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(96789.922052).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(0.5).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.018544).epsilon(0.001));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-4667.418728).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-17.493623).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-6521.785497).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-632.382962).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-2486.749730).epsilon(0.001));

        REQUIRE(get<Temperature>(result[1]) == Approx(230.657526).epsilon(0.001));
        REQUIRE(get<Pressure>(result[1]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[1]) == Approx(0.003697).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[1]) == Approx(96789.922051).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[1]) == Approx(0.5).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[1]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[1]) == Approx(0.000071).epsilon(0.01));
        REQUIRE(get<Enthalpy>(result[1]) == Approx(-23379.777150).epsilon(0.001));
        REQUIRE(get<Entropy>(result[1]) == Approx(-98.619779).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[1]) == Approx(-23386.868147).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[1]) == Approx(-632.382962).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[1]) == Approx(-639.473959).epsilon(0.001));
    }

    SECTION("Tx Flash of propane @ 223.15 K and x = 0.0 - Saturated Liquid")
    {
        auto result = propane.flashTx(223.15, 0.0);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;

        REQUIRE(get<Temperature>(result[0]) == Approx(223.15).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(71038.119363).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.002675).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(69274.509900).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.000070).epsilon(0.01));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-24085.412742).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-101.720751).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-24090.376323).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-1386.427112).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-1391.390693).epsilon(0.001));
    }

    SECTION("Tx Flash of propane @ 223.15 K and x = 1.0 - Saturated Vapor")
    {
        auto result = propane.flashTx(223.15, 1.0);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;

        REQUIRE(get<Temperature>(result[0]) == Approx(223.15).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(71038.119363).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.974596).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(69274.509900).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.025454).epsilon(0.001));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-5075.219460).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-16.530551).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-6883.457001).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-1386.427112).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-3194.664652).epsilon(0.001));
    }

    SECTION("Tx Flash of propane @ 223.15 K and x = 0.5 - Two-Phase")
    {
        auto result = propane.flashTx(223.15, 0.5);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;

        REQUIRE(get<Temperature>(result[0]) == Approx(223.15).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(71038.119363).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.974596).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(69274.509900).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(0.5).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.025454).epsilon(0.001));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-5075.219460).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-16.530551).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-6883.457001).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-1386.427112).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-3194.664652).epsilon(0.001));

        REQUIRE(get<Temperature>(result[1]) == Approx(223.15).epsilon(0.001));
        REQUIRE(get<Pressure>(result[1]) == Approx(71038.119363).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[1]) == Approx(0.002675).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[1]) == Approx(69274.509900).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[1]) == Approx(0.5).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[1]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[1]) == Approx(0.000070).epsilon(0.01));
        REQUIRE(get<Enthalpy>(result[1]) == Approx(-24085.412742).epsilon(0.001));
        REQUIRE(get<Entropy>(result[1]) == Approx(-101.720751).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[1]) == Approx(-24090.376323).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[1]) == Approx(-1386.427112).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[1]) == Approx(-1391.390693).epsilon(0.001));
    }

    SECTION("PH Flash of propane @ 1 bar and H = -14000 - Two-Phase")
    {
        auto result = propane.flashPH(100000, -14000);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;

        REQUIRE(get<Temperature>(result[0]) == Approx(230.657526).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.966927).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(96789.922052).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(0.501261).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.018544).epsilon(0.001));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-4667.418728).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-17.493623).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-6521.785497).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-632.382962).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-2486.749730).epsilon(0.001));

        REQUIRE(get<Temperature>(result[1]) == Approx(230.657526).epsilon(0.001));
        REQUIRE(get<Pressure>(result[1]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[1]) == Approx(0.003697).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[1]) == Approx(96789.922051).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[1]) == Approx(0.498739).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[1]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[1]) == Approx(0.000071).epsilon(0.01));
        REQUIRE(get<Enthalpy>(result[1]) == Approx(-23379.777150).epsilon(0.001));
        REQUIRE(get<Entropy>(result[1]) == Approx(-98.619779).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[1]) == Approx(-23386.868147).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[1]) == Approx(-632.382962).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[1]) == Approx(-639.473959).epsilon(0.001));
    }

    SECTION("PH Flash of propane @ 1 bar and H = 2000 - Superheated Vapor")
    {
        auto result = propane.flashPH(100000, 2000);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;

        REQUIRE(get<Temperature>(result[0]) == Approx(325.627654).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.987043).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(98717.657406).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.026723).epsilon(0.001));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(2000).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(6.533741).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-672.338547).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-127.566650).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-2799.905197).epsilon(0.001));
    }

    SECTION("PH Flash of propane @ 1 bar and H = 16000 - T > Tc")
    {
        auto result = propane.flashPH(100000, 16000);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;

        REQUIRE(get<Temperature>(result[0]) == Approx(473.704635).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.995872).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(99587.960365).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.039223).epsilon(0.001));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(16000).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(41.635185).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(12077.660810).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-3722.780350).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-7645.119541).epsilon(0.001));
    }

    SECTION("PH Flash of propane @ 1 bar and H = -24000 - Compressed Liquid")
    {
        auto result = propane.flashPH(100000, -24000);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;

        REQUIRE(get<Temperature>(result[0]) == Approx(224.055958).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.003757).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(72299.863177).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.000070).epsilon(0.01));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-24000).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-101.347824).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-24006.998906).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-1292.416115).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-1299.415022).epsilon(0.001));
    }

    SECTION("PH Flash of propane @ 50 bar and H = -20000 - P > Pc")
    {
        auto result = propane.flashPH(5000000, -20000);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;

        REQUIRE(get<Temperature>(result[0]) == Approx(263.308763).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(5000000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.171527).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(374944.969272).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.000075).epsilon(0.01));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-20000).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-86.373130).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-20375.518173).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(2742.802087).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(2367.283914).epsilon(0.001));
    }

    SECTION("PH Flash of propane @ 50 bar and H = 12000 - Supercritical")
    {
        auto result = propane.flashPH(5000000, 12000);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;

        REQUIRE(get<Temperature>(result[0]) == Approx(464.976124).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(5000000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.785893).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(4020028.107737).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.000608).epsilon(0.01));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(12000).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(2.302270).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(8961.717250).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(10929.499302).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(7891.216551).epsilon(0.001));
    }

    SECTION("PS Flash of propane @ 1 bar and S = -60 - Two-Phase")
    {
        auto result = propane.flashPS(100000, -60);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;

        REQUIRE(get<Temperature>(result[0]) == Approx(230.657526).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.966927).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(96789.922052).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(0.476046).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.018544).epsilon(0.001));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-4667.418728).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-17.493623).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-6521.785497).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-632.382962).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-2486.749730).epsilon(0.001));

        REQUIRE(get<Temperature>(result[1]) == Approx(230.657526).epsilon(0.001));
        REQUIRE(get<Pressure>(result[1]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[1]) == Approx(0.003697).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[1]) == Approx(96789.922051).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[1]) == Approx(0.523954).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[1]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[1]) == Approx(0.000071).epsilon(0.01));
        REQUIRE(get<Enthalpy>(result[1]) == Approx(-23379.777150).epsilon(0.001));
        REQUIRE(get<Entropy>(result[1]) == Approx(-98.619779).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[1]) == Approx(-23386.868147).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[1]) == Approx(-632.382962).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[1]) == Approx(-639.473959).epsilon(0.001));
    }

    SECTION("PS Flash of propane @ 1 bar and S = 7 - Superheated Vapor")
    {
        auto result = propane.flashPS(100000, 7);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;

        REQUIRE(get<Temperature>(result[0]) == Approx(327.532918).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.987255).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(98738.413477).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.026886).epsilon(0.001));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(2152.271027).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(7.0).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-536.281691).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-140.459399).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-2829.012117).epsilon(0.001));
    }

    SECTION("PS Flash of propane @ 1 bar and S = 40 - T > Tc")
    {
        auto result = propane.flashPS(100000, 40);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;

        REQUIRE(get<Temperature>(result[0]) == Approx(466.573938).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.995657).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(99566.632742).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.038625).epsilon(0.001));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(15231.239828).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(40.0).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(11368.776266).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-3431.717744).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-7294.181306).epsilon(0.001));
    }

    SECTION("PS Flash of propane @ 1 bar and S = -100 - Compressed Liquid")
    {
        auto result = propane.flashPS(100000, -100);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;

        REQUIRE(get<Temperature>(result[0]) == Approx(227.308593).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.003727).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(83688.539684).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.000070).epsilon(0.01));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-23695.822221).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-100.0).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-23702.865692).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-964.962920).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-972.006390).epsilon(0.001));
    }

    SECTION("PS Flash of propane @ 50 bar and S = -85 - P > Pc")
    {
        auto result = propane.flashPS(5000000, -85);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;

        REQUIRE(get<Temperature>(result[0]) == Approx(266.835388).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(5000000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.170734).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(416727.575115).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.000076).epsilon(0.01));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-19636.022175).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-85.0).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-20014.811861).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(3044.985805).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(2666.196118).epsilon(0.001));
    }

    SECTION("PS Flash of propane @ 50 bar and S = 3 - Supercritical")
    {
        auto result = propane.flashPS(5000000, 3.0);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;

        REQUIRE(get<Temperature>(result[0]) == Approx(467.605973).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(5000000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.790814).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(4038164.623638).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.000615).epsilon(0.01));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(12325.343803).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(3.0).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(9250.744824).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(10922.525896).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(7847.926917).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 273.15 K and 1 bar - Copy Constructor")
    {
        auto temp   = propane;
        auto result = temp.flashPT(100000.0, 273.15);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;
        REQUIRE(get<Temperature>(result[0]) == Approx(273.15).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.978939).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(97931.788809).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.022233).epsilon(0.001));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-1899.261857).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-6.496219).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-4122.525469).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-124.819694).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-2348.083305).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 273.15 K and 1 bar - Move Constructor")
    {
        auto temp   = std::move(propane);
        auto result = temp.flashPT(100000.0, 273.15);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;
        REQUIRE(get<Temperature>(result[0]) == Approx(273.15).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.978939).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(97931.788809).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.022233).epsilon(0.001));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-1899.261857).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-6.496219).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-4122.525469).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-124.819694).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-2348.083305).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 273.15 K and 1 bar - Copy Assignment")
    {
        EOSPengRobinson temp {};
        temp        = propane;
        auto result = temp.flashPT(100000.0, 273.15);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;
        REQUIRE(get<Temperature>(result[0]) == Approx(273.15).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.978939).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(97931.788809).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.022233).epsilon(0.001));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-1899.261857).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-6.496219).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-4122.525469).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-124.819694).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-2348.083305).epsilon(0.001));
    }

    SECTION("PT Flash of propane @ 273.15 K and 1 bar - Move Assignment")
    {
        EOSPengRobinson temp {};
        temp        = std::move(propane);
        auto result = temp.flashPT(100000.0, 273.15);

        using PCProps::EquationOfState::PhaseDataElement;
        using std::get;
        REQUIRE(get<Temperature>(result[0]) == Approx(273.15).epsilon(0.001));
        REQUIRE(get<Pressure>(result[0]) == Approx(100000.0).epsilon(0.001));
        REQUIRE(get<Compressibility>(result[0]) == Approx(0.978939).epsilon(0.001));
        REQUIRE(get<Fugacity>(result[0]) == Approx(97931.788809).epsilon(0.001));
        REQUIRE(get<MolarFraction>(result[0]) == Approx(1.0).epsilon(0.001));
        REQUIRE(get<MolecularWeight>(result[0]) == Approx(44.096000).epsilon(0.001));
        REQUIRE(get<Volume>(result[0]) == Approx(0.022233).epsilon(0.001));
        REQUIRE(get<Enthalpy>(result[0]) == Approx(-1899.261857).epsilon(0.001));
        REQUIRE(get<Entropy>(result[0]) == Approx(-6.496219).epsilon(0.001));
        REQUIRE(get<InternalEnergy>(result[0]) == Approx(-4122.525469).epsilon(0.001));
        REQUIRE(get<GibbsEnergy>(result[0]) == Approx(-124.819694).epsilon(0.001));
        REQUIRE(get<HelmholzEnergy>(result[0]) == Approx(-2348.083305).epsilon(0.001));
    }
}