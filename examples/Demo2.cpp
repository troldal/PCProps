#include <iomanip>
#include <iostream>

#include <library/EquationOfState/EOSPengRobinson.hpp>
#include <library/HeatCapacity/AlyLee.hpp>
#include <library/HeatCapacity/PPDSLiquid.hpp>
#include <library/LiquidVolume/Rackett.hpp>
#include <library/PCComponent.hpp>
#include <library/VaporPressure/AntoineExtended.hpp>

using PCProps::EquationOfState::EOSPengRobinson;
using PCProps::HeatCapacity::AlyLee;
using PCProps::HeatCapacity::PPDSLiquid;
using PCProps::LiquidVolume::Rackett;
using PCProps::VaporPressure::AntoineExtended;

using PCProps::PCComponent;
using PCProps::PCComponentData;
using PCProps::PCEquationOfState;

int main()
{
    using std::get;

    PCComponentData data;

    data.name  = "PROPANE";
    data.name  = "C3H8";
    data.casrn = "74-98-6";

    data.molecularWeight         = 44.096;
    data.boilingTemperature      = 231.050008111111;
    data.freezingTemperature     = 85.15;
    data.criticalTemperature     = 369.83;
    data.criticalPressure        = 4.248E6;
    data.criticalVolume          = 0.0002;
    data.criticalDensity         = 220.48;
    data.criticalCompressibility = 0.27629827986994;
    data.acentricFactor          = 0.1523;

    data.equationOfState          = EOSPengRobinson {};
    data.idealGasCpCorrelation    = AlyLee(AlyLee::CreateFromDIPPR { 0.5192E5, 1.9245E5, 1.6265E3, 1.168E5, 723.6 });
    data.liquidCpCorrelation      = PPDSLiquid(PPDSLiquid::CreateFromDIPPR { 62.983, 113630, 633.21, -873.46, 369.83 });
    data.vaporPressureCorrelation = AntoineExtended(AntoineExtended::CreateFromDIPPR { 59.078, -3492.6, -6.0669, 1.0919E-05, 2 });
    // surfaceTensionFunction
    // heatOfVaporizationFunction
    // vaporThermalConductivityFunction
    // liquidThermalConductivityFunction
    // vaporViscosityFunction
    // liquidViscosityFunction
    data.saturatedLiquidVolumeCorrelation = Rackett(Rackett::CreateFromDIPPR { 1.3757, 0.27453, 369.83, 0.29359 });
    // compressedLiquidVolumeFunction

    auto propane = PCComponent(data);

    std::cout << propane.flashPT(1E5, 298.15) << std::endl;

    return 0;
}