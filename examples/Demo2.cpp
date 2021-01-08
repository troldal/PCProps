#include <iomanip>
#include <iostream>

#include <library/EquationOfState/PengRobinson.hpp>
#include <library/HeatCapacity/AlyLee.hpp>
#include <library/HeatCapacity/PPDSLiquid.hpp>
#include <library/LiquidVolume/Rackett.hpp>
#include <library/PCComponent.hpp>
#include <library/VaporPressure/AntoineExtended.hpp>

using PCProps::EquationOfState::PengRobinson;
using PCProps::HeatCapacity::AlyLee;
using PCProps::HeatCapacity::PPDSLiquid;
using PCProps::LiquidVolume::Rackett;
using PCProps::VaporPressure::AntoineExtended;

using PCProps::PCComponent;
using PCProps::PCComponentData;
using PCProps::PCEquationOfState;
using PCProps::PCPhase;
using PCProps::Utilities::Pressure;
using PCProps::Utilities::Temperature;
using PCProps::Utilities::VaporFraction;

using namespace PCProps;

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

    data.equationOfState                      = PengRobinson {};
    data.idealGasCpCorrelation                = AlyLee(AlyLee::CreateFromDIPPR { 0.5192E5, 1.9245E5, 1.6265E3, 1.168E5, 723.6 });
    data.liquidCpCorrelation                  = PPDSLiquid(PPDSLiquid::CreateFromDIPPR { 62.983, 113630, 633.21, -873.46, 369.83 });
    data.vaporPressureCorrelation             = AntoineExtended(AntoineExtended::CreateFromDIPPR { 59.078, -3492.6, -6.0669, 1.0919E-05, 2 });
    data.surfaceTensionCorrelation            = {};
    data.heatOfVaporizationCorrelation        = {};
    data.vaporThermalConductivityCorrelation  = {};
    data.liquidThermalConductivityCorrelation = {};
    data.vaporViscosityCorrelation            = {};
    data.liquidViscosityCorrelation           = {};
    data.saturatedLiquidVolumeCorrelation     = Rackett(Rackett::CreateFromDIPPR { 1.3757, 0.27453, 369.83, 0.29359 });
    data.compressedLiquidVolumeCorrelation    = {};

    auto propane = PCComponent(data);

    std::cout << propane.saturationPressure(369.83) << std::endl;

    for (const auto& phase : propane.flash(Pressure(4.248E6), VaporFraction(0.5))) std::cout << phase << std::endl;

    for (const auto& phase : propane.flash(Temperature(369.8), VaporFraction(0.5))) std::cout << phase << std::endl;
    //
    //    auto psat = AntoineExtended(AntoineExtended::CreateFromDIPPR { 59.078, -3492.6, -6.0669, 1.0919E-05, 2 });
    //    std::cout << psat(500) << std::endl;
    //
    //    for (const auto& phase : propane.flash(Pressure(2.44656E7), Temperature(500)))
    //        std::cout << phase << std::endl;

    return 0;
}