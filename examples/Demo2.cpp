#include <iomanip>
#include <iostream>

#include <EOSLib.hpp>
#include <Fluid.hpp>
#include <PropertyLib.hpp>

using PCProps::EquationOfState::PengRobinson;
using PCProps::HeatCapacity::AlyLee;
using PCProps::HeatCapacity::PPDSLiquid;
using PCProps::LiquidVolume::Rackett;
using PCProps::VaporPressure::AntoineExtended;

using PCProps::PCComponentData;
using PCProps::PCPhase;
using PCProps::Viscosity::Lucas;
using PCProps::Viscosity::DIPPR102;
using PCProps::Viscosity::KirchhoffExtended;

using namespace PCProps;

int main()
{
    using std::get;

    PCComponentData data;

    data.name  = "PROPANE";
    data.name  = "C3H8";
    data.casrn = "74-98-6";

    data.molarWeight             = 44.096;
    data.boilingTemperature      = 231.05;
    data.freezingTemperature     = 85.15;
    data.criticalTemperature     = 369.83;
    data.criticalPressure        = 4.248E6;
    data.criticalVolume          = 0.0002;
    data.criticalDensity         = 220.48;
    data.criticalCompressibility = 0.2763;
    data.acentricFactor          = 0.1523;
    data.dipoleMoment            = 0.083;

    data.idealGasCpCorrelation                   = AlyLee(AlyLee::CreateFromDIPPR { 0.5192E5, 1.9245E5, 1.6265E3, 1.168E5, 723.6 });
    data.liquidCpCorrelation                     = PPDSLiquid(PPDSLiquid::CreateFromDIPPR { 62.983, 113630, 633.21, -873.46, 369.83 });
    data.vaporPressureCorrelation                = AntoineExtended(AntoineExtended::CreateFromDIPPR { 59.078, -3492.6, -6.0669, 1.0919E-05, 2 });
    data.surfaceTensionCorrelation               = {};
    data.heatOfVaporizationCorrelation           = {};
    data.satVaporThermalConductivityCorrelation  = {};
    data.satLiquidThermalConductivityCorrelation = {};
    data.satVaporViscosityCorrelation            = Lucas(369.83, 4.248E6, 0.2763, 44.096, 0.083);    // DIPPR102(4.9054E-08, 0.90125);
    data.satLiquidViscosityCorrelation           = KirchhoffExtended(-17.156, 646.25, 1.1101, -7.3439E-11, 4);
    data.satLiquidVolumeCorrelation              = Rackett(Rackett::CreateFromDIPPR { 1.3757, 0.27453, 369.83, 0.29359 });

    auto propane = PureComponent(data);
    auto fluid = Fluid(propane, PengRobinson{});

//    auto phases = fluid.flashPT(101325, 200.0);
//    for (const auto& phase : phases) std::cout << phase << std::endl;
//
//    double volume = 0.0;
//    for (const auto& phase : phases) {
//        volume += phase[PCMolarVolume] * phase[PCMolarFlow];
//    }

//    for (const auto& phase : fluid.flashPT(100000.0, 330.0)) std::cout << phase << std::endl;
    for (const auto& phase : fluid.flashTV(468.0, 0.000615)) std::cout << phase << std::endl;

    return 0;
}