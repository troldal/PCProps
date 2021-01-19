#include <iomanip>
#include <iostream>

#include <Component/PureComponent.hpp>
#include <Fluid.hpp>
#include <HeatCapacity/AlyLee.hpp>
#include <HeatCapacity/PPDSLiquid.hpp>
#include <PengRobinson/PengRobinson.hpp>
#include <SaturatedLiquidVolume/Rackett.hpp>
#include <VaporPressure/AmbroseWalton.hpp>
#include <VaporPressure/AntoineExtended.hpp>
#include <Viscosity/DIPPR102.hpp>
#include <Viscosity/KirchhoffExtended.hpp>
#include <Viscosity/Lucas.hpp>
#include <common/PropertyData.hpp>

using PCProps::VaporPressure::AmbroseWalton;

using PCProps::EquationOfState::PengRobinson;
using PCProps::HeatCapacity::AlyLee;
using PCProps::HeatCapacity::PPDSLiquid;
using PCProps::LiquidVolume::Rackett;
using PCProps::VaporPressure::AntoineExtended;

using PCProps::PCComponentData;
using PCProps::PCComponentData;
using PCProps::PCComponentData;
using PCProps::PCPhase;
using PCProps::Viscosity::Lucas;
using PCProps::Viscosity::DIPPR102;
using PCProps::Viscosity::KirchhoffExtended;

using PCProps::Enthalpy;
using PCProps::Entropy;
using PCProps::Pressure;
using PCProps::Temperature;

using namespace PCProps;

int main()
{
    using std::get;

    PCComponentData data;

    data.name  = "PROPANE";
    data.name  = "C3H8";
    data.casrn = "74-98-6";

    data.molarWeight             = 44.096;
    data.boilingTemperature      = 231.050008111111;
    data.freezingTemperature     = 85.15;
    data.criticalTemperature     = 369.83;
    data.criticalPressure        = 4.248E6;
    data.criticalVolume          = 0.0002;
    data.criticalDensity         = 220.48;
    data.criticalCompressibility = 0.27629827986994;
    data.acentricFactor          = 0.1523;
    data.dipoleMoment            = 0.083;

    auto psat = AmbroseWalton(369.83, 4.248E6, 0.1523);

    data.idealGasCpCorrelation                = AlyLee(AlyLee::CreateFromDIPPR { 0.5192E5, 1.9245E5, 1.6265E3, 1.168E5, 723.6 });
    data.liquidCpCorrelation                  = PPDSLiquid(PPDSLiquid::CreateFromDIPPR { 62.983, 113630, 633.21, -873.46, 369.83 });
    data.vaporPressureCorrelation             = AntoineExtended(AntoineExtended::CreateFromDIPPR { 59.078, -3492.6, -6.0669, 1.0919E-05, 2 });
    data.surfaceTensionCorrelation            = {};
    data.heatOfVaporizationCorrelation        = {};
    data.saturatedVaporThermalConductivityCorrelation = {};
    data.saturatedLiquidThermalConductivityCorrelation = {};
    data.saturatedVaporViscosityCorrelation   = Lucas(369.83, 4.248E6, 0.2763, 44.096, 0.083);    //DIPPR102(4.9054E-08, 0.90125);
    data.saturatedLiquidViscosityCorrelation  = KirchhoffExtended(-17.156,646.25,1.1101,-7.3439E-11, 4);
    data.saturatedLiquidVolumeCorrelation     = Rackett(Rackett::CreateFromDIPPR { 1.3757, 0.27453, 369.83, 0.29359 });

    auto propane = PureComponent(data);
    auto fluid = Fluid(propane, PengRobinson{});

    std::cout << "Propane at 25 C and 2 bar: " << std::endl;
    auto a = fluid.flash(Pressure(2E5), Temperature(298.15));
    for (const auto& phase : a) std::cout << phase << std::endl;
    std::cout << "==================================================" << std::endl;

    std::cout << "Compression to 10 bar: " << std::endl;
    auto b = fluid.flash(Pressure(10E5), Entropy(PCPhase(a[0]).entropy()));
    for (const auto& phase : b) std::cout << phase << std::endl;
    std::cout << "==================================================" << std::endl;

    std::cout << "Cooling to 25 C: " << std::endl;
    auto c = fluid.flash(Pressure(10E5), Temperature(298.15));
    std::cout << "dT: " << c[0][PCTemperature] - b[0][PCTemperature] << std::endl;
    std::cout << "dH: " << c[0][PCEnthalpy] - b[0][PCEnthalpy] << std::endl;
    for (const auto& phase : c) std::cout << phase << std::endl;
    std::cout << "==================================================" << std::endl;

    std::cout << "Throttling to 2 bar: " << std::endl;
    auto d = fluid.flash(Pressure(2E5), Enthalpy(c[0][PCEnthalpy]));
    for (const auto& phase : d) std::cout << phase << std::endl;
    std::cout << "==================================================" << std::endl;

    return 0;
}