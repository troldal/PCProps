#include <iomanip>
#include <iostream>

#include <EOSLib.hpp>
#include <Fluid.hpp>
#include <PropertyLib.hpp>

#include <common/PropertyData.hpp>

using PCProps::VaporPressure::AmbroseWalton;

using PCProps::EquationOfState::PengRobinson;
using PCProps::HeatCapacity::AlyLee;
using PCProps::HeatCapacity::PPDSLiquid;
using PCProps::HeatCapacity::Polynomial;
using PCProps::LiquidVolume::Rackett;
using PCProps::LiquidVolume::YenWoods;
using PCProps::VaporPressure::AntoineExtended;

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

//using PCProps::UnitOps::Pipe;

using namespace PCProps;

int main()
{
    using std::get;

    PCComponentData data;

    data.name  = "WATER";
    data.name  = "H2=";
    data.casrn = "7732-18-5";

    data.molarWeight             = 18.015;
    data.boilingTemperature      = 373.15;
    data.freezingTemperature     = 273.15;
    data.criticalTemperature     = 647.10;
    data.criticalPressure        = 22064000.0;
    data.criticalVolume          = 0.0000559472;
    data.criticalDensity         = 322.00003;
    data.criticalCompressibility = 0.229435018515262;
    data.acentricFactor          = 0.3449;
    data.dipoleMoment            = 1.8546;

    data.idealGasCpCorrelation                = AlyLee(AlyLee::CreateFromDIPPR { 0.33363E5, 0.2679E5, 2.6105E3, 0.08896E5, 1169.0 });
    data.liquidCpCorrelation                  = Polynomial(Polynomial::CreateFromDIPPR {276370, -2090.1, 8.125, -0.014116, 9.3701E-06});
    data.vaporPressureCorrelation             = AntoineExtended(AntoineExtended::CreateFromDIPPR { 73.649, -7258.2, -7.3037, 4.1653E-06, 2 });
    data.surfaceTensionCorrelation            = {};
    data.heatOfVaporizationCorrelation        = {};
    data.satVaporThermalConductivityCorrelation        = {};
    data.satLiquidThermalConductivityCorrelation       = {};
    data.satVaporViscosityCorrelation                  = DIPPR102(1.7096E-08, 1.1146);
    data.satLiquidViscosityCorrelation                 = KirchhoffExtended(-52.843,3703.6,5.866,-5.879E-29, 10);
    data.satLiquidVolumeCorrelation                    = YenWoods(YenWoods::CreateFromYenWoodsEstimation{647.10, 0.0000559472, 0.229435018515262});

    auto propane = PureComponent(data);
    auto fluid = Fluid(propane, PengRobinson{});

    std::cout << "Water at 25 C and 1 atm: " << std::endl;
    auto a = fluid.flash(Pressure(101325), Temperature(298.15));
    for (const auto& phase : a) std::cout << phase << std::endl;
    std::cout << "==================================================" << std::endl;

//    auto pipe = Pipe(9000, 152.0/1000, 7, 0);
//    std::cout << pipe.computeOutletPressure(a, 2620) << std::endl;


    return 0;
}