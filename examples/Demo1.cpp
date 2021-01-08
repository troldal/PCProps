#include <iomanip>
#include <iostream>

#include <library/EquationOfState/PengRobinson.hpp>
#include <library/HeatCapacity/AlyLee.hpp>
#include <library/HeatCapacity/PPDSLiquid.hpp>
#include <library/LiquidVolume/Rackett.hpp>
#include <library/PCComponent.hpp>
#include <library/PCPropsData.hpp>
#include <library/VaporPressure/AmbroseWalton.hpp>
#include <library/VaporPressure/AntoineExtended.hpp>

using PCProps::EquationOfState::PengRobinson;
using PCProps::HeatCapacity::AlyLee;
using PCProps::HeatCapacity::PPDSLiquid;
using PCProps::LiquidVolume::Rackett;
using PCProps::VaporPressure::AmbroseWalton;
using PCProps::VaporPressure::AntoineExtended;

using PCProps::PCComponent;
using PCProps::PCComponentData;
using PCProps::PCEquationOfState;
using PCProps::PCPhase;

using PCProps::Utilities::Enthalpy;
using PCProps::Utilities::Entropy;
using PCProps::Utilities::Pressure;
using PCProps::Utilities::Temperature;

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

    std::cout << "Propane at 25 C and 2 bar: " << std::endl;
    auto a = propane.flash(Pressure(2E5), Temperature(298.15));
    for (const auto& phase : a) std::cout << phase << std::endl;
    std::cout << "==================================================" << std::endl;

    std::cout << "Compression to 10 bar: " << std::endl;
    auto b = propane.flash(Pressure(10E5), Entropy(PCPhase(a[0]).entropy()));
    for (const auto& phase : b) std::cout << phase << std::endl;
    std::cout << "==================================================" << std::endl;

    std::cout << "Cooling to 25 C: " << std::endl;
    auto c = propane.flash(Pressure(10E5), Temperature(298.15));
    std::cout << "dT: " << c[0][PCTemperature] - b[0][PCTemperature] << std::endl;
    std::cout << "dH: " << c[0][PCEnthalpy] - b[0][PCEnthalpy] << std::endl;
    for (const auto& phase : c) std::cout << phase << std::endl;
    std::cout << "==================================================" << std::endl;

    std::cout << "Throttling to 2 bar: " << std::endl;
    auto d = propane.flash(Pressure(2E5), Enthalpy(c[0][PCEnthalpy]));
    for (const auto& phase : d) std::cout << phase << std::endl;
    std::cout << "==================================================" << std::endl;

    return 0;
}