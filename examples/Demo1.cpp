#include <iomanip>
#include <iostream>

#include <EOSLib.hpp>
#include <Fluid.hpp>
#include <PropertyLib.hpp>

#include <common/PhaseProperties.hpp>
#include <json/json.hpp>

using PCProps::VaporPressure::AmbroseWalton;

using PCProps::EquationOfState::PengRobinson;
using PCProps::HeatCapacity::AlyLee;
using PCProps::HeatCapacity::PPDSLiquid;
using PCProps::LiquidVolume::Rackett;
using PCProps::VaporPressure::AntoineExtended;

using PCProps::PCComponentData;
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

    using PCProps::LiquidVolume::Thomson;
    data.compressedLiquidVolume = Thomson(data.criticalTemperature.value(),
                                          data.criticalPressure.value(),
                                          data.acentricFactor.value());

    using namespace PCProps::CompressedLiquidViscosity;
    data.compressedLiquidViscosity = CompressedLiquidViscosity::Lucas(data.criticalTemperature.value(),
                                                                      data.criticalPressure.value(),
                                                                      data.acentricFactor.value());

    using namespace PCProps::CompressedVaporViscosity;
    data.compressedVaporViscosity = CompressedVaporViscosity::Lucas(data.criticalTemperature.value(),
                                                                    data.criticalPressure.value(),
                                                                    data.criticalCompressibility.value(),
                                                                    data.molarWeight.value(),
                                                                    data.dipoleMoment.value());

    auto fluid = Fluid(PureComponent{data}, PengRobinson{});

    std::cout << "Propane at 25 C and 2 bar: " << std::endl;
    auto a = nlohmann::json::parse(fluid.flashPT((2E5), (298.15)));
    for (const auto& phase : a) std::cout << PhaseProperties {phase} << std::endl;
    std::cout << "==================================================" << std::endl;

    std::cout << "Compression to 10 bar: " << std::endl;
    auto b = nlohmann::json::parse(fluid.flashPS((10E5), (a[0]["Entropy"])));
    for (const auto& phase : b) std::cout << PhaseProperties {phase} << std::endl;
    std::cout << "==================================================" << std::endl;

    std::cout << "Cooling to 25 C: " << std::endl;
    auto c = nlohmann::json::parse(fluid.flashPT((10E5), (298.15)));
    std::cout << "dT: " << c[0]["Temperature"].get<double>() - b[0]["Temperature"].get<double>() << std::endl;
    std::cout << "dH: " << c[0]["Enthalpy"].get<double>() - b[0]["Enthalpy"].get<double>() << std::endl;
    for (const auto& phase : c) std::cout << PhaseProperties {phase} << std::endl;
    std::cout << "==================================================" << std::endl;

    std::cout << "Throttling to 2 bar: " << std::endl;
    auto d = nlohmann::json::parse(fluid.flashPH((2E5), (c[0]["Enthalpy"])));
    for (const auto& phase : d) std::cout << PhaseProperties {phase} << std::endl;
    std::cout << "==================================================" << std::endl;

    return 0;
}