#include <iomanip>
#include <iostream>

#include <EOSLib.hpp>
#include <Fluid.hpp>
#include <PropertyLib.hpp>
#include <PureComponentFactory.hpp>
#include <DataSource.hpp>

#include <FluidProperties.hpp>

using PCProps::VaporPressure::AmbroseWalton;

using PCProps::EquationOfState::PengRobinson;
using PCProps::HeatCapacity::AlyLee;
using PCProps::HeatCapacity::PPDSLiquid;
using PCProps::LiquidVolume::Rackett;
using PCProps::VaporPressure::AntoineExtended;

using PCProps::Viscosity::Lucas;
using PCProps::Viscosity::DIPPR102;
using PCProps::Viscosity::KirchhoffExtended;

using namespace PCProps;

int main()
{
    auto ds = DataSource("Mini PCD.xlsx");
    auto pcf = PureComponentFactory(ds.load());
    auto pc = pcf.makeComponent("74-98-6");

    pc.addDataItem("IdealGasCp", AlyLee(AlyLee::CreateFromDIPPR { 0.5192E5, 1.9245E5, 1.6265E3, 1.168E5, 723.6 }));
    pc.addDataItem("LiquidCp", PPDSLiquid(PPDSLiquid::CreateFromDIPPR { 62.983, 113630, 633.21, -873.46, 369.83 }));
    pc.addDataItem("VaporPressure", AntoineExtended(AntoineExtended::CreateFromDIPPR { 59.078, -3492.6, -6.0669, 1.0919E-05, 2 }));
    pc.addDataItem("SaturatedVaporViscosity", Viscosity::Lucas(369.83, 4.248E6, 0.2763, 44.096, 0.083));
    pc.addDataItem("SaturatedLiquidViscosity", KirchhoffExtended(-17.156, 646.25, 1.1101, -7.3439E-11, 4));
    pc.addDataItem("SaturatedLiquidVolume", Rackett(Rackett::CreateFromDIPPR { 1.3757, 0.27453, 369.83, 0.29359 }));

    auto fluid = Fluid(pc, PengRobinson{});

    std::cout << "Propane at 25 C and 2 bar: " << std::endl;
    auto a = FluidProperties(fluid.flashPT((2E5), (298.15)));
    //std::cout << a << std::endl;
    a.print(std::cout);
    std::cout << "==============================================================================" << std::endl;

    std::cout << "Compression to 10 bar: " << std::endl;
    auto b = FluidProperties(fluid.flashPS((10E5), (a[0].Entropy)));
    //std::cout << b << std::endl;
    b.print(std::cout);
    std::cout << "==============================================================================" << std::endl;

    std::cout << "Cooling to 25 C: " << std::endl;
    auto c = FluidProperties(fluid.flashPT((10E5), (298.15)));
    //std::cout << c << std::endl;
    c.print(std::cout);
    std::cout << "==============================================================================" << std::endl;

    std::cout << "Throttling to 2 bar: " << std::endl;
    auto d = FluidProperties(fluid.flashPH((2E5), (c[0].Enthalpy)));
    //std::cout << d << std::endl;
    d.print(std::cout);
    std::cout << "==============================================================================" << std::endl;

    return 0;
}