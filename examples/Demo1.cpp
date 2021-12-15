#include <iomanip>
#include <iostream>

#include <EOSLib.hpp>
#include <Fluid.hpp>
#include <PureComponentFactory.hpp>
#include <DataSource.hpp>
#include <FluidProperties.hpp>

using PCProps::EquationOfState::PengRobinson;
using namespace PCProps;

int main()
{
    auto ds = DataSource("Mini PCD.xlsx");
    auto pcf = PureComponentFactory(ds.load());
    auto pc = pcf.makeComponent("132259-10-0");
    auto fluid = Fluid(pc, PengRobinson{});

//    FluidProperties(fluid.flash("Tx", 447.3, 0.5)).print(std::cout);
//    FluidProperties(fluid.flash("Px", 4247999.0, 0.5)).print(std::cout);

    auto nbp = 59.0;
    auto tc = 132.44;
    auto diff = (tc - nbp) / 99.0;
    for (int i = 0; i < 100; ++i) {
        std::cout << nbp + i * diff << ";" << FluidProperties(fluid.flash("Tx", nbp + i * diff, 0.5))[0].Pressure << std::endl;
    }


//
//    std::cout << "Propane at 25 C and 2 bar: " << std::endl;
//    auto a = FluidProperties(fluid.flash("PT", 2e5, 298.15));
////    std::cout << a << std::endl;
//    a.print(std::cout);
//    std::cout << "==============================================================================" << std::endl;
//
//    std::cout << "Compression to 10 bar: " << std::endl;
//    auto b = FluidProperties(fluid.flash("PS", (10E5), (a[0].Entropy)));
////    std::cout << b << std::endl;
//    b.print(std::cout);
//    std::cout << "==============================================================================" << std::endl;
//
//    std::cout << "Cooling to 25 C: " << std::endl;
//    auto c = FluidProperties(fluid.flash("PT", (10E5), (298.15)));
////    std::cout << c << std::endl;
//    c.print(std::cout);
//    std::cout << "==============================================================================" << std::endl;
//
//    std::cout << "Throttling to 2 bar: " << std::endl;
//    auto d = FluidProperties(fluid.flash("PH", (2E5), (c[0].Enthalpy)));
////    std::cout << d << std::endl;
//    d.print(std::cout);
//    std::cout << "==============================================================================" << std::endl;

    return 0;
}