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
    auto pc = pcf.makeComponent("74-98-6");
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