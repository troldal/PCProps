#include <iomanip>
#include <iostream>

#include <EquationOfState.hpp>
#include <PropertyPackage.hpp>
#include <PureComponentFactory.hpp>
#include <DataSource.hpp>
#include <FluidProperties.hpp>
#include <UnitOps.hpp>

#include <sciplot/sciplot.hpp>

using PCProps::EquationOfState::PengRobinson;
using namespace PCProps;
using namespace sciplot;

int main()
{

    auto ds = DataSource("Mini PCD.xlsx");
    auto pcf = PureComponentFactory(ds.load());

    auto propane    = pcf.makeComponent("74-98-6");
    auto propertypackage = PropertyPackage(propane, PengRobinson {});

    auto fluid   = UnitOps::Stream(propertypackage, 10.0);


    std::cout << "Water at 25 C and 5 bar: " << std::endl;
    auto a = fluid.flashPT(5E5, 298.15);
    for (const auto& phase : a) std::cout << phase << std::endl;
    std::cout << "==================================================" << std::endl;

    auto cp = UnitOps::CentrifugalPump();
    cp.setDifferentialPressure(20E5);
    cp.setInletStream(&fluid);
    auto b = cp().properties();

    //    std::cout << "Water at 25 C and 25 bar: " << std::endl;
    //    auto b = fluid.flashPT(25E6, 298.15);
    for (const auto& phase : b) std::cout << phase << std::endl;
    std::cout << "==================================================" << std::endl;

    //    auto pipe = Pipe(9000, 152.0 / 1000, 7, 0.0);
    //    std::cout << pipe.computeInletPressure(fluid, 2600) << std::endl;

    return 0;
}