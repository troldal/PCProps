#include <iomanip>
#include <iostream>

#include <DataSource.hpp>
#include <EquationOfState.hpp>
#include <FluidProperties.hpp>
#include <PropertyPackage.hpp>
#include <PureComponentFactory.hpp>
#include <UnitOps.hpp>

#include <sciplot/sciplot.hpp>

using PCProps::EquationOfState::PengRobinson;
using namespace PCProps;
using namespace sciplot;

int main()
{
    auto ds  = DataSource("Mini PCD.xlsx");
    auto pcf = PureComponentFactory(ds.load());

    auto water           = pcf.makeComponent("7732-18-5");
    auto propertypackage = PropertyPackage(water, PengRobinson {});

    auto inletStream = UnitOps::Stream(propertypackage, 5.0);

    std::cout << "Water at 25 C and 5 bar: " << std::endl;
    inletStream.flash("PT", 5E5, 298.15);
    FluidProperties(inletStream.properties()).print(std::cout);
    std::cout << "===============================================================" << std::endl;

    auto pump = UnitOps::CentrifugalPump("{\"OutletPressure\": 10E5, \"PumpEfficiency\": 0.75}");
    pump.setInletStream(&inletStream);

    auto& outletStream = pump.outputStream();
    pump.compute();

    std::cout << "Water at 25 C and 10 bar: " << std::endl;
    FluidProperties(outletStream.properties()).print(std::cout);
    std::cout << "===============================================================" << std::endl;

    std::cout << pump.results() << std::endl;

    return 0;
}