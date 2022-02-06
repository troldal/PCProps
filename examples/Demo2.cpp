#include <iomanip>
#include <iostream>
#include <vector>
#include <variant>

#include <DataSource.hpp>
#include <EquationOfState.hpp>
#include <FluidProperties.hpp>
#include <PropertyPackage.hpp>
#include <PureComponentFactory.hpp>
#include <UnitOps.hpp>
#include <IUnitOperation.hpp>

#include <sciplot/sciplot.hpp>

using PCProps::EquationOfState::PengRobinson;
using namespace PCProps;
using namespace sciplot;

int main()
{
    //    auto ds  = DataSource("Mini PCD.xlsx");
    //    auto pcf = PureComponentFactory(ds.load());
    //
    //    auto water           = pcf.makeComponent("7732-18-5");
    //    auto propertypackage = PropertyPackage(water, PengRobinson {});

    //    auto inletStream = UnitOps::Stream(propertypackage, 5.0);
    //
    //    std::cout << "Water at 25 C and 5 bar: " << std::endl;
    //    inletStream.flash("PT", 5E5, 298.15);
    //    FluidProperties(inletStream.properties()).print(std::cout);
    //    std::cout << "===============================================================" << std::endl;
    //
    //    auto pump = UnitOps::CentrifugalPump("{\"OutletPressure\": 10E5, \"PumpEfficiency\": 0.7 }");
    //    pump.setInletStreams({ &inletStream });
    //
    //    auto& outletStream = pump.outletStreams().front();
    //    pump.compute();
    //
    //    std::cout << "Water at 25 C and 10 bar: " << std::endl;
    //    FluidProperties(outletStream.properties()).print(std::cout);
    //    std::cout << "===============================================================" << std::endl;

    //    std::cout << pump.results() << std::endl;


//    using Element = std::variant<IUnitOperation, UnitOps::Stream>;
//    std::vector<Element> unitOps;
//
//
    auto ds  = DataSource("Mini PCD.xlsx");
    auto pcf = PureComponentFactory(ds.load());

    auto propane           = pcf.makeComponent("74-98-6");
    auto propertypackage = PropertyPackage(propane, PengRobinson {});
//
//    unitOps.emplace_back(UnitOps::Stream(propertypackage, 5.0));
//    unitOps.back() setSpecification("{\"FlashSpecification\": \"PT\", \"Pressure\": 10E5, \"Temperature\": 298.15}");
//
//    auto valve = UnitOps::Valve();
//    valve.setInletStreams({ unitOps.back() });
//    unitOps.emplace_back(UnitOps::Valve());
//    unitOps.back().setInletStreams({ &inletStream });


    auto inletStream = UnitOps::Stream(propertypackage, 5.0);

    std::cout << "Propane at 25 C and 10 bar: " << std::endl;
//    inletStream.flash("PT", 10E5, 298.15);
    inletStream.setSpecification("{\"FlashSpecification\": \"PT\", \"Pressure\": 10E5, \"Temperature\": 298.15}");
    inletStream.compute();
    FluidProperties(inletStream.results()).print(std::cout);
    std::cout << "===============================================================" << std::endl;

//    auto valve = UnitOps::Valve("{\"OutletPressure\": 2E5}");
    auto valve = UnitOps::Valve();
    valve.setInletStreams({ &inletStream });

    auto& outletStream = valve.outletStreams().front();
    valve.setSpecification("{\"OutletPressure\": 2E5}");
//    unitOps.emplace_back(std::move(valve));
//    for (auto& uo : unitOps) uo.compute();
    valve.compute();

    std::cout << "Water at X C and 2 bar: " << std::endl;
    FluidProperties(outletStream.results()).print(std::cout);
    std::cout << "===============================================================" << std::endl;

    return 0;
}