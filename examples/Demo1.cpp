#include <iomanip>
#include <iostream>

#include <EOSLib.hpp>
#include <Fluid.hpp>
#include <PureComponentFactory.hpp>
#include <DataSource.hpp>
#include <FluidProperties.hpp>

#include <sciplot/sciplot.hpp>

using PCProps::EquationOfState::PengRobinson;
using namespace PCProps;
using namespace sciplot;

int main()
{

    auto ds = DataSource("Mini PCD.xlsx");
    auto pcf = PureComponentFactory(ds.load());

    auto plot_psat = [&](const std::string& CAS, const std::string& Name) {
        std::cout << "Plotting Vapor Pressure Curve: " << Name << std::endl;
        auto pc = pcf.makeComponent(CAS);
        auto fluid = Fluid(pc, PengRobinson{});

        auto nbp = pc.property("NormalFreezingPoint");
        auto tc = pc.property("CriticalTemperature");
        auto diff = (tc - nbp) / 99.0;
        std::vector<double> x, y;
        for (int i = 0; i < 100; ++i) {
            x.emplace_back(nbp + i * diff);
            y.emplace_back(FluidProperties(fluid.flash("Tx", nbp + i * diff, 0.5))[0].Pressure/100000);
        }

        Plot plot;
        plot.fontName("Avenir");
        plot.autoclean(true);
        plot.gnuplot("set title \"Vapor Pressure\" font \"Avenir\"");

        // Plot the data
        plot.drawCurve(x, y).label(Name);
        plot.legend().atTopLeft();
        plot.grid().show();
        plot.ylabel("Pressure (Pa)");
        plot.xlabel("Temperature (K)");

        plot.save("./Psat (T) - " + Name + ".pdf");

    };

    auto plot_psat2 = [&](const std::string& CAS, const std::string& Name) {
        std::cout << "Plotting Vapor Pressure Curve: " << Name << std::endl;
        auto pc = pcf.makeComponent(CAS);
        auto fluid = Fluid(pc, PengRobinson{});

        auto cp = pc.property("CriticalPressure");
        auto diff = (cp - 100.0) / 99.0;
        std::vector<double> x, y;
        for (int i = 0; i < 100; ++i) {
            x.emplace_back(FluidProperties(fluid.flash("Px", 100.0 + i * diff, 0.5))[0].Temperature);
            y.emplace_back((100.0 + i * diff)/100000);
        }

        Plot plot;
        plot.fontName("Avenir");
        plot.autoclean(true);
        plot.gnuplot("set title \"Vapor Pressure\" font \"Avenir\"");

        // Plot the data
        plot.drawCurve(x, y).label(Name);
        plot.legend().atTopLeft();
        plot.grid().show();
        plot.ylabel("Pressure (Pa)");
        plot.xlabel("Temperature (K)");

        plot.save("./Psat (P) - " + Name + ".pdf");

    };

    auto plot_pv = [&](const std::string& CAS, const std::string& Name) {
        std::cout << "Plotting PV Isotherm: " << Name << std::endl;
        auto pcomp = pcf.makeComponent(CAS);
        auto fluid = Fluid(pcomp, PengRobinson{});
        auto pc = pcomp.property("CriticalPressure");
        auto tc = pcomp.property("CriticalTemperature");

        auto b = (0.07779607 * 8.31446261815324 * tc / pc);

        double temp = tc;
        auto f = 1.2;
        auto diff = (b * 50 - b*f) / 299.0;
        auto psat = FluidProperties(fluid.flash("Tx", temp, 0.5))[0].Pressure/100000;

        Plot plot;
        plot.fontName("Avenir");
        plot.autoclean(true);
        plot.gnuplot("set title \"" +Name + " PV Isotherm\" font \"Avenir\"");

        for (int j = 0; j < 2; ++j) {
            std::vector<double> x, y;
            for (int i = 0; i < 300; ++i) {
                temp = tc - j * 20.0;
                auto p = FluidProperties(fluid.flash("TV", temp, b * f + i * diff))[0].Pressure / 100000;
                if (p > psat * 3) continue;
                x.emplace_back(b * f + i * diff);
                y.emplace_back(p);
            }
            plot.drawCurve(x, y).labelNone();

        }

        //plot.legend().atOutsideBottom();
        plot.grid().show();
        plot.ylabel("Pressure (bar)");
        plot.xlabel("Volume (m3/mol)");

        plot.save("./PV Isotherm - " + Name + ".pdf");

    };

    plot_pv("132259-10-0", "AIR");
    plot_pv("7664-41-7", "AMMONIA");
    plot_pv("7440-37-1", "ARGON");
    plot_pv("106-97-8", "BUTANE");
    plot_pv("124-38-9", "CARBON DIOXIDE");
    plot_pv("630-08-0", "CARBON MONOXIDE");
    plot_pv("124-18-5", "DECANE");
    plot_pv("74-84-0", "ETHANE");
    plot_pv("50-00-0", "FORMALDEHYDE");
    plot_pv("142-82-5", "HEPTANE");
    plot_pv("110-54-3", "HEXANE");
    plot_pv("7783-06-4", "HYDROGEN SULFIDE");
    plot_pv("74-82-8", "METHANE");
    plot_pv("7727-37-9", "NITROGEN");
    plot_pv("10102-43-9", "NITRIC OXIDE");
    plot_pv("111-84-2", "NONANE");
    plot_pv("111-65-9", "OCTANE");
    plot_pv("7782-44-7", "OXYGEN");
    plot_pv("109-66-0", "PENTANE");
    plot_pv("74-98-6", "PROPANE");
    plot_pv("7446-09-5", "SULFUR DIOXIDE");
    plot_pv("7732-18-5", "WATER");

    plot_psat("132259-10-0", "AIR");
    plot_psat("7664-41-7", "AMMONIA");
    plot_psat("7440-37-1", "ARGON");
    plot_psat("106-97-8", "BUTANE");
    plot_psat("124-38-9", "CARBON DIOXIDE");
    plot_psat("630-08-0", "CARBON MONOXIDE");
    plot_psat("124-18-5", "DECANE");
    plot_psat("74-84-0", "ETHANE");
    plot_psat("50-00-0", "FORMALDEHYDE");
    plot_psat("142-82-5", "HEPTANE");
    plot_psat("110-54-3", "HEXANE");
    plot_psat("7783-06-4", "HYDROGEN SULFIDE");
    plot_psat("74-82-8", "METHANE");
    plot_psat("7727-37-9", "NITROGEN");
    plot_psat("10102-43-9", "NITRIC OXIDE");
    plot_psat("111-84-2", "NONANE");
    plot_psat("111-65-9", "OCTANE");
    plot_psat("7782-44-7", "OXYGEN");
    plot_psat("109-66-0", "PENTANE");
    plot_psat("74-98-6", "PROPANE");
    plot_psat("7446-09-5", "SULFUR DIOXIDE");
    plot_psat("7732-18-5", "WATER");

    plot_psat2("132259-10-0", "AIR");
    plot_psat2("7664-41-7", "AMMONIA");
    plot_psat2("7440-37-1", "ARGON");
    plot_psat2("106-97-8", "BUTANE");
    plot_psat2("124-38-9", "CARBON DIOXIDE");
    plot_psat2("630-08-0", "CARBON MONOXIDE");
    plot_psat2("124-18-5", "DECANE");
    plot_psat2("74-84-0", "ETHANE");
    plot_psat2("50-00-0", "FORMALDEHYDE");
    plot_psat2("142-82-5", "HEPTANE");
    plot_psat2("110-54-3", "HEXANE");
    plot_psat2("7783-06-4", "HYDROGEN SULFIDE");
    plot_psat2("74-82-8", "METHANE");
    plot_psat2("7727-37-9", "NITROGEN");
    plot_psat2("10102-43-9", "NITRIC OXIDE");
    plot_psat2("111-84-2", "NONANE");
    plot_psat2("111-65-9", "OCTANE");
    plot_psat2("7782-44-7", "OXYGEN");
    plot_psat2("109-66-0", "PENTANE");
    plot_psat2("74-98-6", "PROPANE");
    plot_psat2("7446-09-5", "SULFUR DIOXIDE");
    plot_psat2("7732-18-5", "WATER");

    //    FluidProperties(fluid.flash("Tx", 447.3, 0.5)).print(std::cout);
//    FluidProperties(fluid.flash("Px", 4247999.0, 0.5)).print(std::cout);


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