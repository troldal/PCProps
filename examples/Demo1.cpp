#include <iomanip>
#include <iostream>

#include <EOSLib.hpp>
#include <Fluid.hpp>
#include <PropertyLib.hpp>
#include <PureComponent.hpp>
#include <PureComponentFactory.hpp>

#include <PhaseProperties.hpp>
#include <FluidProperties.hpp>
#include <json/json.hpp>
#include <OpenXLSX.hpp>

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
using namespace OpenXLSX;
using namespace nlohmann;

int main()
{

    XLDocument doc;
    doc.open("./Mini PCD.xlsx");
    auto wbk = doc.workbook();
    auto wks = wbk.worksheet("Constants");

    std::vector<json> objects;
    for (const auto &row: wks.rows(2, wks.rowCount())) {
        std::vector<XLCellValue> values = row.values();
        json object;
        object["Name"] = values[1].type() == XLValueType::Empty ? "" : values[1];
        object["Formula"] = values[3].type() == XLValueType::Empty ? "" : values[3];
        object["CAS"] = values[5].type() == XLValueType::Empty ? "" : values[5];
        object["MolarWeight"] = values[7].type() == XLValueType::Empty ? 0.0 : values[7].get<double>();
        object["MeltingTemperature"] = values[16].type() == XLValueType::Empty ? 0.0 : values[16].get<double>();
        object["BoilingTemperature"] = values[8].type() == XLValueType::Empty ? 0.0 : values[8].get<double>();
        object["CriticalTemperature"] = values[9].type() == XLValueType::Empty ? 0.0 : values[9].get<double>();
        object["CriticalPressure"] = values[10].type() == XLValueType::Empty ? 0.0 : values[10].get<double>();
        object["CriticalDensity"] = values[12].type() == XLValueType::Empty ? 0.0 : values[12].get<double>();
        object["AcentricFactor"] = values[14].type() == XLValueType::Empty ? 0.0 : values[14].get<double>();

        objects.emplace_back(object);
    }

    auto js = json(objects);
    auto pcfactory = PureComponentFactory(js.dump());
    auto obj = pcfactory.makeComponent("132259-10-0");

    std::cout << obj.property("MolarWeight") << std::endl;

//    data.name  = "PROPANE";
//    data.name  = "C3H8";
//    data.casrn = "74-98-6";

    auto pc = PureComponent{};

    using PCProps::LiquidVolume::Thomson;
    using namespace PCProps::CompressedLiquidViscosity;
    using namespace PCProps::CompressedVaporViscosity;

    pc.addDataItem("MolarWeight", 44.0956);
    pc.addDataItem("BoilingTemperature", 231.05);
    pc.addDataItem("FreezingTemperature", 85.15);
    pc.addDataItem("CriticalTemperature", 369.83);
    pc.addDataItem("CriticalPressure", 4.248E6);
    pc.addDataItem("CriticalVolume", 0.0002);
    pc.addDataItem("CriticalDensity", 220.48);
    pc.addDataItem("CriticalCompressibility", 0.2763);
    pc.addDataItem("AcentricFactor", 0.1523);
    pc.addDataItem("DipoleMoment", 0.083);

    pc.addDataItem("IdealGasCp", AlyLee(AlyLee::CreateFromDIPPR { 0.5192E5, 1.9245E5, 1.6265E3, 1.168E5, 723.6 }));
    pc.addDataItem("LiquidCp", PPDSLiquid(PPDSLiquid::CreateFromDIPPR { 62.983, 113630, 633.21, -873.46, 369.83 }));
    pc.addDataItem("VaporPressure", AntoineExtended(AntoineExtended::CreateFromDIPPR { 59.078, -3492.6, -6.0669, 1.0919E-05, 2 }));
    pc.addDataItem("SaturatedVaporViscosity", Viscosity::Lucas(369.83, 4.248E6, 0.2763, 44.096, 0.083));
    pc.addDataItem("SaturatedLiquidViscosity", KirchhoffExtended(-17.156, 646.25, 1.1101, -7.3439E-11, 4));
    pc.addDataItem("SaturatedLiquidVolume", Rackett(Rackett::CreateFromDIPPR { 1.3757, 0.27453, 369.83, 0.29359 }));

    pc.addDataItem("CompressedLiquidVolume", Thomson(369.83,4.248E6,0.1523));
    pc.addDataItem("CompressedLiquidViscosity", CompressedLiquidViscosity::Lucas(369.83,4.248E6,0.1523));
    pc.addDataItem("CompressedVaporViscosity", CompressedVaporViscosity::Lucas(369.83,4.248E6,0.2763,44.096,0.083));

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