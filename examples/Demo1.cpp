#include <CDJoback.hpp>
#include <PCComponent.hpp>
#include <PCPropsException.hpp>
#include <VPAmbroseWalton.hpp>
#include <VPHoffmannFlorin.hpp>
#include <VPRiedel.hpp>
#include <iostream>
#include <vector>

int main()
{
    PCProps::ConstantData::CDJoback acetone(std::vector<std::pair<int, int>> { std::make_pair(2, 1), std::make_pair(1, 24) }, 58.08, 10);

    //    PCProps::ConstantData::CDJoback acetone(groups, 58.08, 10);
    std::cout << "Acetone Tb: " << acetone.boilingTemperature() << std::endl;
    std::cout << "Acetone Tm: " << acetone.meltingTemperature() << std::endl;
    std::cout << "Acetone Tc: " << acetone.criticalTemperature() << std::endl;
    std::cout << "Acetone Pc: " << acetone.criticalPressure() << std::endl;
    std::cout << "Acetone Vc: " << acetone.criticalVolume() << std::endl;
    std::cout << "Acetone Hform: " << acetone.enthalpyOfFormation() << std::endl;
    std::cout << "Acetone Gform: " << acetone.gibbsEnergyOfFormation() << std::endl;
    std::cout << "Acetone Hfus: " << acetone.enthalpyOfFusion() << std::endl;
    std::cout << "Acetone Hvap: " << acetone.enthalpyOfVaporization() << std::endl;
    std::cout << "Acetone igCp: " << acetone.idealGasCp(300.0) << std::endl;
    std::cout << "Acetone liquid viscosity: " << acetone.liquidViscosity(300.0) << std::endl;

    PCProps::PCComponentData data;
    data.vaporPressureFunction = PCProps::VaporPressure::VPHoffmannFlorin(404.87, 101325.0, 632.35, 45.1911E5);
    data.molecularWeight       = 16.043;

    try {
        PCProps::PCComponent comp(data);
        std::cout << "Psat: " << comp.vaporPressure(500.0) << std::endl;
        std::cout << "MW: " << comp.molecularWeight() << std::endl;
    }

    catch (const PCProps::PCPropsException& e) {
        std::cout << e.what();
    }


//    auto psat = PCProps::VaporPressure::VPRiedel();
//    std::cout << "Psat: " << psat(500.0) << std::endl;
//
//    psat = PCProps::VaporPressure::VPRiedel(404.87, 632.35, 45.1911E5);
//    std::cout << "Psat: " << psat(500.0) << std::endl;
//
//    psat = PCProps::VaporPressure::VPRiedel(R"({"tboil": 404.87, "tcrit": 632.35, "pcrit": 45.1911E5})");
//    std::cout << "Psat: " << psat(500.0) << std::endl;
//    std::cout << psat.coefficients() << std::endl;

    return 0;
}