#include <VPRiedel.hpp>
#include <PCComponent.hpp>
#include <PCPropsException.hpp>
#include <iostream>

int main()
{
    PCProps::PCComponentData data;
    data.vaporPressureFunction = PCProps::VaporPressure::VPRiedel(273.15+20.3, 466.05, 55.5E5);
    data.molecularWeight = 16.043;

    try {
        PCProps::PCComponent comp(data);
        std::cout << "Psat: " << comp.vaporPressure(273.15+19.9) << std::endl;
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