#include <CDJoback.hpp>
#include <PCComponent.hpp>
#include <PCPropsException.hpp>
#include <SLVRackett.hpp>
#include <VPHoffmannFlorin.hpp>
#include <iostream>
#include <list>

using PCProps::ConstantData::CDJoback;
using PCProps::ConstantData::CDJobackGroup;

using PCProps::LiquidVolumes::SLVRackett;

int main()
{
    auto rackett = SLVRackett();
    std::cout << rackett(300.0) << std::endl;

    auto R143a_1 = SLVRackett::createFromCriticalProperties(346.30, 37.92E5, 0.255);
    std::cout << "R143a Density @ 300K: " << R143a_1(300.0) * 1000000 << std::endl;

    auto R143a_2 = SLVRackett::createFromAcentricFactor(346.30, 37.92E5, 0.259);
    std::cout << "R143a Density @ 300K: " << R143a_2(300.0) * 1000000 << std::endl;

    auto R143a_3 = SLVRackett::createFromReferencePointA(346.3, 245.0, 75.38E-6, 0.259);
    std::cout << "R143a Density @ 300K: " << R143a_3(300.0) * 1000000 << std::endl;

    auto R143a_4 = SLVRackett::createFromReferencePointB(346.3, 245.0, 75.38E-6, 0.255);
    std::cout << "R143a Density @ 300K: " << R143a_4(300.0) * 1000000 << std::endl;

    CDJoback acetone(std::list<CDJobackGroup> { { 1, 2 }, { 24, 1 } }, 58.08, 10);
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

    return 0;
}