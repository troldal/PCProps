#include <ConstantData/CDJoback.hpp>
#include <LiquidVolume/SLVRackett.hpp>
#include <VaporPressure/VPHoffmannFlorin.hpp>
#include <VaporPressure/VPRiedel.hpp>
#include <iostream>
#include <library/LiquidVolume/CLVThomson.hpp>
#include <library/LiquidVolume/SLVElbro.hpp>
#include <library/LiquidVolume/SLVHankinsonThomson.hpp>
#include <library/LiquidVolume/SLVYenWoods.hpp>
#include <library/PCComponent.hpp>
#include <library/PCPropsException.hpp>
#include <list>
#include <map>
#include <vector>

using PCProps::ConstantData::CDJoback;
using PCProps::ConstantData::CDJobackGroup;

using PCProps::LiquidVolume::SLVElbro;
using PCProps::LiquidVolume::SLVElbroGroup;
using PCProps::LiquidVolume::SLVHankinsonThomson;
using PCProps::LiquidVolume::SLVRackett;
using PCProps::LiquidVolume::SLVYenWoods;

using namespace PCProps::VaporPressure;
using namespace PCProps::LiquidVolume;

int main()
{
    //    auto ammonia_psat = VPRiedel(239.706, 405.65, 11280000.000);
    //    auto ammonia_vsat = SLVHankinsonThomson::createFromCharacteristicVolume(405.65, 0.00007, 0.2526);

    auto ammonia_psat = [](double _) { return 102.97E5; };
    auto ammonia_vsat = [](double _) { return 49.15E-6; };
    auto ammonia      = CLVThomson::create(405.4, 113.53E5, 0.256, ammonia_vsat, ammonia_psat);
    std::cout << ammonia(398.0, 400E5) * 1000000 << std::endl;

    auto hexadecane = SLVElbro::create(std::vector<SLVElbroGroup> { { 1, 2 }, { 2, 14 } });
    std::cout << hexadecane(298.15) << std::endl;

    auto isobutane = SLVHankinsonThomson::createFromCharacteristicVolume(408.04, 256.8E-6, 0.1825);
    std::cout << isobutane(310.93) << std::endl;

    auto isobutane2 = SLVHankinsonThomson::createFromEstimatedProperties(408.04, 3640000, 0.1825);
    std::cout << isobutane2(310.93) << std::endl;

    auto yw = SLVYenWoods::createFromYenWoodsEstimation(647.14, 55.45E-6, 0.245);
    std::cout << yw(300.0) << std::endl;

    auto ppds = SLVYenWoods::createFromPPDSCoefficients(647.14, 55.9472E-6, 18.02, 1094.0233, -1813.2295, 3863.9557, -2479.813);
    std::cout << ppds(300.0) << std::endl;

    auto dippr = SLVYenWoods::createFromDIPPR116Coefficients(647.14, 55.9472E-6, 58.606, -95.396, 213.89, -141.26);
    std::cout << dippr(300.0) << std::endl;

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