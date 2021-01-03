#include <ConstantData/CDJoback.hpp>
#include <LiquidVolume/SLVRackett.hpp>
#include <VaporPressure/VPAmbroseWalton.hpp>
#include <VaporPressure/VPAntoineExt.hpp>
#include <VaporPressure/VPHoffmannFlorin.hpp>
#include <VaporPressure/VPRiedel.hpp>
#include <iomanip>
#include <iostream>
#include <library/EquationOfState/EOSPengRobinson.hpp>
#include <library/EquationOfState/EOSUtilities.hpp>
#include <library/HeatCapacity/IGAlyLee.hpp>
#include <library/LiquidVolume/CLVAalto.hpp>
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

using PCProps::HeatCapacity::IGAlyLee;
using PCProps::LiquidVolume::SLVElbro;
using PCProps::LiquidVolume::SLVElbroGroup;
using PCProps::LiquidVolume::SLVHankinsonThomson;
using PCProps::LiquidVolume::SLVRackett;
using PCProps::LiquidVolume::SLVYenWoods;

using namespace PCProps::VaporPressure;
using namespace PCProps::LiquidVolume;

using namespace PCProps::EquationOfState;

void print(const PhaseData& data)
{
    std::cout << std::setprecision(6) << std::fixed;
    std::cout << "Temperature      : " << std::right << std::setw(15) << std::get<Temperature>(data) << " K" << std::endl;
    std::cout << "Pressure         : " << std::right << std::setw(15) << std::get<Pressure>(data) << " Pa" << std::endl;
    std::cout << "Compressibility  : " << std::right << std::setw(15) << std::get<Compressibility>(data) << std::endl;
    std::cout << "Fugacity         : " << std::right << std::setw(15) << std::get<Fugacity>(data) << " Pa" << std::endl;
    std::cout << "Moles            : " << std::right << std::setw(15) << std::get<Moles>(data) << std::endl;
    std::cout << "Molecular Weight : " << std::right << std::setw(15) << std::get<MolecularWeight>(data) << " g/mol" << std::endl;
    std::cout << "Molar Volume     : " << std::right << std::setw(15) << std::get<Volume>(data) << " m3/mol" << std::endl;
    std::cout << "Enthalpy         : " << std::right << std::setw(15) << std::get<Enthalpy>(data) << " J/mol" << std::endl;
    std::cout << "Entropy          : " << std::right << std::setw(15) << std::get<Entropy>(data) << " J/mol-K" << std::endl;
    std::cout << "Internal Energy  : " << std::right << std::setw(15) << std::get<InternalEnergy>(data) << " J/mol" << std::endl;
    std::cout << "Gibbs Energy     : " << std::right << std::setw(15) << std::get<GibbsEnergy>(data) << " J/mol" << std::endl;
    std::cout << "Helmholz Energy  : " << std::right << std::setw(15) << std::get<HelmholzEnergy>(data) << " J/mol" << std::endl;
    std::cout << std::endl;
}

int main()
{
    //            auto            PSat = VPAmbroseWalton(190.6, 4.604E6, 0.011);
    //            auto            igCp = IGAlyLee(0.33298E5, 0.79933E5, 2.0869E3, 0.41602E5, 991.69);
    //            EOSPengRobinson methane(190.6, 4.604E6, 0.011, 16.043, PSat, igCp);
    //            double          t = 112.4; // K
    //            double          p = 0.1E6; // Pa
    //
    //            for (const auto& phase : methane.flashPT(p, t))
    //                print(phase);

    auto tc    = 369.83;
    auto pc    = 4.248E6;
    auto omega = 0.1523;
    auto mw    = 44.096;

    auto            PSat = VPAmbroseWalton(tc, pc, omega);
    auto            igCp = IGAlyLee(0.5192E5, 1.9245E5, 1.6265E3, 1.168E5, 723.6);
    EOSPengRobinson propane(tc, pc, omega, mw, PSat, igCp);

    for (const auto& phase : propane.flashPS(5000000, 3)) print(phase);

    //    for (const auto& phase : propane.flashPT(101325 * 10, 298.15)) print(phase);
    //    for (const auto& phase : propane.flashPH(101325 * 2, -16139.662736)) print(phase);

    return 0;
}