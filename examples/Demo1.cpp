#include <ConstantData/CDJoback.hpp>
#include <LiquidVolume/Rackett.hpp>
#include <VaporPressure/AmbroseWalton.hpp>
#include <VaporPressure/AntoineExtended.hpp>
#include <VaporPressure/HoffmannFlorin.hpp>
#include <VaporPressure/Riedel.hpp>
#include <iomanip>
#include <iostream>
#include <library/EquationOfState/EOSPengRobinson.hpp>
#include <library/EquationOfState/EOSUtilities.hpp>
#include <library/HeatCapacity/AlyLee.hpp>
#include <library/LiquidVolume/Aalto.hpp>
#include <library/LiquidVolume/Elbro.hpp>
#include <library/LiquidVolume/HankinsonThomson.hpp>
#include <library/LiquidVolume/Thomson.hpp>
#include <library/LiquidVolume/YenWoods.hpp>
#include <library/PCComponent.hpp>
#include <library/PCPropsException.hpp>
#include <list>
#include <map>
#include <vector>

using PCProps::ConstantData::CDJoback;
using PCProps::ConstantData::CDJobackGroup;

using PCProps::HeatCapacity::AlyLee;
using PCProps::LiquidVolume::Elbro;
using PCProps::LiquidVolume::SLVElbroGroup;

using namespace PCProps::VaporPressure;
using namespace PCProps::LiquidVolume;

using namespace PCProps::EquationOfState;

void print(const PhaseData& data)
{
    std::cout << std::setprecision(10) << std::fixed;
    std::cout << "Temperature      : " << std::right << std::setw(20) << std::get<Temperature>(data) << " K" << std::endl;
    std::cout << "Pressure         : " << std::right << std::setw(20) << std::get<Pressure>(data) << " Pa" << std::endl;
    std::cout << "Compressibility  : " << std::right << std::setw(20) << std::get<Compressibility>(data) << std::endl;
    std::cout << "Fugacity         : " << std::right << std::setw(20) << std::get<Fugacity>(data) << " Pa" << std::endl;
    std::cout << "Moles            : " << std::right << std::setw(20) << std::get<MolarFraction>(data) << std::endl;
    std::cout << "Molecular Weight : " << std::right << std::setw(20) << std::get<MolecularWeight>(data) << " g/mol" << std::endl;
    std::cout << "Molar Volume     : " << std::right << std::setw(20) << std::get<Volume>(data) << " m3/mol" << std::endl;
    std::cout << "Enthalpy         : " << std::right << std::setw(20) << std::get<Enthalpy>(data) << " J/mol" << std::endl;
    std::cout << "Entropy          : " << std::right << std::setw(20) << std::get<Entropy>(data) << " J/mol-K" << std::endl;
    std::cout << "Internal Energy  : " << std::right << std::setw(20) << std::get<InternalEnergy>(data) << " J/mol" << std::endl;
    std::cout << "Gibbs Energy     : " << std::right << std::setw(20) << std::get<GibbsEnergy>(data) << " J/mol" << std::endl;
    std::cout << "Helmholz Energy  : " << std::right << std::setw(20) << std::get<HelmholzEnergy>(data) << " J/mol" << std::endl;
    std::cout << std::endl;
}

int main()
{
    using std::get;

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

    auto vp = AntoineExtended(AntoineExtended::CreateFromDIPPR { 59.078, -3492.6, -6.0669, 1.0919E-05, 2 });
    std::cout << vp(85.47) << std::endl;
    std::cout << vp(369.83) << std::endl;

    auto            PSat = AmbroseWalton(tc, pc, omega);
    auto            igCp = AlyLee(AlyLee::CreateFromDIPPR { 0.5192E5, 1.9245E5, 1.6265E3, 1.168E5, 723.6 });
    EOSPengRobinson propane(tc, pc, omega, mw, vp, igCp);

    //    for (const auto& phase : propane.flashPS(5000000, 3)) print(phase);

    //    for (const auto& phase : propane.flashPT(101325 * 10, 298.15)) print(phase);
    //    for (const auto& phase : propane.flashPH(101325 * 2, -16139.662736)) print(phase);

    std::cout << "Propane at 25 C and 2 bar: " << std::endl;
    auto a = propane.flashPT(2E5, 298.15);
    for (const auto& phase : a) print(phase);

    std::cout << "Compression to 10 bar: " << std::endl;
    auto b = propane.flashPS(10E5, get<Entropy>(a[0]));
    for (const auto& phase : b) print(phase);

    std::cout << "Cooling to 25 C: " << std::endl;
    auto c = propane.flashPT(10E5, 298.15);
    std::cout << "dT: " << get<Temperature>(c[0]) - get<Temperature>(b[0]) << std::endl;
    std::cout << "dH: " << get<Enthalpy>(c[0]) - get<Enthalpy>(b[0]) << std::endl;
    for (const auto& phase : c) print(phase);

    std::cout << "Throttling to 2 bar: " << std::endl;
    auto d = propane.flashPH(2E5, get<Enthalpy>(c[0]));
    for (const auto& phase : d) print(phase);

    return 0;
}