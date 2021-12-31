#include <iomanip>
#include <iostream>

#include <PropertyPackage.hpp>
#include <Fluid.hpp>
#include <PropertyLib.hpp>

#include <common/PhaseProperties.hpp>

using PCProps::VaporPressure::AmbroseWalton;

using PCProps::EquationOfState::PengRobinson;
using PCProps::HeatCapacity::AlyLee;
using PCProps::HeatCapacity::PPDSLiquid;
using PCProps::LiquidVolume::Rackett;
using PCProps::VaporPressure::AntoineExtended;

using PCProps::PCComponentData;
using PCProps::Viscosity::Lucas;
using PCProps::Viscosity::DIPPR102;
using PCProps::Viscosity::KirchhoffExtended;

using PCProps::Enthalpy;
using PCProps::Entropy;
using PCProps::Pressure;
using PCProps::Temperature;

using namespace PCProps;

int main()
{

    auto mu_lp = KirchhoffExtended(-11.358, 1213.1, 0.0, 0.0, 0.0);
    auto psat = AmbroseWalton(572.19, 34.71E5, 0.235);
    auto mu_hp = LucasHPL(572.19, 34.71E5, 0.235, psat, mu_lp);
    std::cout << mu_lp(146.58) << std::endl;
    std::cout << mu_lp(457.68) << std::endl << std::endl;

    std::cout << mu_hp(300.0, 500E5) << std::endl;

    return 0;
}