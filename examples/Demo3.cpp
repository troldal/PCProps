#include <iomanip>
#include <iostream>

#include <library/Viscosity/KirchhoffExtended.hpp>
#include <library/Viscosity/LucasHPL.hpp>
#include <library/VaporPressure/AmbroseWalton.hpp>

using PCProps::Viscosity::LucasHPL;
using PCProps::Viscosity::KirchhoffExtended;
using PCProps::VaporPressure::AmbroseWalton;
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