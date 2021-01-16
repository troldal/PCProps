#include <iomanip>
#include <iostream>

#include <library/Utilities/Calculus.hpp>

using namespace PCProps::Numerics;

int main()
{
    std::cout << std::setprecision(20);
    std::cout << integrate([&](double x) { return x * x ; }, 0.0, 1.0) << std::endl;


    return 0;
}