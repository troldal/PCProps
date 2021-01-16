#include <iomanip>
#include <iostream>

#include <library/Utilities/Calculus.hpp>
#include <external/numeric/integration.hpp>

using namespace PCProps::Numerics;

int main()
{
    std::cout << numeric::integrate([&](double x) { return x * x; }, -1.0, 1.0) << std::endl;

    return 0;
}