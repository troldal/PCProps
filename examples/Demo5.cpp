#include <iomanip>
#include <iostream>

#include <external/numeric/integration.hpp>

int main()
{
    std::cout << numeric::integrate([&](double x) { return x * x; }, -1.0, 1.0) << std::endl;

    return 0;
}