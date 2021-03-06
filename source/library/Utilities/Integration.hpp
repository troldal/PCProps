//
// Created by Kenneth Balslev on 01/01/2021.
//

#ifndef PCPROPS_INTEGRATION_HPP
#define PCPROPS_INTEGRATION_HPP

namespace PCProps::Numerics
{
    double integrate(const std::function<double(double)>& func, double x1, double x2)
    {
        using std::abs;
        using std::min;
        auto trapez = [&](double first, double second) -> double {
            auto val1 = func(first);
            auto val2 = func(second);

            auto area1 = (second - first) * val1;
            auto area2 = (second - first) * val2;

            auto diff = abs(area1 - area2);

            return min(area1, area2) + diff / 2.0;
        };

        auto diff      = x2 - x1;
        auto intervals = 100;
        auto result    = 0.0;

        for (int i = 0; i < intervals; ++i) {
            result += trapez(x1 + i * (diff / intervals), x1 + (i + 1) * (diff / intervals));
        }

        return result;
    }

}    // namespace PCProps::Numerics

#endif    // PCPROPS_INTEGRATION_HPP
