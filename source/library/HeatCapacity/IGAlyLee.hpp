//
// Created by Kenneth Balslev on 31/12/2020.
//

#ifndef PCPROPS_IGALYLEE_HPP
#define PCPROPS_IGALYLEE_HPP

namespace PCProps::HeatCapacity
{
    class IGAlyLee
    {
        double m_A { 0.0 };
        double m_B { 0.0 };
        double m_C { 0.0 };
        double m_D { 0.0 };
        double m_E { 0.0 };

    public:
        IGAlyLee();

        IGAlyLee(double coeffA, double coeffB, double coeffC, double coeffD, double coeffE);

        IGAlyLee(const IGAlyLee& other);

        IGAlyLee(IGAlyLee&& other) noexcept;

        ~IGAlyLee();

        IGAlyLee& operator=(const IGAlyLee& other);

        IGAlyLee& operator=(IGAlyLee&& other) noexcept;

        double operator()(double temperature);

        double evaluate(double temperature);

        double differential(double temperature);

        double integral(double temperature);
    };

}    // namespace PCProps::HeatCapacity

#endif    // PCPROPS_IGALYLEE_HPP
