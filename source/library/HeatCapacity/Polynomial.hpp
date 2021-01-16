//
// Created by Kenneth Balslev on 04/01/2021.
//

#ifndef PCPROPS_POLYNOMIAL_HPP
#define PCPROPS_POLYNOMIAL_HPP

namespace PCProps::HeatCapacity
{
    class Polynomial
    {
        double m_A {};
        double m_B {};
        double m_C {};
        double m_D {};
        double m_E {};
        double m_F {};
        double m_G {};

    public:
        struct CreateFromDIPPR
        {
            double A {};
            double B {};
            double C {};
            double D {};
            double E {};
        };

        struct CreateFromYaws
        {
            double A {};
            double B {};
            double C {};
            double D {};
        };

        Polynomial();

        Polynomial(const CreateFromDIPPR& coefficients);

        Polynomial(const CreateFromYaws& coefficients);

        Polynomial(const Polynomial& other);

        Polynomial(Polynomial&&) noexcept;

        ~Polynomial();

        Polynomial& operator=(const Polynomial& other);

        Polynomial& operator=(Polynomial&& other) noexcept;

        double operator()(double temperature) const;
    };
}    // namespace PCProps::HeatCapacity

#endif    // PCPROPS_POLYNOMIAL_HPP
