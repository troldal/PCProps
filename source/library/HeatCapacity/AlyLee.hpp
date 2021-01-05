//
// Created by Kenneth Balslev on 31/12/2020.
//

#ifndef PCPROPS_ALYLEE_HPP
#define PCPROPS_ALYLEE_HPP

namespace PCProps::HeatCapacity
{
    class AlyLee
    {
        double m_A { 0.0 };
        double m_B { 0.0 };
        double m_C { 0.0 };
        double m_D { 0.0 };
        double m_E { 0.0 };

    public:
        struct CreateFromDIPPR
        {
            double A;
            double B;
            double C;
            double D;
            double E;
        };

        /**
         * @brief
         */
        AlyLee();

        /**
         * @brief
         * @param coeffA
         * @param coeffB
         * @param coeffC
         * @param coeffD
         * @param coeffE
         */
        AlyLee(double coeffA, double coeffB, double coeffC, double coeffD, double coeffE);

        /**
         * @brief
         * @param coefficients
         */
        explicit AlyLee(const CreateFromDIPPR& coefficients);

        /**
         * @brief
         * @param other
         */
        AlyLee(const AlyLee& other);

        /**
         * @brief
         * @param other
         */
        AlyLee(AlyLee&& other) noexcept;

        /**
         * @brief
         */
        ~AlyLee();

        /**
         * @brief
         * @param other
         * @return
         */
        AlyLee& operator=(const AlyLee& other);

        /**
         * @brief
         * @param other
         * @return
         */
        AlyLee& operator=(AlyLee&& other) noexcept;

        /**
         * @brief
         * @param temperature
         * @return
         */
        double operator()(double temperature);

        /**
         * @brief
         * @param temperature
         * @return
         */
        double evaluateCp(double temperature) const;

        /**
         * @brief
         * @param temperature
         * @return
         */
        double derivativeOfCp(double temperature) const;

        /**
         * @brief
         * @param temperature
         * @return
         */
        double integralOfCp(double temperature) const;

        /**
         * @brief
         * @param temperature
         * @return
         */
        double integralOfCpOverT(double temperature) const;
    };

}    // namespace PCProps::HeatCapacity

#endif    // PCPROPS_ALYLEE_HPP
