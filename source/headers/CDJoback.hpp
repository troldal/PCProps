//
// Created by Kenneth Balslev on 18/12/2020.
//

#ifndef PCPROPS_PDJOBACK_HPP
#define PCPROPS_PDJOBACK_HPP

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace PCProps::ConstantData
{
    struct CDJobackGroup
    {
        std::string           description;
        std::optional<double> tc;
        std::optional<double> pc;
        std::optional<double> vc;
        std::optional<double> tb;
        std::optional<double> tm;
        std::optional<double> hform;
        std::optional<double> gform;
        std::optional<double> igCp_a;
        std::optional<double> igCp_b;
        std::optional<double> igCp_c;
        std::optional<double> igCp_d;
        std::optional<double> hfus;
        std::optional<double> hvap;
        std::optional<double> liqVis_a;
        std::optional<double> liqVis_b;
    };

    /**
     * @brief
     */
    class CDJoback
    {
        std::vector<std::pair<int, CDJobackGroup>> m_groups {};

        std::optional<double> m_sumTc;
        std::optional<double> m_sumPc;
        std::optional<double> m_sumVc;
        std::optional<double> m_sumTb;
        std::optional<double> m_sumTm;
        std::optional<double> m_sumHform;
        std::optional<double> m_sumGform;
        std::optional<double> m_sumIgCp_a;
        std::optional<double> m_sumIgCp_b;
        std::optional<double> m_sumIgCp_c;
        std::optional<double> m_sumIgCp_d;
        std::optional<double> m_sumHfus;
        std::optional<double> m_sumHvap;
        std::optional<double> m_sumLiqVis_a;
        std::optional<double> m_sumLiqVis_b;

        double m_boilingTemperature;
        double m_molecularWeight = 0.0;
        int    m_atomCount       = 0;

    public:
        // ===== Constructors and Destructors =====

        /**
         * @brief
         */
        CDJoback();

        /**
         * @brief
         * @param groups
         * @param molecularWeight
         * @param atomCount
         */
        CDJoback(const std::vector<std::pair<int, int>>& groups, double molecularWeight, int atomCount);

        /**
         * @brief
         * @param other
         */
        CDJoback(const CDJoback& other);

        /**
         * @brief
         * @param other
         */
        CDJoback(CDJoback&& other) noexcept;

        /**
         * @brief
         */
        ~CDJoback();

        // ===== Manipulators =====

        /**
         * @brief
         * @param other
         * @return
         */
        CDJoback& operator=(const CDJoback& other);

        /**
         * @brief
         * @param other
         * @return
         */
        CDJoback& operator=(CDJoback&& other) noexcept;

        // ===== Accessors =====

        /**
         * @brief
         * @return
         */
        double boilingTemperature() const;

        /**
         * @brief
         * @return
         */
        double meltingTemperature() const;

        /**
         * @brief
         * @param boilingTemperature
         * @return
         */
        double criticalTemperature(double boilingTemperature = 0.0) const;

        /**
         * @brief
         * @return
         */
        double criticalPressure() const;

        /**
         * @brief
         * @return
         */
        double criticalVolume() const;

        /**
         * @brief
         * @return
         */
        double enthalpyOfFormation() const;

        /**
         * @brief
         * @return
         */
        double gibbsEnergyOfFormation() const;

        /**
         * @brief
         * @return
         */
        double enthalpyOfFusion() const;

        /**
         * @brief
         * @return
         */
        double enthalpyOfVaporization() const;

        /**
         * @brief
         * @param temperature
         * @return
         */
        double idealGasCp(double temperature) const;

        /**
         * @brief
         * @param temperature
         * @return
         */
        double liquidViscosity(double temperature) const;
    };

}    // namespace PCProps::ConstantData

#endif    // PCPROPS_PDJOBACK_HPP
