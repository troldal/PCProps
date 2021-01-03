/*

8888888b.   .d8888b.  8888888b.
888   Y88b d88P  Y88b 888   Y88b
888    888 888    888 888    888
888   d88P 888        888   d88P 888d888 .d88b.  88888b.  .d8888b
8888888P"  888        8888888P"  888P"  d88""88b 888 "88b 88K
888        888    888 888        888    888  888 888  888 "Y8888b.
888        Y88b  d88P 888        888    Y88..88P 888 d88P      X88
888         "Y8888P"  888        888     "Y88P"  88888P"   88888P'
                                                 888
                                                 888
                                                 888

Copyright (c) 2020 Kenneth Troldal Balslev

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#ifndef PCPROPS_PDJOBACK_HPP
#define PCPROPS_PDJOBACK_HPP

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace PCProps::ConstantData
{
    /**
     * @brief The CDJobackGroup is an alias for a std::pair<int, int>, used to define a single Joback group
     * for a component. The first element is the Joback group index, and the second element is the count. For
     * example, Acetone consists of two methyl groups and one keton group. Methyl has an index of 1 and keton
     * has an index 24. Hence the CDJobackGroup definitions are <1, 2> and <24, 1>
     */
    using CDJobackGroup = std::pair<int, int>;

    /**
     * @brief The CDJoback class implements the Joback group contribution method for estimating pure component properties.
     * @details The Joback method can be used to estimate critical temperature, critical pressure, critical volume, normal
     * boiling point, melting point, enthalpy of formation, Gibbs energy of formation, computeEnthalpy of fusion, enthalpy of vaporization,
     * ideal gas Cp and liquid viscosity. Documentation of how to use can be found on https://en.wikipedia.org/wiki/Joback_method
     * and other places. The original M.Sc. thesis by Joback can be found on MIT's website: https://dspace.mit.edu/handle/1721.1/15374.
     *
     * The group definition are tabulated below:
     * Index  | Group
     * :----- | :-----
     * NON-RING GROUPS ||
     * 1  | −CH3
     * 2  | −CH2−
     * 3  | >CH−
     * 4  | >C<
     * 5  | =CH2
     * 6  | =CH−
     * 7  | =C<
     * 8  | =C=
     * 9  | ≡CH
     * 10 | ≡C−
     * RING GROUPS ||
     * 11 | −CH2−
     * 12 | >CH−
     * 13 | >C<
     * 14 | =CH−
     * 15 | =C<
     * HALOGEN GROUPS ||
     * 16 | −F
     * 17 | −Cl
     * 18 | −Br
     * 19 | −I
     * OXYGEN GROUPS ||
     * 20 | −OH (alcohol)
     * 21 | −OH (phenol)
     * 22 | −O− (non-ring)
     * 23 | −O− (ring)
     * 24 | >C=O (non-ring)
     * 25 | >C=O (ring)
     * 26 | −CH=O (aldehyde)
     * 27 | −COOH (acid)
     * 28 | −COO− (ester)
     * 29 | =O (other than above)
     * NITROGEN GROUPS ||
     * 30 | −NH2
     * 31 | >NH (non-ring)
     * 32 | >NH (ring)
     * 33 | >N− (non-ring)
     * 34 | −N= (non-ring)
     * 35 | −N= (ring)
     * 36 | =NH
     * 37 | −CN
     * 38 | −NO2
     * SULFUR GROUPS ||
     * 39 | −SH
     * 40 | −S− (non-ring)
     * 41 | −S− (ring)
     *
     */
    class CDJoback
    {
        std::optional<double> m_sumTc;       /**< The sum of Joback terms for Tc estimation. */
        std::optional<double> m_sumPc;       /**< The sum of Joback terms for Pc estimation. */
        std::optional<double> m_sumVc;       /**< The sum of Joback terms for Vc estimation. */
        std::optional<double> m_sumTb;       /**< The sum of Joback terms for Tb estimation. */
        std::optional<double> m_sumTm;       /**< The sum of Joback terms for Tm, estimation. */
        std::optional<double> m_sumHform;    /**< The sum of Joback terms for computeEnthalpy for formation estimation. */
        std::optional<double> m_sumGform;    /**< The sum of Joback terms for Gibbs energy of formation estimation. */
        std::optional<double> m_sumIgCp_a;   /**< The sum of Joback terms for ideal gas Cp (a) estimation. */
        std::optional<double> m_sumIgCp_b;   /**< The sum of Joback terms for ideal gas Cp (b) estimation. */
        std::optional<double> m_sumIgCp_c;   /**< The sum of Joback terms for ideal gas Cp (c) estimation. */
        std::optional<double> m_sumIgCp_d;   /**< The sum of Joback terms for ideal gas Cp (d) estimation. */
        std::optional<double> m_sumHfus;     /**< The sum of Joback terms for computeEnthalpy of fusion estimation. */
        std::optional<double> m_sumHvap;     /**< The sum of Joback terms for computeEnthalpy of vaporization estimation. */
        std::optional<double> m_sumLiqVis_a; /**< The sum of Joback terms for liquid viscosity (a) estimation. */
        std::optional<double> m_sumLiqVis_b; /**< The sum of Joback terms for liquid viscosity (b) estimation. */

        double m_molecularWeight    = 0.0; /**< The molecular weight of the component. */
        int    m_atomCount          = 0;   /**< The number of atoms in the component. */
        double m_boilingTemperature = 0.0; /**< The normal boiling point of the component. */

    public:
        // ===== Constructors and Destructors ===== //

        /**
         * @brief Constructor, default
         */
        CDJoback();

        /**
         * @brief Constructor, taking Joback groups, molecular weight and atom count as arguments.
         * @param groups A std::vector with CDJobackGroups (aka std::pair<int, int>) with group index and count.
         * @param molecularWeight The molecular weight of the compound.
         * @param atomCount The number of atoms in the compound.
         * @param boilingTemperature The normal boiling point of the component (needed for Tc estimation).
         * If not available, the normal boiling point can be estimated using the boilingTemperature() member function.
         */
        CDJoback(const std::vector<CDJobackGroup>& groups, double molecularWeight, int atomCount, double boilingTemperature = 0.0);

        /**
         * @brief Constructor, taking Joback groups, molecular weight and atom count as arguments.
         * @details This constructor template can take any container supporting range-based for loops as an argument.
         * If using a map (e.g. std::map or std::unordered_map), the Joback group index is the key, and the count is the value.
         * @tparam Container The type of container holding the CDJobackGroups.
         * @param groups A container (e.g. std::list or std::map) with CDJobackGroups (aka std::pair<int, int>) with group index and count.
         * @param molecularWeight The molecular weight of the compound.
         * @param atomCount The number of atoms in the compound.
         */
        template<typename Container>
        CDJoback(const Container& groups, double molecularWeight, int atomCount)
        {
            std::vector<CDJobackGroup> groupvec;
            for (const auto& item : groups) groupvec.template emplace_back(item);

            *this = CDJoback(groupvec, molecularWeight, atomCount, 0.0);
        }

        /**
         * @brief Copy constructor.
         */
        CDJoback(const CDJoback& other);

        /**
         * @brief Move constructor.
         */
        CDJoback(CDJoback&& other) noexcept;

        /**
         * @brief Destructor.
         */
        ~CDJoback();

        // ===== Manipulators ===== //

        /**
         * @brief Copy assignment operator.
         */
        CDJoback& operator=(const CDJoback& other);

        /**
         * @brief Move assignment operator.
         */
        CDJoback& operator=(CDJoback&& other) noexcept;

        // ===== Accessors ===== //

        /**
         * @brief Estimate the normal boiling point temperature [K].
         * @return The normal boiling point temperature [K]
         */
        double boilingTemperature() const;

        /**
         * @brief
         * @return
         */
        bool boilingTemperatureIsValid() const;

        /**
         * @brief Estimate the melting temperature [K].
         * @return The melting temperature [K].
         */
        double meltingTemperature() const;

        /**
         * @brief
         * @return
         */
        bool meltingTemperatureIsValid() const;

        /**
         * @brief Estimate the critical temperature [K].
         * @return The critical temperature [K].
         */
        double criticalTemperature() const;

        /**
         * @brief
         * @return
         */
        bool criticalTemperatureIsValid() const;

        /**
         * @brief Estimate the critical pressure [Pa].
         * @return The critical pressure [Pa].
         */
        double criticalPressure() const;

        /**
         * @brief
         * @return
         */
        bool criticalPressureIsValid() const;

        /**
         * @brief Estimate the critical volume [m3/mol].
         * @return The critical volume [m3/mol].
         */
        double criticalVolume() const;

        /**
         * @brief
         * @return
         */
        bool criticalVolumeIsValid() const;

        /**
         * @brief Estimate the computeEnthalpy of formation [J/mol].
         * @return The computeEnthalpy of formation [J/mol].
         */
        double enthalpyOfFormation() const;

        /**
         * @brief
         * @return
         */
        bool enthalpyOfFormationIsValid() const;

        /**
         * @brief Estimate the Gibbs energy of formation [J/mol].
         * @return The Gibbs energy of formation [J/mol].
         */
        double gibbsEnergyOfFormation() const;

        /**
         * @brief
         * @return
         */
        bool gibbsEnergyOfFormationIsValid() const;

        /**
         * @brief Estimate the computeEnthalpy of fusion [J/mol].
         * @return The computeEnthalpy of fusion [J/mol].
         */
        double enthalpyOfFusion() const;

        /**
         * @brief
         * @return
         */
        bool enthalpyOfFusionIsValid() const;

        /**
         * @brief Estimate the computeEnthalpy of vaporization [J/mol].
         * @return The computeEnthalpy of vaporization [J/mol].
         */
        double enthalpyOfVaporization() const;

        /**
         * @brief
         * @return
         */
        bool enthalpyOfVaporizationIsValid() const;

        /**
         * @brief Estimate the ideal gas Cp at a given temperature.
         * @param temperature The temperature [K] at which to estimate the ideal gas Cp.
         * @return The ideal gas Cp [J/(mol K)].
         */
        double idealGasCp(double temperature) const;

        /**
         * @brief
         * @return
         */
        bool idealGasCpIsValid() const;

        /**
         * @brief Estimate the liquid viscosity at a given temperature.
         * @param temperature The temperature [K] at which to estimate the liquid viscosity.
         * @return The liquid viscosity [Pa s]
         */
        double liquidViscosity(double temperature) const;

        /**
         * @brief
         * @return
         */
        bool liquidViscosityIsValid() const;
    };

}    // namespace PCProps::ConstantData

#endif    // PCPROPS_PDJOBACK_HPP
