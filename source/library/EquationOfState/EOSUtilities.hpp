//
// Created by Kenneth Balslev on 29/12/2020.
//

#ifndef PCPROPS_EOSUTILITIES_HPP
#define PCPROPS_EOSUTILITIES_HPP

#include <string>

namespace PCProps::EquationOfState
{
    enum PhaseDataElement {
        Moles,
        MolecularWeight,
        Temperature,
        Pressure,
        Volume,
        Fugacity,
        Compressibility,
        Enthalpy,
        Entropy,
        InternalEnergy,
        GibbsEnergy,
        HelmholzEnergy
    };

    using PhaseData = std::tuple<
        double,  /* moles [-] */
        double,  /* molecular weight [g/mol] */
        double,  /* temperature [K] */
        double,  /* pressure [Pa] */
        double,  /* volume [m3/mol] */
        double,  /* fugacity [Pa] */
        double,  /* compressibility [-] */
        double,  /* enthalpy */
        double,  /* entropy */
        double,  /* internal energy */
        double,  /* Gibbs energy */
        double>; /* Helmholz energy */

    using Phases = std::vector<PhaseData>;

}    // namespace PCProps::EquationOfState

#endif    // PCPROPS_EOSUTILITIES_HPP
