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

#ifndef PCPROPS_EOSUTILITIES_HPP
#define PCPROPS_EOSUTILITIES_HPP

#include <string>
#include <tuple>
#include <vector>

namespace PCProps::EquationOfState
{
    enum PhaseDataElement {
        MolarFraction,
        Temperature,
        Pressure,
        Volume,
        FugacityCoefficient,
        Fugacity,
        Compressibility,
        Enthalpy,
        Entropy,
        InternalEnergy,
        GibbsEnergy,
        HelmholzEnergy
    };

    using PhaseData = std::tuple<
        double,  /* mole fraction [-] */
        double,  /* temperature [K] */
        double,  /* pressure [Pa] */
        double,  /* volume [m3/mol] */
        double,  /* fugacity coefficient [-] */
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