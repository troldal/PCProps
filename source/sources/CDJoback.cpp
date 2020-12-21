//
// Created by Kenneth Balslev on 18/12/2020.
//

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <optional>
#include <tuple>

#include "CDJoback.hpp"

namespace
{
    using CDJobackGroupDef = std::tuple<std::string,
                                        std::optional<double>,
                                        std::optional<double>,
                                        std::optional<double>,
                                        std::optional<double>,
                                        std::optional<double>,
                                        std::optional<double>,
                                        std::optional<double>,
                                        std::optional<double>,
                                        std::optional<double>,
                                        std::optional<double>,
                                        std::optional<double>,
                                        std::optional<double>,
                                        std::optional<double>,
                                        std::optional<double>,
                                        std::optional<double>>;

    enum CDJobackProperty { Tc = 1, Pc, Vc, Tb, Tm, Hform, Gform, igCp_a, igCp_b, igCp_c, igCp_d, Hfus, Hvap, liqVis_a, liqVis_b };

    const std::array<CDJobackGroupDef, 41> JobackGroups {
        CDJobackGroupDef { "-CH3 (non-ring)", 0.0141, -0.0012, 65, 23.58, -5.10, -76.45, -43.96, 1.95E+1, -8.08E-3, 1.53E-4, -9.67E-8, 0.908, 2.373, 548.29, -1.719 },
        CDJobackGroupDef { "-CH2- (non-ring)", 0.0189, 0.0000, 56, 22.88, 11.27, -20.64, 8.42, -9.09E-1, 9.50E-2, -5.44E-5, 1.19E-8, 2.590, 2.226, 94.16, -0.199 },
        CDJobackGroupDef { ">CH- (non-ring)", 0.0164, 0.0020, 41, 21.74, 12.64, 29.89, 58.36, -2.30E+1, 2.04E-1, -2.65E-4, 1.20E-7, 0.749, 1.691, -322.15, 1.187 },
        CDJobackGroupDef { ">C< (non-ring)", 0.0067, 0.0043, 27, 18.25, 46.43, 82.23, 116.02, -6.62E+1, 4.27E-1, -6.41E-4, 3.01E-7, -1.460, 0.636, -573.56, 2.307 },
        CDJobackGroupDef { "=CH2 (non-ring)", 0.0113, -0.0028, 56, 18.18, -4.32, -9.630, 3.77, 2.36E+1, -3.81E-2, 1.72E-4, -1.03E-7, -0.473, 1.724, 495.01, -1.539 },
        CDJobackGroupDef { "=CH- (non-ring)", 0.0129, -0.0006, 46, 24.96, 8.73, 37.97, 48.53, -8.00, 1.05E-1, -9.63E-5, 3.56E-8, 2.691, 2.205, 82.28, -0.242 },
        CDJobackGroupDef { "=C< (non-ring)", 0.0117, 0.0011, 38, 24.14, 11.14, 83.99, 92.36, -2.81E+1, 2.08E-1, -3.06E-4, 1.46E-7, 3.063, 2.138, std::nullopt, std::nullopt },
        CDJobackGroupDef { "=C= (non-ring)", 0.0026, 0.0028, 36, 26.15, 17.78, 142.14, 136.70, 2.74E+1, -5.57E-2, 1.01E-4, -5.02E-8, 4.720, 2.661, std::nullopt, std::nullopt },
        CDJobackGroupDef { "≡CH (non-ring)", 0.0027, -0.0008, 46, 9.20, -11.18, 79.30, 77.71, 2.45E+1, -2.71E-2, 1.11E-4, -6.78E-8, 2.322, 1.155, std::nullopt, std::nullopt },
        CDJobackGroupDef { "≡C- (non-ring)", 0.0020, 0.0016, 37, 27.38, 64.32, 115.51, 109.82, 7.87, 2.01E-2, -8.33E-6, 1.39E-9, 4.151, 3.302, std::nullopt, std::nullopt },
        CDJobackGroupDef { "-CH2- (ring)", 0.0100, 0.0025, 48, 27.15, 7.75, -26.80, -3.68, -6.03, 8.54E-2, -8.00E-6, -1.80E-8, 0.490, 2.398, 307.53, -0.798 },
        CDJobackGroupDef { ">CH- (ring)", 0.0122, 0.0004, 38, 21.78, 19.88, 8.67, 40.99, -2.05E+1, 1.62E-1, -1.60E-4, 6.24E-8, 3.243, 1.942, -394.29, 1.251 },
        CDJobackGroupDef { ">C< (ring)", 0.0042, 0.0061, 27, 21.32, 60.15, 79.72, 87.88, -9.09E+1, 5.57E-1, -9.00E-4, 4.69E-7, -1.373, 0.644, std::nullopt, std::nullopt },
        CDJobackGroupDef { "=CH- (ring)", 0.0082, 0.0011, 41, 26.73, 8.13, 2.09, 11.30, -2.14, 5.74E-2, -1.64E-6, -1.59E-8, 1.101, 2.544, 259.65, -0.702 },
        CDJobackGroupDef { "=C< (ring)", 0.0143, 0.0008, 32, 31.01, 37.02, 46.43, 54.05, -8.25, 1.01E-1, -1.42E-4, 6.78E-8, 2.394, 3.059, -245.74, 0.912 },
        CDJobackGroupDef { "-F (halogen)", 0.0111, -0.0057, 27, -0.03, -15.78, -251.92, -247.19, 2.65E+1, -9.13E-2, 1.91E-4, -1.03E-7, 1.398, -0.670, std::nullopt, std::nullopt },
        CDJobackGroupDef { "-Cl (halogen)", 0.0105, -0.0049, 58, 38.13, 13.55, -71.55, -64.31, 3.33E+1, -9.63E-2, 1.87E-4, -9.96E-8, 2.515, 4.532, 625.45, -1.814 },
        CDJobackGroupDef { "-Br (halogen)", 0.0133, 0.0057, 71, 66.86, 43.43, -29.48, -38.06, 2.86E+1, -6.49E-2, 1.36E-4, -7.45E-8, 3.603, 6.582, 738.91, -2.038 },
        CDJobackGroupDef { "-I (halogen)", 0.0068, -0.0034, 97, 93.84, 41.69, 21.06, 5.74, 3.21E+1, -6.41E-2, 1.26E-4, -6.87E-8, 2.724, 9.520, 809.55, -2.224 },
        CDJobackGroupDef { "-OH (alcohol)", 0.0741, 0.0112, 28, 92.88, 44.45, -208.04, -189.20, 2.57E+1, -6.91E-2, 1.77E-4, -9.88E-8, 2.406, 16.826, 2173.72, -5.057 },
        CDJobackGroupDef { "-OH (phenol)", 0.0240, 0.0184, -25, 76.34, 82.83, -221.65, -197.37, -2.81, 1.11E-1, -1.16E-4, 4.94E-8, 4.490, 12.499, 3018.17, -7.314 },
        CDJobackGroupDef { "-O- (non-ring)", 0.0168, 0.0015, 18, 22.42, 22.23, -132.22, -105.00, 2.55E+1, -6.32E-2, 1.11E-4, -5.48E-8, 1.188, 2.410, 122.09, -0.386 },
        CDJobackGroupDef { "-O- (ring)", 0.0098, 0.0048, 13, 31.22, 23.05, -138.16, -98.22, 1.22E+1, -1.26E-2, 6.03E-5, -3.86E-8, 5.879, 4.682, 440.24, -0.953 },
        CDJobackGroupDef { ">C=O (non-ring)", 0.0380, 0.0031, 62, 76.75, 61.20, -133.22, -120.50, 6.45, 6.70E-2, -3.57E-5, 2.86E-9, 4.189, 8.972, 340.35, -0.350 },
        CDJobackGroupDef { ">C=O (ring)", 0.0284, 0.0028, 55, 94.97, 75.97, -164.50, -126.27, 3.04E+1, -8.29E-2, 2.36E-4, -1.31E-7, 0., 6.645, std::nullopt, std::nullopt },
        CDJobackGroupDef { "O=CH- (aldehyde)", 0.0379, 0.0030, 82, 72.24, 36.90, -162.03, -143.48, 3.09E+1, -3.36E-2, 1.60E-4, -9.88E-8, 3.197, 9.093, 740.92, -1.713 },
        CDJobackGroupDef { "-COOH (acid)", 0.0791, 0.0077, 89, 169.09, 155.50, -426.72, -387.87, 2.41E+1, 4.27E-2, 8.04E-5, -6.87E-8, 11.051, 19.537, 1317.23, -2.578 },
        CDJobackGroupDef { "-COO- (ester)", 0.0481, 0.0005, 82, 81.10, 53.60, -337.92, -301.95, 2.45E+1, 4.02E-2, 4.02E-5, -4.52E-8, 6.959, 9.633, 483.88, -0.966 },
        CDJobackGroupDef { "=O (general)", 0.0143, 0.0101, 36, -10.50, 2.08, -247.61, -250.83, 6.82, 1.96E-2, 1.27E-5, -1.78E-8, 3.624, 5.909, 675.24, -1.340 },
        CDJobackGroupDef { "-NH2 (general)", 0.0243, 0.0109, 38, 73.23, 66.89, -22.02, 14.07, 2.69E+1, -4.12E-2, 1.64E-4, -9.76E-8, 3.515, 10.788, std::nullopt, std::nullopt },
        CDJobackGroupDef { ">NH (non-ring)", 0.0295, 0.0077, 35, 50.17, 52.66, 53.47, 89.39, -1.21, 7.62E-2, -4.86E-5, 1.05E-8, 5.099, 6.436, std::nullopt, std::nullopt },
        CDJobackGroupDef { ">NH (ring)", 0.0130, 0.0114, 29, 52.82, 101.51, 31.65, 75.61, 1.18E+1, -2.30E-2, 1.07E-4, -6.28E-8, 7.490, 6.930, std::nullopt, std::nullopt },
        CDJobackGroupDef { ">N- (non-ring)", 0.0169, 0.0074, 9, 11.74, 48.84, 123.34, 163.16, -3.11E+1, 2.27E-1, -3.20E-4, 1.46E-7, 4.703, 1.896, std::nullopt, std::nullopt },
        CDJobackGroupDef { "-N= (non-ring)", 0.0255, -0.0099, std::nullopt, 74.60, std::nullopt, 23.61, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, 3.335, std::nullopt, std::nullopt },
        CDJobackGroupDef { "-N= (ring)", 0.0085, 0.0076, 34, 57.55, 68.40, 55.52, 79.93, 8.83, -3.84E-3, 4.35E-5, -2.60E-8, 3.649, 6.528, std::nullopt, std::nullopt },
        CDJobackGroupDef { "=NH (general)", std::nullopt, std::nullopt, std::nullopt, 83.08, 68.91, 93.70, 119.66, 5.69, -4.12E-3, 1.28E-4, -8.88E-8, std::nullopt, 12.169, std::nullopt, std::nullopt },
        CDJobackGroupDef { "-CN (general)", 0.0496, -0.0101, 91, 125.66, 59.89, 88.43, 89.22, 3.65E+1, -7.33E-2, 1.84E-4, -1.03E-7, 2.414, 12.851, std::nullopt, std::nullopt },
        CDJobackGroupDef { "-NO2 (general)", 0.0437, 0.0064, 91, 152.54, 127.24, -66.57, -16.83, 2.59E+1, -3.74E-3, 1.29E-4, -8.88E-8, 9.679, 16.738, std::nullopt, std::nullopt },
        CDJobackGroupDef { "-SH (general)", 0.0031, 0.0084, 63, 63.56, 20.09, -17.33, -22.99, 3.53E+1, -7.58E-2, 1.85E-4, -1.03E-7, 2.360, 6.884, std::nullopt, std::nullopt },
        CDJobackGroupDef { "-S- (non-ring)", 0.0119, 0.0049, 54, 68.78, 34.40, 41.87, 33.12, 1.96E+1, -5.61E-3, 4.02E-5, -2.76E-8, 4.130, 6.817, std::nullopt, std::nullopt },
        CDJobackGroupDef { "-S- (ring)", 0.0019, 0.0051, 38, 52.10, 79.93, 39.10, 27.76, 1.67E+1, 4.81E-3, 2.77E-5, -2.11E-8, 1.557, 5.984, std::nullopt, std::nullopt }
    };

    template<int index>
    std::optional<double> accumulate(const std::vector<std::pair<int, int>>& groups)
    {
        std::vector<std::optional<double>> items {};
        items.reserve(groups.size());

        for (const auto& item : groups)
            items.emplace_back(std::get<index>(JobackGroups.at(static_cast<uint64_t>(item.first - 1))).has_value()
                                   ? std::optional(static_cast<double>(item.second) * std::get<index>(JobackGroups.at(static_cast<uint64_t>(item.first - 1))).value())
                                   : std::nullopt);

        return (std::find_if(items.begin(), items.end(), [](const std::optional<double>& item) { return !item.has_value(); }) == items.end()
                    ? std::optional(std::accumulate(items.begin(), items.end(), 0.0, [](double result, const std::optional<double>& item) { return result + item.value(); }))
                    : std::nullopt);
    }

}    // namespace

namespace PCProps::ConstantData
{
    // ===== Constructor, default
    CDJoback::CDJoback() = default;

    // ===== Constructor, taking groups as arguments
    CDJoback::CDJoback(const std::vector<CDJobackGroup>& groups, double molecularWeight, int atomCount, double boilingTemperature)
        : m_sumTc { accumulate<Tc>(groups) },
          m_sumPc { accumulate<Pc>(groups) },
          m_sumVc { accumulate<Vc>(groups) },
          m_sumTb { accumulate<Tb>(groups) },
          m_sumTm { accumulate<Tm>(groups) },
          m_sumHform { accumulate<Hform>(groups) },
          m_sumGform { accumulate<Gform>(groups) },
          m_sumIgCp_a { accumulate<igCp_a>(groups) },
          m_sumIgCp_b { accumulate<igCp_b>(groups) },
          m_sumIgCp_c { accumulate<igCp_c>(groups) },
          m_sumIgCp_d { accumulate<igCp_d>(groups) },
          m_sumHfus { accumulate<Hfus>(groups) },
          m_sumHvap { accumulate<Hvap>(groups) },
          m_sumLiqVis_a { accumulate<liqVis_a>(groups) },
          m_sumLiqVis_b { accumulate<liqVis_b>(groups) },
          m_molecularWeight { molecularWeight },
          m_atomCount { atomCount },
          m_boilingTemperature { boilingTemperature }
    {}

    // ===== Copy constructor
    CDJoback::CDJoback(const CDJoback& other) = default;

    // ===== Move constructor
    CDJoback::CDJoback(CDJoback&& other) noexcept = default;

    // ===== Destructor
    CDJoback::~CDJoback() = default;

    // ===== Copy assignment operator
    CDJoback& CDJoback::operator=(const CDJoback& other) = default;

    // ===== Move assignment operator
    CDJoback& CDJoback::operator=(CDJoback&& other) noexcept = default;

    // ===== Compute normal boiling point temperature [K] for substance
    double CDJoback::boilingTemperature() const
    {
        if (static_cast<bool>(m_boilingTemperature)) return m_boilingTemperature;
        return 198.2 + m_sumTb.value();
    }

    // ===== Check that normal boiling point temperature can be computed
    bool CDJoback::boilingTemperatureIsValid() const
    {
        return m_sumTb.has_value();
    }

    // ===== Compute melting temperature [K] for substance
    double CDJoback::meltingTemperature() const
    {
        return 122.5 + m_sumTm.value();
    }

    // ===== Check that melting temperature can be computed
    bool CDJoback::meltingTemperatureIsValid() const
    {
        return m_sumTm.has_value();
    }

    // ===== Compute critical temperature [K] for substance
    double CDJoback::criticalTemperature() const
    {
        using std::pow;
        return boilingTemperature() / (0.584 + 0.965 * m_sumTc.value() - pow(m_sumTc.value(), 2));
    }

    // ===== Check that critical temperature can be computed
    bool CDJoback::criticalTemperatureIsValid() const
    {
        return m_sumTc.has_value();
    }

    // ===== Compute critical pressure [Pa] for substance
    double CDJoback::criticalPressure() const
    {
        using std::pow;
        return pow(0.113 + 0.0032 * m_atomCount - m_sumPc.value(), -2) * 100000;
    }

    // ===== Check that critical pressure can be computed
    bool CDJoback::criticalPressureIsValid() const
    {
        return m_sumPc.has_value();
    }

    // ===== Compute critical volume [m3/mol] for substance
    double CDJoback::criticalVolume() const
    {
        return (17.5 + m_sumVc.value()) / 1000000;
    }

    // ===== Check that critical volume can be computed
    bool CDJoback::criticalVolumeIsValid() const
    {
        return m_sumVc.has_value();
    }

    // ===== Compute heat of formation [J/mol] for substance
    double CDJoback::enthalpyOfFormation() const
    {
        return (68.29 + m_sumHform.value()) * 1000;
    }

    // ===== Check that heat of formation can be computed
    bool CDJoback::enthalpyOfFormationIsValid() const
    {
        return m_sumHform.has_value();
    }

    // ===== Compute Gibbs energy of formation [J/mol] for substance
    double CDJoback::gibbsEnergyOfFormation() const
    {
        return (53.88 + m_sumGform.value()) * 1000;
    }

    // ===== Check that Gibbs energy of formation can be computed
    bool CDJoback::gibbsEnergyOfFormationIsValid() const
    {
        return m_sumGform.has_value();
    }

    // ===== Compute heat of fusion [J/mol] for substance
    double CDJoback::enthalpyOfFusion() const
    {
        return (-0.88 + m_sumHfus.value()) * 1000;
    }

    // ===== Check that enthalpy of fusion can be computed
    bool CDJoback::enthalpyOfFusionIsValid() const
    {
        return m_sumHfus.has_value();
    }

    // ===== Compute heat of vaporization [J/mol] for substance
    double CDJoback::enthalpyOfVaporization() const
    {
        return (15.30 + m_sumHvap.value()) * 1000;
    }

    // ===== Check that enthalpy of vaporization can be computed
    bool CDJoback::enthalpyOfVaporizationIsValid() const
    {
        return m_sumHvap.has_value();
    }

    // ===== Compute ideal gas Cp [J/(mol K)] for substance
    double CDJoback::idealGasCp(double temperature) const
    {
        return m_sumIgCp_a.value() - 37.93 + (m_sumIgCp_b.value() + 0.21) * temperature + (m_sumIgCp_c.value() - 3.91E-4) * pow(temperature, 2) + (m_sumIgCp_d.value() + 2.06E-7) * pow(temperature, 3);
    }

    // ===== Check that ideal gas Cp can be computed
    bool CDJoback::idealGasCpIsValid() const
    {
        return m_sumIgCp_a.has_value() && m_sumIgCp_b.has_value() && m_sumIgCp_c.has_value() && m_sumIgCp_d.has_value();
    }

    // ====== Compute liquid viscosity [Pa s] for substance.
    double CDJoback::liquidViscosity(double temperature) const
    {
        using std::exp;
        return m_molecularWeight * exp((m_sumLiqVis_a.value() - 597.82) / temperature + m_sumLiqVis_b.value() - 11.202);
    }

    // ===== Check that liquid viscosity can be computed
    bool CDJoback::liquidViscosityIsValid() const
    {
        return m_sumLiqVis_a.has_value() && m_sumLiqVis_b.has_value();
    }

}    // namespace PCProps::ConstantData