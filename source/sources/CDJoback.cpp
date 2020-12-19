//
// Created by Kenneth Balslev on 18/12/2020.
//

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <optional>
#include <string>

#include "CDJoback.hpp"

namespace
{
    using PCProps::ConstantData::CDJobackGroup;
    const std::array<CDJobackGroup, 41> JobackGroups {
        CDJobackGroup { "-CH3 (non-ring)", 0.0141, -0.0012, 65, 23.58, -5.10, -76.45, -43.96, 1.95E+1, -8.08E-3, 1.53E-4, -9.67E-8, 0.908, 2.373, 548.29, -1.719 },
        CDJobackGroup { "-CH2- (non-ring)", 0.0189, 0.0000, 56, 22.88, 11.27, -20.64, 8.42, -9.09E-1, 9.50E-2, -5.44E-5, 1.19E-8, 2.590, 2.226, 94.16, -0.199 },
        CDJobackGroup { ">CH- (non-ring)", 0.0164, 0.0020, 41, 21.74, 12.64, 29.89, 58.36, -2.30E+1, 2.04E-1, -2.65E-4, 1.20E-7, 0.749, 1.691, -322.15, 1.187 },
        CDJobackGroup { ">C< (non-ring)", 0.0067, 0.0043, 27, 18.25, 46.43, 82.23, 116.02, -6.62E+1, 4.27E-1, -6.41E-4, 3.01E-7, -1.460, 0.636, -573.56, 2.307 },
        CDJobackGroup { "=CH2 (non-ring)", 0.0113, -0.0028, 56, 18.18, -4.32, -9.630, 3.77, 2.36E+1, -3.81E-2, 1.72E-4, -1.03E-7, -0.473, 1.724, 495.01, -1.539 },
        CDJobackGroup { "=CH- (non-ring)", 0.0129, -0.0006, 46, 24.96, 8.73, 37.97, 48.53, -8.00, 1.05E-1, -9.63E-5, 3.56E-8, 2.691, 2.205, 82.28, -0.242 },
        CDJobackGroup { "=C< (non-ring)", 0.0117, 0.0011, 38, 24.14, 11.14, 83.99, 92.36, -2.81E+1, 2.08E-1, -3.06E-4, 1.46E-7, 3.063, 2.138, std::nullopt, std::nullopt },
        CDJobackGroup { "=C= (non-ring)", 0.0026, 0.0028, 36, 26.15, 17.78, 142.14, 136.70, 2.74E+1, -5.57E-2, 1.01E-4, -5.02E-8, 4.720, 2.661, std::nullopt, std::nullopt },
        CDJobackGroup { "≡CH (non-ring)", 0.0027, -0.0008, 46, 9.20, -11.18, 79.30, 77.71, 2.45E+1, -2.71E-2, 1.11E-4, -6.78E-8, 2.322, 1.155, std::nullopt, std::nullopt },
        CDJobackGroup { "≡C- (non-ring)", 0.0020, 0.0016, 37, 27.38, 64.32, 115.51, 109.82, 7.87, 2.01E-2, -8.33E-6, 1.39E-9, 4.151, 3.302, std::nullopt, std::nullopt },
        CDJobackGroup { "-CH2- (ring)", 0.0100, 0.0025, 48, 27.15, 7.75, -26.80, -3.68, -6.03, 8.54E-2, -8.00E-6, -1.80E-8, 0.490, 2.398, 307.53, -0.798 },
        CDJobackGroup { ">CH- (ring)", 0.0122, 0.0004, 38, 21.78, 19.88, 8.67, 40.99, -2.05E+1, 1.62E-1, -1.60E-4, 6.24E-8, 3.243, 1.942, -394.29, 1.251 },
        CDJobackGroup { ">C< (ring)", 0.0042, 0.0061, 27, 21.32, 60.15, 79.72, 87.88, -9.09E+1, 5.57E-1, -9.00E-4, 4.69E-7, -1.373, 0.644, std::nullopt, std::nullopt },
        CDJobackGroup { "=CH- (ring)", 0.0082, 0.0011, 41, 26.73, 8.13, 2.09, 11.30, -2.14, 5.74E-2, -1.64E-6, -1.59E-8, 1.101, 2.544, 259.65, -0.702 },
        CDJobackGroup { "=C< (ring)", 0.0143, 0.0008, 32, 31.01, 37.02, 46.43, 54.05, -8.25, 1.01E-1, -1.42E-4, 6.78E-8, 2.394, 3.059, -245.74, 0.912 },
        CDJobackGroup { "-F (halogen)", 0.0111, -0.0057, 27, -0.03, -15.78, -251.92, -247.19, 2.65E+1, -9.13E-2, 1.91E-4, -1.03E-7, 1.398, -0.670, std::nullopt, std::nullopt },
        CDJobackGroup { "-Cl (halogen)", 0.0105, -0.0049, 58, 38.13, 13.55, -71.55, -64.31, 3.33E+1, -9.63E-2, 1.87E-4, -9.96E-8, 2.515, 4.532, 625.45, -1.814 },
        CDJobackGroup { "-Br (halogen)", 0.0133, 0.0057, 71, 66.86, 43.43, -29.48, -38.06, 2.86E+1, -6.49E-2, 1.36E-4, -7.45E-8, 3.603, 6.582, 738.91, -2.038 },
        CDJobackGroup { "-I (halogen)", 0.0068, -0.0034, 97, 93.84, 41.69, 21.06, 5.74, 3.21E+1, -6.41E-2, 1.26E-4, -6.87E-8, 2.724, 9.520, 809.55, -2.224 },
        CDJobackGroup { "-OH (alcohol)", 0.0741, 0.0112, 28, 92.88, 44.45, -208.04, -189.20, 2.57E+1, -6.91E-2, 1.77E-4, -9.88E-8, 2.406, 16.826, 2173.72, -5.057 },
        CDJobackGroup { "-OH (phenol)", 0.0240, 0.0184, -25, 76.34, 82.83, -221.65, -197.37, -2.81, 1.11E-1, -1.16E-4, 4.94E-8, 4.490, 12.499, 3018.17, -7.314 },
        CDJobackGroup { "-O- (non-ring)", 0.0168, 0.0015, 18, 22.42, 22.23, -132.22, -105.00, 2.55E+1, -6.32E-2, 1.11E-4, -5.48E-8, 1.188, 2.410, 122.09, -0.386 },
        CDJobackGroup { "-O- (ring)", 0.0098, 0.0048, 13, 31.22, 23.05, -138.16, -98.22, 1.22E+1, -1.26E-2, 6.03E-5, -3.86E-8, 5.879, 4.682, 440.24, -0.953 },
        CDJobackGroup { ">C=O (non-ring)", 0.0380, 0.0031, 62, 76.75, 61.20, -133.22, -120.50, 6.45, 6.70E-2, -3.57E-5, 2.86E-9, 4.189, 8.972, 340.35, -0.350 },
        CDJobackGroup { ">C=O (ring)", 0.0284, 0.0028, 55, 94.97, 75.97, -164.50, -126.27, 3.04E+1, -8.29E-2, 2.36E-4, -1.31E-7, 0., 6.645, std::nullopt, std::nullopt },
        CDJobackGroup { "O=CH- (aldehyde)", 0.0379, 0.0030, 82, 72.24, 36.90, -162.03, -143.48, 3.09E+1, -3.36E-2, 1.60E-4, -9.88E-8, 3.197, 9.093, 740.92, -1.713 },
        CDJobackGroup { "-COOH (acid)", 0.0791, 0.0077, 89, 169.09, 155.50, -426.72, -387.87, 2.41E+1, 4.27E-2, 8.04E-5, -6.87E-8, 11.051, 19.537, 1317.23, -2.578 },
        CDJobackGroup { "-COO- (ester)", 0.0481, 0.0005, 82, 81.10, 53.60, -337.92, -301.95, 2.45E+1, 4.02E-2, 4.02E-5, -4.52E-8, 6.959, 9.633, 483.88, -0.966 },
        CDJobackGroup { "=O (general)", 0.0143, 0.0101, 36, -10.50, 2.08, -247.61, -250.83, 6.82, 1.96E-2, 1.27E-5, -1.78E-8, 3.624, 5.909, 675.24, -1.340 },
        CDJobackGroup { "-NH2 (general)", 0.0243, 0.0109, 38, 73.23, 66.89, -22.02, 14.07, 2.69E+1, -4.12E-2, 1.64E-4, -9.76E-8, 3.515, 10.788, std::nullopt, std::nullopt },
        CDJobackGroup { ">NH (non-ring)", 0.0295, 0.0077, 35, 50.17, 52.66, 53.47, 89.39, -1.21, 7.62E-2, -4.86E-5, 1.05E-8, 5.099, 6.436, std::nullopt, std::nullopt },
        CDJobackGroup { ">NH (ring)", 0.0130, 0.0114, 29, 52.82, 101.51, 31.65, 75.61, 1.18E+1, -2.30E-2, 1.07E-4, -6.28E-8, 7.490, 6.930, std::nullopt, std::nullopt },
        CDJobackGroup { ">N- (non-ring)", 0.0169, 0.0074, 9, 11.74, 48.84, 123.34, 163.16, -3.11E+1, 2.27E-1, -3.20E-4, 1.46E-7, 4.703, 1.896, std::nullopt, std::nullopt },
        CDJobackGroup { "-N= (non-ring)",
                        0.0255,
                        -0.0099,
                        std::nullopt,
                        74.60,
                        std::nullopt,
                        23.61,
                        std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        3.335,
                        std::nullopt,
                        std::nullopt },
        CDJobackGroup { "-N= (ring)", 0.0085, 0.0076, 34, 57.55, 68.40, 55.52, 79.93, 8.83, -3.84E-3, 4.35E-5, -2.60E-8, 3.649, 6.528, std::nullopt, std::nullopt },
        CDJobackGroup { "=NH (general)",
                        std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        83.08,
                        68.91,
                        93.70,
                        119.66,
                        5.69,
                        -4.12E-3,
                        1.28E-4,
                        -8.88E-8,
                        std::nullopt,
                        12.169,
                        std::nullopt,
                        std::nullopt },
        CDJobackGroup { "-CN (general)", 0.0496, -0.0101, 91, 125.66, 59.89, 88.43, 89.22, 3.65E+1, -7.33E-2, 1.84E-4, -1.03E-7, 2.414, 12.851, std::nullopt, std::nullopt },
        CDJobackGroup { "-NO2 (general)", 0.0437, 0.0064, 91, 152.54, 127.24, -66.57, -16.83, 2.59E+1, -3.74E-3, 1.29E-4, -8.88E-8, 9.679, 16.738, std::nullopt, std::nullopt },
        CDJobackGroup { "-SH (general)", 0.0031, 0.0084, 63, 63.56, 20.09, -17.33, -22.99, 3.53E+1, -7.58E-2, 1.85E-4, -1.03E-7, 2.360, 6.884, std::nullopt, std::nullopt },
        CDJobackGroup { "-S- (non-ring)", 0.0119, 0.0049, 54, 68.78, 34.40, 41.87, 33.12, 1.96E+1, -5.61E-3, 4.02E-5, -2.76E-8, 4.130, 6.817, std::nullopt, std::nullopt },
        CDJobackGroup { "-S- (ring)", 0.0019, 0.0051, 38, 52.10, 79.93, 39.10, 27.76, 1.67E+1, 4.81E-3, 2.77E-5, -2.11E-8, 1.557, 5.984, std::nullopt, std::nullopt }
    };

}    // namespace

namespace PCProps::ConstantData
{
    // ===== Constructor, default
    CDJoback::CDJoback() = default;

    // ===== Constructor, taking groups as arguments
    CDJoback::CDJoback(const std::vector<std::pair<int, int>>& groups, double molecularWeight, int atomCount) : m_molecularWeight(molecularWeight), m_atomCount(atomCount)
    {
        std::vector<std::optional<double>> itemsTc;
        std::vector<std::optional<double>> itemsPc;
        std::vector<std::optional<double>> itemsVc;
        std::vector<std::optional<double>> itemsTb;
        std::vector<std::optional<double>> itemsTm;
        std::vector<std::optional<double>> itemsHform;
        std::vector<std::optional<double>> itemsGform;
        std::vector<std::optional<double>> itemsIgCp_a;
        std::vector<std::optional<double>> itemsIgCp_b;
        std::vector<std::optional<double>> itemsIgCp_c;
        std::vector<std::optional<double>> itemsIgCp_d;
        std::vector<std::optional<double>> itemsHfus;
        std::vector<std::optional<double>> itemsHvap;
        std::vector<std::optional<double>> itemsLiqVis_a;
        std::vector<std::optional<double>> itemsLiqVis_b;

        for (const auto& item : groups) {
            itemsTc.emplace_back(JobackGroups.at(item.second - 1).tc.has_value() ? std::optional(static_cast<double>(item.first) * JobackGroups.at(item.second - 1).tc.value())
                                                                                 : std::nullopt);
            itemsPc.emplace_back(JobackGroups.at(item.second - 1).pc.has_value() ? std::optional(static_cast<double>(item.first) * JobackGroups.at(item.second - 1).pc.value())
                                                                                 : std::nullopt);
            itemsVc.emplace_back(JobackGroups.at(item.second - 1).vc.has_value() ? std::optional(static_cast<double>(item.first) * JobackGroups.at(item.second - 1).vc.value())
                                                                                 : std::nullopt);
            itemsTb.emplace_back(JobackGroups.at(item.second - 1).tb.has_value() ? std::optional(static_cast<double>(item.first) * JobackGroups.at(item.second - 1).tb.value())
                                                                                 : std::nullopt);
            itemsTm.emplace_back(JobackGroups.at(item.second - 1).tm.has_value() ? std::optional(static_cast<double>(item.first) * JobackGroups.at(item.second - 1).tm.value())
                                                                                 : std::nullopt);
            itemsHform.emplace_back(JobackGroups.at(item.second - 1).hform.has_value()
                                        ? std::optional(static_cast<double>(item.first) * JobackGroups.at(item.second - 1).hform.value())
                                        : std::nullopt);
            itemsGform.emplace_back(JobackGroups.at(item.second - 1).gform.has_value()
                                        ? std::optional(static_cast<double>(item.first) * JobackGroups.at(item.second - 1).gform.value())
                                        : std::nullopt);
            itemsIgCp_a.emplace_back(JobackGroups.at(item.second - 1).igCp_a.has_value()
                                         ? std::optional(static_cast<double>(item.first) * JobackGroups.at(item.second - 1).igCp_a.value())
                                         : std::nullopt);
            itemsIgCp_b.emplace_back(JobackGroups.at(item.second - 1).igCp_b.has_value()
                                         ? std::optional(static_cast<double>(item.first) * JobackGroups.at(item.second - 1).igCp_b.value())
                                         : std::nullopt);
            itemsIgCp_c.emplace_back(JobackGroups.at(item.second - 1).igCp_c.has_value()
                                         ? std::optional(static_cast<double>(item.first) * JobackGroups.at(item.second - 1).igCp_c.value())
                                         : std::nullopt);
            itemsIgCp_d.emplace_back(JobackGroups.at(item.second - 1).igCp_d.has_value()
                                         ? std::optional(static_cast<double>(item.first) * JobackGroups.at(item.second - 1).igCp_d.value())
                                         : std::nullopt);
            itemsHfus.emplace_back(
                JobackGroups.at(item.second - 1).hfus.has_value() ? std::optional(static_cast<double>(item.first) * JobackGroups.at(item.second - 1).hfus.value()) : std::nullopt);
            itemsHvap.emplace_back(
                JobackGroups.at(item.second - 1).hvap.has_value() ? std::optional(static_cast<double>(item.first) * JobackGroups.at(item.second - 1).hvap.value()) : std::nullopt);
            itemsLiqVis_a.emplace_back(JobackGroups.at(item.second - 1).liqVis_a.has_value()
                                           ? std::optional(static_cast<double>(item.first) * JobackGroups.at(item.second - 1).liqVis_a.value())
                                           : std::nullopt);
            itemsLiqVis_b.emplace_back(JobackGroups.at(item.second - 1).liqVis_b.has_value()
                                           ? std::optional(static_cast<double>(item.first) * JobackGroups.at(item.second - 1).liqVis_b.value())
                                           : std::nullopt);
        }

        m_sumTc = (std::find_if(itemsTc.begin(), itemsTc.end(), [](const std::optional<double>& item) { return !item.has_value(); }) == itemsTc.end()
                       ? std::optional(std::accumulate(itemsTc.begin(), itemsTc.end(), 0.0, [](double result, const std::optional<double>& item) { return result + item.value(); }))
                       : std::nullopt);

        m_sumPc = (std::find_if(itemsPc.begin(), itemsPc.end(), [](const std::optional<double>& item) { return !item.has_value(); }) == itemsPc.end()
                       ? std::optional(std::accumulate(itemsPc.begin(), itemsPc.end(), 0.0, [](double result, const std::optional<double>& item) { return result + item.value(); }))
                       : std::nullopt);

        m_sumVc = (std::find_if(itemsVc.begin(), itemsVc.end(), [](const std::optional<double>& item) { return !item.has_value(); }) == itemsVc.end()
                       ? std::optional(std::accumulate(itemsVc.begin(), itemsVc.end(), 0.0, [](double result, const std::optional<double>& item) { return result + item.value(); }))
                       : std::nullopt);

        m_sumTb = (std::find_if(itemsTb.begin(), itemsTb.end(), [](const std::optional<double>& item) { return !item.has_value(); }) == itemsTb.end()
                       ? std::optional(std::accumulate(itemsTb.begin(), itemsTb.end(), 0.0, [](double result, const std::optional<double>& item) { return result + item.value(); }))
                       : std::nullopt);

        m_sumTm = (std::find_if(itemsTm.begin(), itemsTm.end(), [](const std::optional<double>& item) { return !item.has_value(); }) == itemsTm.end()
                       ? std::optional(std::accumulate(itemsTm.begin(), itemsTm.end(), 0.0, [](double result, const std::optional<double>& item) { return result + item.value(); }))
                       : std::nullopt);

        m_sumHform =
            (std::find_if(itemsHform.begin(), itemsHform.end(), [](const std::optional<double>& item) { return !item.has_value(); }) == itemsHform.end()
                 ? std::optional(std::accumulate(itemsHform.begin(), itemsHform.end(), 0.0, [](double result, const std::optional<double>& item) { return result + item.value(); }))
                 : std::nullopt);

        m_sumGform =
            (std::find_if(itemsGform.begin(), itemsGform.end(), [](const std::optional<double>& item) { return !item.has_value(); }) == itemsGform.end()
                 ? std::optional(std::accumulate(itemsGform.begin(), itemsGform.end(), 0.0, [](double result, const std::optional<double>& item) { return result + item.value(); }))
                 : std::nullopt);

        m_sumIgCp_a =
            (std::find_if(itemsIgCp_a.begin(), itemsIgCp_a.end(), [](const std::optional<double>& item) { return !item.has_value(); }) == itemsIgCp_a.end()
                 ? std::optional(
                       std::accumulate(itemsIgCp_a.begin(), itemsIgCp_a.end(), 0.0, [](double result, const std::optional<double>& item) { return result + item.value(); }))
                 : std::nullopt);

        m_sumIgCp_b =
            (std::find_if(itemsIgCp_b.begin(), itemsIgCp_b.end(), [](const std::optional<double>& item) { return !item.has_value(); }) == itemsIgCp_b.end()
                 ? std::optional(
                       std::accumulate(itemsIgCp_b.begin(), itemsIgCp_b.end(), 0.0, [](double result, const std::optional<double>& item) { return result + item.value(); }))
                 : std::nullopt);

        m_sumIgCp_c =
            (std::find_if(itemsIgCp_c.begin(), itemsIgCp_c.end(), [](const std::optional<double>& item) { return !item.has_value(); }) == itemsIgCp_c.end()
                 ? std::optional(
                       std::accumulate(itemsIgCp_c.begin(), itemsIgCp_c.end(), 0.0, [](double result, const std::optional<double>& item) { return result + item.value(); }))
                 : std::nullopt);

        m_sumIgCp_d =
            (std::find_if(itemsIgCp_d.begin(), itemsIgCp_d.end(), [](const std::optional<double>& item) { return !item.has_value(); }) == itemsIgCp_d.end()
                 ? std::optional(
                       std::accumulate(itemsIgCp_d.begin(), itemsIgCp_d.end(), 0.0, [](double result, const std::optional<double>& item) { return result + item.value(); }))
                 : std::nullopt);

        m_sumHfus =
            (std::find_if(itemsHfus.begin(), itemsHfus.end(), [](const std::optional<double>& item) { return !item.has_value(); }) == itemsHfus.end()
                 ? std::optional(std::accumulate(itemsHfus.begin(), itemsHfus.end(), 0.0, [](double result, const std::optional<double>& item) { return result + item.value(); }))
                 : std::nullopt);

        m_sumHvap =
            (std::find_if(itemsHvap.begin(), itemsHvap.end(), [](const std::optional<double>& item) { return !item.has_value(); }) == itemsHvap.end()
                 ? std::optional(std::accumulate(itemsHvap.begin(), itemsHvap.end(), 0.0, [](double result, const std::optional<double>& item) { return result + item.value(); }))
                 : std::nullopt);

        m_sumLiqVis_a =
            (std::find_if(itemsLiqVis_a.begin(), itemsLiqVis_a.end(), [](const std::optional<double>& item) { return !item.has_value(); }) == itemsLiqVis_a.end()
                 ? std::optional(
                       std::accumulate(itemsLiqVis_a.begin(), itemsLiqVis_a.end(), 0.0, [](double result, const std::optional<double>& item) { return result + item.value(); }))
                 : std::nullopt);

        m_sumLiqVis_b =
            (std::find_if(itemsLiqVis_b.begin(), itemsLiqVis_b.end(), [](const std::optional<double>& item) { return !item.has_value(); }) == itemsLiqVis_b.end()
                 ? std::optional(
                       std::accumulate(itemsLiqVis_b.begin(), itemsLiqVis_b.end(), 0.0, [](double result, const std::optional<double>& item) { return result + item.value(); }))
                 : std::nullopt);

        for (const auto& group : groups) {
            m_groups.emplace_back(std::make_pair(group.first, JobackGroups.at(group.second - 1)));
        }
    }

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
        //        double result = 198.2;

        //        for (const auto& item : m_groups) result += item.second.tb.value() * item.first;

        return 198.2 + m_sumTb.value();
    }

    // ===== Compute melting temperature [K] for substance
    double CDJoback::meltingTemperature() const
    {
        //        double result = 122.5;
        //
        //        for (const auto& item : m_groups) {
        //            result += item.second.tm.value() * item.first;
        //        }

        return 122.5 + m_sumTm.value();
    }

    // ===== Compute critical temperature [K] for substance
    double CDJoback::criticalTemperature(double boilingTemperature) const
    {
        using std::pow;
        if (!boilingTemperature) boilingTemperature = this->boilingTemperature();
        return boilingTemperature / (0.584 + 0.965 * m_sumTc.value() - pow(m_sumTc.value(), 2));
    }

    // ===== Compute critical pressure [Pa] for substance
    double CDJoback::criticalPressure() const
    {
        using std::pow;
        //        double contrib = 0.0;
        //        for (const auto& item : m_groups) {
        //            contrib += item.second.pc.value() * item.first;
        //        }

        return pow(0.113 + 0.0032 * m_atomCount - m_sumPc.value(), -2) * 100000;
    }

    // ===== Compute critical volume [m3/mol] for substance
    double CDJoback::criticalVolume() const
    {
        //        double contrib = 17.5;
        //        for (const auto& item : m_groups) {
        //            contrib += item.second.vc.value() * item.first;
        //        }

        return (17.5 + m_sumVc.value()) / 1000000;
    }

    // ===== Compute heat of formation [J/mol] for substance
    double CDJoback::enthalpyOfFormation() const
    {
        //        double contrib = 68.29;
        //        for (const auto& item : m_groups) {
        //            contrib += item.second.hform.value() * item.first;
        //        }

        return (68.29 + m_sumHform.value()) * 1000;
    }

    // ===== Compute Gibbs energy of formation [J/mol] for substance
    double CDJoback::gibbsEnergyOfFormation() const
    {
        //        double contrib = 53.88;
        //        for (const auto& item : m_groups) {
        //            contrib += item.second.gform.value() * item.first;
        //        }

        return (53.88 + m_sumGform.value()) * 1000;
    }

    // ===== Compute heat of fusion [J/mol] for substance
    double CDJoback::enthalpyOfFusion() const
    {
        //        double contrib = -0.88;
        //        for (const auto& item : m_groups) {
        //            contrib += item.second.hfus.value() * item.first;
        //        }

        return (-0.88 + m_sumHfus.value()) * 1000;
    }

    // ===== Compute heat of vaporization [J/mol] for substance
    double CDJoback::enthalpyOfVaporization() const
    {
        //        double contrib = 15.30;
        //        for (const auto& item : m_groups) {
        //            contrib += item.second.hvap.value() * item.first;
        //        }

        return (15.30 + m_sumHvap.value()) * 1000;
    }

    // ===== Compute ideal gas Cp [J/(mol K)] for substance
    double CDJoback::idealGasCp(double temperature) const
    {
        //        double a = 0.0;
        //        for (const auto& item : m_groups) {
        //            a += item.second.igCp_a.value() * item.first;
        //        }
        //
        //        double b = 0.0;
        //        for (const auto& item : m_groups) {
        //            b += item.second.igCp_b.value() * item.first;
        //        }
        //
        //        double c = 0.0;
        //        for (const auto& item : m_groups) {
        //            c += item.second.igCp_c.value() * item.first;
        //        }
        //
        //        double d = 0.0;
        //        for (const auto& item : m_groups) {
        //            d += item.second.igCp_d.value() * item.first;
        //        }

        return m_sumIgCp_a.value() - 37.93 + (m_sumIgCp_b.value() + 0.21) * temperature + (m_sumIgCp_c.value() - 3.91E-4) * pow(temperature, 2) +
               (m_sumIgCp_d.value() + 2.06E-7) * pow(temperature, 3);
    }

    // ====== Compute liquid viscosity [Pa s] for substance.
    double CDJoback::liquidViscosity(double temperature) const
    {
        using std::exp;
        //
        //        double a = 0.0;
        //        for (const auto& item : m_groups) {
        //            a += item.second.liqVis_a.value() * item.first;
        //        }
        //
        //        double b = 0.0;
        //        for (const auto& item : m_groups) {
        //            b += item.second.liqVis_b.value() * item.first;
        //        }

        return m_molecularWeight * exp((m_sumLiqVis_a.value() - 597.82) / temperature + m_sumLiqVis_b.value() - 11.202);
    }
}    // namespace PCProps::ConstantData