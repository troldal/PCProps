//
// Created by Kenneth Balslev on 26/11/2021.
//

#include "DataSource.hpp"

#include <json/json.hpp>
#include <OpenXLSX.hpp>
#include <algorithm>

using namespace nlohmann;
using namespace OpenXLSX;

namespace {

    /**
     * @brief
     * @param record
     * @return
     */
    json extractConstants(const std::vector<XLCellValue>& record) {

        json result;
        result["Name"]                    = record[1].type() == XLValueType::Empty ? "" : record[1];
        result["CAS"]                     = record[2].type() == XLValueType::Empty ? "" : record[2];
        result["MolarWeight"]             = record[3].type() == XLValueType::Empty ? 0.0 : record[3].get<double>();
        result["NormalFreezingPoint"]      = record[4].type() == XLValueType::Empty ? 0.0 : record[4].get<double>();
        result["NormalBoilingPoint"]      = record[5].type() == XLValueType::Empty ? 0.0 : record[5].get<double>();
        result["CriticalTemperature"]     = record[6].type() == XLValueType::Empty ? 0.0 : record[6].get<double>();
        result["CriticalPressure"]        = record[7].type() == XLValueType::Empty ? 0.0 : record[7].get<double>();
        result["CriticalVolume"]          = record[8].type() == XLValueType::Empty ? 0.0 : record[8].get<double>();
        result["CriticalDensity"]         = record[9].type() == XLValueType::Empty ? 0.0 : record[9].get<double>();
        result["CriticalCompressibility"] = record[10].type() == XLValueType::Empty ? 0.0 : record[10].get<double>();
        result["AcentricFactor"]          = record[11].type() == XLValueType::Empty ? 0.0 : record[11].get<double>();
        result["DipoleMoment"]            = record[12].type() == XLValueType::Empty ? 0.0 : record[12].get<double>();

        return result;
    }

    /**
     * @brief
     * @param record
     * @return
     */
    json extractCorrelationData(const std::vector<XLCellValue>& record) {

        json result;
        result["Equation"] = record[3].type() == XLValueType::Empty ? "" : record[3];
        result["C1"] = record[4].type() == XLValueType::Empty ? 0.0 : record[4].get<double>();
        result["C2"] = record[5].type() == XLValueType::Empty ? 0.0 : record[5].get<double>();
        result["C3"] = record[6].type() == XLValueType::Empty ? 0.0 : record[6].get<double>();
        result["C4"] = record[7].type() == XLValueType::Empty ? 0.0 : record[7].get<double>();
        result["C5"] = record[8].type() == XLValueType::Empty ? 0.0 : record[8].get<double>();

        return result;
    }
}

namespace PCProps {

    /**
     *
     */
    class DataSource::impl
    {
    public:

        /**
         * @brief Default constructor.
         */
        impl() = default;

        /**
         * @brief Constructor taking path to Excel spreadsheet as argument.
         * @param filename Path to Excel spreadsheet.
         */
        explicit impl(const std::string& filename) : m_filename(filename) {}

        /**
         * @brief Destructor;
         */
        ~impl() = default;

        /**
         *
         * @param other
         */
        impl(const impl& other) = default;

        /**
         *
         * @param other
         */
        impl(impl&& other) noexcept = default;

        /**
         *
         * @param other
         * @return
         */
        impl& operator=(const impl& other) = default;

        /**
         *
         * @param other
         * @return
         */
        impl& operator=(impl&& other) noexcept = default;

        /**
         *
         * @return
         */
        json load() const {

            XLDocument doc;
            doc.open(m_filename);
            auto wbk = doc.workbook();

            // ===== Load pure component constants (e.g. critical properties)
            auto wks = wbk.worksheet("Constants");
            std::vector<json> objects;
            for (const auto &row: wks.rows(2, wks.rowCount()))
                objects.emplace_back(extractConstants(row.values()));

            // ===== Load ideal gas heat capacity correlation coefficients
            wks = wbk.worksheet("IdealGasCp");
            for (const auto &row: wks.rows(2, wks.rowCount())) {
                std::vector<XLCellValue> values = row.values();
                auto& object = *std::find_if(objects.begin(), objects.end(), [&](const json& obj){return obj["CAS"].get<std::string>() == values[2].get<std::string>();} );
                if (object == *objects.end()) continue;
                object["IdealGasCp"] = extractCorrelationData(values);
            }

            // ===== Load liquid heat capacity correlation coefficients
            wks = wbk.worksheet("LiquidCp");
            for (const auto &row: wks.rows(2, wks.rowCount())) {
                std::vector<XLCellValue> values = row.values();
                auto& object = *std::find_if(objects.begin(), objects.end(), [&](const json& obj){return obj["CAS"].get<std::string>() == values[2].get<std::string>();} );
                if (object == *objects.end()) continue;
                object["LiquidCp"] = extractCorrelationData(values);
            }

            // ===== Load vapor pressure correlation coefficients
            wks = wbk.worksheet("VaporPressure");
            for (const auto &row: wks.rows(2, wks.rowCount())) {
                std::vector<XLCellValue> values = row.values();
                auto& object = *std::find_if(objects.begin(), objects.end(), [&](const json& obj){return obj["CAS"].get<std::string>() == values[2].get<std::string>();} );
                if (object == *objects.end()) continue;
                object["VaporPressure"] = extractCorrelationData(values);
            }

            // ===== Load saturated vapor viscosity correlation coefficients
            wks = wbk.worksheet("SaturatedVaporViscosity");
            for (const auto &row: wks.rows(2, wks.rowCount())) {
                std::vector<XLCellValue> values = row.values();
                auto& object = *std::find_if(objects.begin(), objects.end(), [&](const json& obj){return obj["CAS"].get<std::string>() == values[2].get<std::string>();} );
                if (object == *objects.end()) continue;
                object["SaturatedVaporViscosity"] = extractCorrelationData(values);
            }

            // ===== Load saturated liquid viscosity correlation coefficients
            wks = wbk.worksheet("SaturatedLiquidViscosity");
            for (const auto &row: wks.rows(2, wks.rowCount())) {
                std::vector<XLCellValue> values = row.values();
                auto& object = *std::find_if(objects.begin(), objects.end(), [&](const json& obj){return obj["CAS"].get<std::string>() == values[2].get<std::string>();} );
                if (object == *objects.end()) continue;
                object["SaturatedLiquidViscosity"] = extractCorrelationData(values);
            }

            // ===== Load saturated liquid volume (density) correlation coefficients
            wks = wbk.worksheet("SaturatedLiquidVolume");
            for (const auto &row: wks.rows(2, wks.rowCount())) {
                std::vector<XLCellValue> values = row.values();
                auto& object = *std::find_if(objects.begin(), objects.end(), [&](const json& obj){return obj["CAS"].get<std::string>() == values[2].get<std::string>();} );
                if (object == *objects.end()) continue;
                object["SaturatedLiquidVolume"] = extractCorrelationData(values);
            }

            // ===== Load vapor thermal conductivity correlation coefficients
            wks = wbk.worksheet("VaporThermalConductivity");
            for (const auto &row: wks.rows(2, wks.rowCount())) {
                std::vector<XLCellValue> values = row.values();
                auto& object = *std::find_if(objects.begin(), objects.end(), [&](const json& obj){return obj["CAS"].get<std::string>() == values[2].get<std::string>();} );
                if (object == *objects.end()) continue;
                object["VaporThermalConductivity"] = extractCorrelationData(values);
            }

            // ===== Load liquid thermal conductivity correlation coefficients
            wks = wbk.worksheet("LiquidThermalConductivity");
            for (const auto &row: wks.rows(2, wks.rowCount())) {
                std::vector<XLCellValue> values = row.values();
                auto& object = *std::find_if(objects.begin(), objects.end(), [&](const json& obj){return obj["CAS"].get<std::string>() == values[2].get<std::string>();} );
                if (object == *objects.end()) continue;
                object["LiquidThermalConductivity"] = extractCorrelationData(values);
            }

            doc.close();
            return json(objects); // NOLINT
        }

    private:
        std::string m_filename;
    };

    // =====
    DataSource::DataSource() : m_impl(nullptr) {}

    // =====
    DataSource::DataSource(const std::string& filename) : m_impl(std::make_unique<impl>(filename)) {}

    // =====
    DataSource::~DataSource() = default;

    // =====
    DataSource::DataSource(const DataSource& other) : m_impl(other.m_impl ? std::make_unique<impl>(*other.m_impl) : nullptr) {}

    // =====
    DataSource::DataSource(DataSource&& other) noexcept = default;

    // =====
    DataSource& DataSource::operator=(const DataSource& other)
    {
        DataSource copy = other;
        *this             = std::move(copy);
        return *this;
    }

    // =====
    DataSource& DataSource::operator=(DataSource&& other) noexcept = default;

    // =====
    DataSource::JSONString DataSource::load() const
    {
        return m_impl->load().dump();
    }

} // namespace PCProps
