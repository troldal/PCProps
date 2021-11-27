//
// Created by Kenneth Balslev on 26/11/2021.
//

#include "DataSource.hpp"

#include <json/json.hpp>
#include <OpenXLSX.hpp>
#include <algorithm>


using namespace nlohmann;
using namespace OpenXLSX;

namespace PCProps {

    /**
     *
     */
    class DataSource::impl
    {
    public:

        /**
         *
         */
        impl() = default;

        /**
         *
         * @param filename
         */
        explicit impl(const std::string& filename) : m_filename(filename) {}

        /**
         *
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
            auto wks = wbk.worksheet("Constants");

            std::vector<json> objects;
            for (const auto &row: wks.rows(2, wks.rowCount())) {
                std::vector<XLCellValue> values = row.values();
                json object;
                object["Name"]                    = values[1].type() == XLValueType::Empty ? "" : values[1];
                object["CAS"]                     = values[2].type() == XLValueType::Empty ? "" : values[2];
                object["MolarWeight"]             = values[4].type() == XLValueType::Empty ? 0.0 : values[4].get<double>();
                object["MeltingTemperature"]      = values[13].type() == XLValueType::Empty ? 0.0 : values[13].get<double>();
                object["BoilingTemperature"]      = values[5].type() == XLValueType::Empty ? 0.0 : values[5].get<double>();
                object["CriticalTemperature"]     = values[6].type() == XLValueType::Empty ? 0.0 : values[6].get<double>();
                object["CriticalPressure"]        = values[7].type() == XLValueType::Empty ? 0.0 : values[7].get<double>();
                object["CriticalVolume"]          = values[8].type() == XLValueType::Empty ? 0.0 : values[8].get<double>();
                object["CriticalDensity"]         = values[9].type() == XLValueType::Empty ? 0.0 : values[9].get<double>();
                object["CriticalCompressibility"] = values[10].type() == XLValueType::Empty ? 0.0 : values[10].get<double>();
                object["AcentricFactor"]          = values[11].type() == XLValueType::Empty ? 0.0 : values[11].get<double>();
                object["DipoleMoment"]            = values[12].type() == XLValueType::Empty ? 0.0 : values[12].get<double>();

                objects.emplace_back(object);
            }

            wks = wbk.worksheet("Ideal Gas Cp");
            for (const auto &row: wks.rows(2, wks.rowCount())) {
                std::vector<XLCellValue> values = row.values();
                auto& object = *std::find_if(objects.begin(), objects.end(), [&](const json& obj){return obj["CAS"].get<std::string>() == values[2].get<std::string>();} );
                json igcp;
                igcp["Equation"] = values[4].type() == XLValueType::Empty ? "" : values[4];
                igcp["C1"] = values[5].type() == XLValueType::Empty ? 0.0 : values[5].get<double>();
                igcp["C2"] = values[6].type() == XLValueType::Empty ? 0.0 : values[6].get<double>();
                igcp["C3"] = values[7].type() == XLValueType::Empty ? 0.0 : values[7].get<double>();
                igcp["C4"] = values[8].type() == XLValueType::Empty ? 0.0 : values[8].get<double>();
                igcp["C5"] = values[9].type() == XLValueType::Empty ? 0.0 : values[9].get<double>();

                object["IdealGasCp"] = igcp;
            }

            wks = wbk.worksheet("Liquid Cp");
            for (const auto &row: wks.rows(2, wks.rowCount())) {
                std::vector<XLCellValue> values = row.values();
                auto& object = *std::find_if(objects.begin(), objects.end(), [&](const json& obj){return obj["CAS"].get<std::string>() == values[2].get<std::string>();} );
                json liqcp;
                liqcp["Equation"] = values[4].type() == XLValueType::Empty ? "" : values[4];
                liqcp["C1"] = values[5].type() == XLValueType::Empty ? 0.0 : values[5].get<double>();
                liqcp["C2"] = values[6].type() == XLValueType::Empty ? 0.0 : values[6].get<double>();
                liqcp["C3"] = values[7].type() == XLValueType::Empty ? 0.0 : values[7].get<double>();
                liqcp["C4"] = values[8].type() == XLValueType::Empty ? 0.0 : values[8].get<double>();
                liqcp["C5"] = values[9].type() == XLValueType::Empty ? 0.0 : values[9].get<double>();

                object["LiquidCp"] = liqcp;
            }

            wks = wbk.worksheet("Vapor Pressure");
            for (const auto &row: wks.rows(2, wks.rowCount())) {
                std::vector<XLCellValue> values = row.values();
                auto& object = *std::find_if(objects.begin(), objects.end(), [&](const json& obj){return obj["CAS"].get<std::string>() == values[2].get<std::string>();} );
                json vp;
                vp["Equation"] = values[5].type() == XLValueType::Empty ? "" : values[5];
                vp["C1"] = values[6].type() == XLValueType::Empty ? 0.0 : values[6].get<double>();
                vp["C2"] = values[7].type() == XLValueType::Empty ? 0.0 : values[7].get<double>();
                vp["C3"] = values[8].type() == XLValueType::Empty ? 0.0 : values[8].get<double>();
                vp["C4"] = values[9].type() == XLValueType::Empty ? 0.0 : values[9].get<double>();
                vp["C5"] = values[10].type() == XLValueType::Empty ? 0.0 : values[10].get<double>();

                object["VaporPressure"] = vp;
            }

            wks = wbk.worksheet("Sat Vap Visc");
            for (const auto &row: wks.rows(2, wks.rowCount())) {
                std::vector<XLCellValue> values = row.values();
                auto& object = *std::find_if(objects.begin(), objects.end(), [&](const json& obj){return obj["CAS"].get<std::string>() == values[2].get<std::string>();} );
                json svv;
                svv["Equation"] = values[4].type() == XLValueType::Empty ? "" : values[4];
                svv["C1"] = values[5].type() == XLValueType::Empty ? 0.0 : values[5].get<double>();
                svv["C2"] = values[6].type() == XLValueType::Empty ? 0.0 : values[6].get<double>();
                svv["C3"] = values[7].type() == XLValueType::Empty ? 0.0 : values[7].get<double>();
                svv["C4"] = values[8].type() == XLValueType::Empty ? 0.0 : values[8].get<double>();

                object["SaturatedVaporViscosity"] = svv;
            }

            wks = wbk.worksheet("Sat Liq Visc");
            for (const auto &row: wks.rows(2, wks.rowCount())) {
                std::vector<XLCellValue> values = row.values();
                auto& object = *std::find_if(objects.begin(), objects.end(), [&](const json& obj){return obj["CAS"].get<std::string>() == values[2].get<std::string>();} );
                json slv;
                slv["Equation"] = values[4].type() == XLValueType::Empty ? "" : values[4];
                slv["C1"] = values[5].type() == XLValueType::Empty ? 0.0 : values[5].get<double>();
                slv["C2"] = values[6].type() == XLValueType::Empty ? 0.0 : values[6].get<double>();
                slv["C3"] = values[7].type() == XLValueType::Empty ? 0.0 : values[7].get<double>();
                slv["C4"] = values[8].type() == XLValueType::Empty ? 0.0 : values[8].get<double>();
                slv["C5"] = values[9].type() == XLValueType::Empty ? 0.0 : values[9].get<double>();

                object["SaturatedLiquidViscosity"] = slv;
            }

            wks = wbk.worksheet("Sat Liq Dens");
            for (const auto &row: wks.rows(2, wks.rowCount())) {
                std::vector<XLCellValue> values = row.values();
                auto& object = *std::find_if(objects.begin(), objects.end(), [&](const json& obj){return obj["CAS"].get<std::string>() == values[2].get<std::string>();} );
                json sld;
                sld["Equation"] = values[4].type() == XLValueType::Empty ? "" : values[4];
                sld["C1"] = values[5].type() == XLValueType::Empty ? 0.0 : values[5].get<double>();
                sld["C2"] = values[6].type() == XLValueType::Empty ? 0.0 : values[6].get<double>();
                sld["C3"] = values[7].type() == XLValueType::Empty ? 0.0 : values[7].get<double>();
                sld["C4"] = values[8].type() == XLValueType::Empty ? 0.0 : values[8].get<double>();

                object["SaturatedLiquidVolume"] = sld;
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
