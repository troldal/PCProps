//
// Created by Kenneth Balslev on 26/11/2021.
//

#include "DataSource.hpp"

#include <json/json.hpp>
#include <OpenXLSX.hpp>


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
