//
// Created by Kenneth Balslev on 26/11/2021.
//

#ifndef PCPROPS_DATASOURCE_HPP
#define PCPROPS_DATASOURCE_HPP

#include <string>
#include <memory>

namespace PCProps
{
    /**
     *
     */
    class DataSource
    {
        using JSONString = std::string;
    public:

        /**
         *
         */
        DataSource();

        /**
         *
         * @param filename
         */
        explicit DataSource(const std::string& filename);

        /**
         *
         */
        ~DataSource();

        /**
         *
         * @param other
         */
        DataSource(const DataSource& other);

        /**
         *
         * @param other
         */
        DataSource(DataSource&& other) noexcept;

        /**
         *
         * @param other
         * @return
         */
        DataSource& operator=(const DataSource& other);

        /**
         *
         * @param other
         * @return
         */
        DataSource& operator=(DataSource&& other) noexcept;

        /**
         *
         * @return
         */
        JSONString load() const;

    private:
        class impl;
        std::unique_ptr<impl> m_impl;
    };

} // namespace PCProps

#endif    // PCPROPS_DATASOURCE_HPP
