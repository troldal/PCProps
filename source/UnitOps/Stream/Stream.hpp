//
// Created by Kenneth Balslev on 24/01/2021.
//

#ifndef PCPROPS_STREAM_HPP
#define PCPROPS_STREAM_HPP

#include <IPropertyPackage.hpp>

#include <memory>
#include <string>

namespace PCProps::UnitOps
{
    /**
     * @brief
     */
    class Stream
    {
        using JSONString = std::string;

    public:

        /**
         * @brief
         */
        Stream();

        /**
         * @brief
         * @param fluid
         * @param quantity
         */
        Stream(const IPropertyPackage& fluid, double quantity);

        /**
         * @brief
         * @param other
         */
        Stream(const Stream& other);

        /**
         * @brief
         * @param other
         */
        Stream(Stream&& other) noexcept;

        /**
         * @brief
         */
        ~Stream();

        /**
         * @brief
         * @param other
         * @return
         */
        Stream& operator=(const Stream& other);

        /**
         * @brief
         * @param other
         * @return
         */
        Stream& operator=(Stream&& other) noexcept;

        /**
         * @brief
         * @param spec
         * @param s1
         * @param s2
         * @return
         */
        JSONString flash(const std::string& spec, double s1, double s2);

        /**
         * @brief
         * @return
         */
        JSONString properties() const;

    private:
        class impl;
        std::unique_ptr<impl> m_impl;

    };
}    // namespace PCProps::UnitOps

#endif    // PCPROPS_STREAM_HPP
