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
         * @param pressure
         * @param temperature
         * @return
         */
        JSONString flashPT(double pressure, double temperature) const;

        /**
         * @brief
         * @param pressure
         * @param vaporFraction
         * @return
         */
        JSONString flashPx(double pressure, double vaporFraction) const;

        /**
         * @brief
         * @param temperature
         * @param vaporFraction
         * @return
         */
        JSONString flashTx(double temperature, double vaporFraction) const;

        /**
         * @brief
         * @param pressure
         * @param enthalpy
         * @return
         */
        JSONString flashPH(double pressure, double enthalpy) const;

        /**
         * @brief
         * @param pressure
         * @param entropy
         * @return
         */
        JSONString flashPS(double pressure, double entropy) const;

        /**
         * @brief
         * @param temperature
         * @param volume
         * @return
         */
        JSONString flashTV(double temperature, double volume) const;

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
