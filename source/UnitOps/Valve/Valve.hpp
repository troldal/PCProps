//
// Created by Kenneth Balslev on 23/01/2022.
//

#ifndef PCPROPS_VALVE_HPP
#define PCPROPS_VALVE_HPP

#include <Stream/Stream.hpp>

#include <string>

namespace PCProps::UnitOps
{

    class Valve
    {

        using JSONString = std::string;

    public:

        /**
         * @brief
         * @param specification
         */
        Valve(const JSONString& specification);

        /**
         * @brief
         * @param other
         */
        Valve(const Valve& other);

        /**
         * @brief
         * @param other
         */
        Valve(Valve&& other) noexcept;

        /**
         * @brief
         */
        ~Valve();

        /**
         * @brief
         * @param other
         * @return
         */
        Valve& operator=(const Valve& other);

        /**
         * @brief
         * @param other
         * @return
         */
        Valve& operator=(Valve&& other) noexcept;

        /**
         * @brief
         * @param stream
         */
        void setInletStreams(const std::vector<const Stream*> streams);

        /**
         * @brief
         * @param specification
         */
        void setSpecification(const JSONString& specification);

        /**
         * @brief
         * @param streamName
         * @return
         */
        const std::vector<Stream>& outletStreams();

        /**
         * @brief
         */
        void compute();

        /**
         * @brief
         * @return
         */
        JSONString results() const;

    private:

        class impl;
        std::unique_ptr<impl> m_impl;

    };
}    // namespace PCProps::UnitOps

#endif    // PCPROPS_VALVE_HPP
