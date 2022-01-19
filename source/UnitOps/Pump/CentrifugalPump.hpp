//
// Created by Kenneth Balslev on 22/01/2021.
//

#ifndef PCPROPS_CENTRIFUGALPUMP_HPP
#define PCPROPS_CENTRIFUGALPUMP_HPP

#include <IPropertyPackage.hpp>
#include <Stream/Stream.hpp>

#include <array>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace PCProps::UnitOps
{

    class CentrifugalPump
    {

        using JSONString = std::string;

    public:

        /**
         * @brief
         * @param specification
         */
        CentrifugalPump(const JSONString& specification);

        /**
         * @brief
         * @param other
         */
        CentrifugalPump(const CentrifugalPump& other);

        /**
         * @brief
         * @param other
         */
        CentrifugalPump(CentrifugalPump&& other) noexcept;

        /**
         * @brief
         */
        ~CentrifugalPump();

        /**
         * @brief
         * @param other
         * @return
         */
        CentrifugalPump& operator=(const CentrifugalPump& other);

        /**
         * @brief
         * @param other
         * @return
         */
        CentrifugalPump& operator=(CentrifugalPump&& other) noexcept;

        /**
         * @brief
         * @return
         */
        const Stream& operator()() const;

        /**
         * @brief
         * @return
         */
        std::string results() const;

        /**
         * @brief
         * @param stream
         */
        void setInletStream(Stream* stream);

        /**
         * @brief
         * @param specification
         */
        void setSpecification(const JSONString& specification);


    private:

        class impl;
        std::unique_ptr<impl> m_impl;

    };
}    // namespace PCProps::UnitOps

#endif    // PCPROPS_CENTRIFUGALPUMP_HPP
