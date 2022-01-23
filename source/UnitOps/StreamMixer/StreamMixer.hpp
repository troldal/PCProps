//
// Created by Kenneth Balslev on 22/01/2022.
//

#ifndef PCPROPS_STREAMMIXER_HPP
#define PCPROPS_STREAMMIXER_HPP

#include <Stream/Stream.hpp>

#include <string>

namespace PCProps::UnitOps
{

    class StreamMixer
    {

        using JSONString = std::string;

    public:

        /**
         * @brief
         * @param specification
         */
        StreamMixer(const JSONString& specification);

        /**
         * @brief
         * @param other
         */
        StreamMixer(const StreamMixer& other);

        /**
         * @brief
         * @param other
         */
        StreamMixer(StreamMixer&& other) noexcept;

        /**
         * @brief
         */
        ~StreamMixer();

        /**
         * @brief
         * @param other
         * @return
         */
        StreamMixer& operator=(const StreamMixer& other);

        /**
         * @brief
         * @param other
         * @return
         */
        StreamMixer& operator=(StreamMixer&& other) noexcept;

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

#endif    // PCPROPS_STREAMMIXER_HPP
