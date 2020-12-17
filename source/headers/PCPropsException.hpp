//
// Created by Kenneth Balslev on 17/12/2020.
//

#ifndef PCPROPS_PCPROPSEXCEPTION_HPP
#define PCPROPS_PCPROPSEXCEPTION_HPP

// ===== External Includes ===== //
#include <stdexcept>

namespace PCProps
{
    class PCPropsException : public std::runtime_error
    {
    public:
        inline explicit PCPropsException(const std::string& err) : runtime_error(err) {};
    };

}

#endif    // PCPROPS_PCPROPSEXCEPTION_HPP
