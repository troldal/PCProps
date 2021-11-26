//
// Created by Kenneth Balslev on 25/11/2021.
//

#ifndef PCPROPS_PURECOMPONENTFACTORY_HPP
#define PCPROPS_PURECOMPONENTFACTORY_HPP

#include "PureComponent.hpp"

#include <string>
#include <memory>

namespace PCProps
{

    class PureComponentFactory
    {
        using JSONString = std::string;
    public:
        PureComponentFactory();

        PureComponentFactory(const JSONString& pcdata);

        ~PureComponentFactory();

        PureComponentFactory(const PureComponentFactory& other);

        PureComponentFactory(PureComponentFactory&& other) noexcept;

        PureComponentFactory& operator=(const PureComponentFactory& other);

        PureComponentFactory& operator=(PureComponentFactory&& other) noexcept;

        void init(const JSONString& pcdata);

        PureComponent makeComponent(const std::string& CAS) const;


    private:

        class impl;
        std::unique_ptr<impl> m_impl;



    };

} // namespace PCProps

#endif    // PCPROPS_PURECOMPONENTFACTORY_HPP
