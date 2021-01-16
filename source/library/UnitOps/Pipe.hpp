//
// Created by Kenneth Balslev on 15/01/2021.
//

#ifndef PCPROPS_PIPE_HPP
#define PCPROPS_PIPE_HPP

#include <PCPropsData.hpp>

namespace PCProps::UnitOps
{
    class Pipe
    {

        double m_length {};
        double m_diameter {};
        double m_inclination {};
        double m_roughness {};

    public:

        Pipe();

        Pipe(double length, double diameter, double inclination, double roughness);

        double computeOutletPressure(const PCPhases& fluid, double molarFlow);


    };
}    // namespace PCProps::UnitOps

#endif    // PCPROPS_PIPE_HPP
