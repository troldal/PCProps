//
// Created by Kenneth Balslev on 20/01/2021.
//

#ifndef PCPROPS_PIPE_HPP
#define PCPROPS_PIPE_HPP

#include <vector>

#include "PipeSegmentLiquid.hpp"
#include <IFluid.hpp>

namespace PCProps::UnitOps
{
    class Pipe
    {
        double m_length {};
        double m_diameter {};
        double m_inclination {};
        double m_roughness {};

        std::vector<PipeSegmentLiquid> m_pipeSegments {};

    public:

        Pipe();

        Pipe(double length, double diameter, double inclination, double roughness);

        double computeOutletPressure(const IFluid& inletFluid, double molarFlow);

        double computeInletPressure(const IFluid& outletFluid, double molarFlow);

        double computeMolarFlow(const IFluid& inletFluid, const IFluid& outletFluid);


    };
}    // namespace PCProps::UnitOps

#endif    // PCPROPS_PIPE_HPP
