//
// Created by Kenneth Balslev on 20/01/2021.
//

#include "Pipe.hpp"

namespace PCProps::UnitOps {
    Pipe::Pipe() = default;

    Pipe::Pipe(double length, double diameter, double inclination, double roughness) {
        m_pipeSegments.emplace_back(PipeSegmentLiquid(length, diameter, inclination, roughness));
    }

    double Pipe::computeOutletPressure(const IFluid& inletFluid, double molarFlow)
    {
        return m_pipeSegments[0].computeOutletPressure(inletFluid, molarFlow);
    }

    double Pipe::computeInletPressure(const IFluid& outletFluid, double molarFlow)
    {
        return m_pipeSegments[0].computeInletPressure(outletFluid, molarFlow);;
    }

    double Pipe::computeMolarFlow(const IFluid& inletFluid, const IFluid& outletFluid)
    {
        return 0;
    }

}