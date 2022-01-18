//
// Created by Kenneth Balslev on 20/01/2021.
//

#include "Pipe.hpp"

namespace PCProps::UnitOps {
    Pipe::Pipe() = default;

    Pipe::Pipe(double length, double diameter, double inclination, double roughness) {
        m_pipeSegments.emplace_back(PipeSegmentLiquid(length, diameter, inclination, roughness));
    }

    double Pipe::computeOutletPressure(const IPropertyPackage& inletFluid, double molarFlow)
    {
        return m_pipeSegments[0].computeOutletPressure(inletFluid, molarFlow);
    }

    double Pipe::computeInletPressure(const IPropertyPackage& outletFluid, double molarFlow)
    {
        return m_pipeSegments[0].computeInletPressure(outletFluid, molarFlow);;
    }

    double Pipe::computeMolarFlow(const IPropertyPackage& inletFluid, const IPropertyPackage& outletFluid)
    {
        return 0;
    }

}