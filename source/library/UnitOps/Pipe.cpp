//
// Created by Kenneth Balslev on 15/01/2021.
//

#include "Pipe.hpp"

#include <library/PCGlobals.hpp>
#include <library/PCPropsException.hpp>

#include <cmath>

namespace PCProps::UnitOps
{
    Pipe::Pipe() = default;

    Pipe::Pipe(double length, double diameter, double inclination, double roughness)
        : m_length(length),
          m_diameter(diameter),
          m_inclination(inclination),
          m_roughness(roughness)
    {}

    double Pipe::computeOutletPressure(const PCPhases& fluid, double molarFlow)
    {
        using std::sin;
        using std::pow;
        if (fluid.size() > 1) throw PCPropsException("Invalid fluid");
        double elevationGradient = -(1.0 / fluid[0][PCMolarVolume]) * (fluid[0][PCMolarWeight] / 1000.0) * Globals::G_ACCL * sin(m_inclination * 0.01745329252);
        double flowArea = Globals::PI * pow(m_diameter,2) / 4.0;
        double volumeFlow = molarFlow * fluid[0][PCMolarVolume];
        double velocity = volumeFlow / flowArea;
        double reynoldsNumber = (1.0 / fluid[0][PCMolarVolume]) * (fluid[0][PCMolarWeight] / 1000.0) * velocity * m_diameter / fluid[0][PCViscosity];

        return reynoldsNumber;
    }

} // namespace PCProps::UnitOps