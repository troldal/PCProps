//
// Created by Kenneth Balslev on 19/01/2021.
//

#include "PipeSegmentLiquid.hpp"

#include <common/Globals.hpp>

#include <cmath>
#include <stdexcept>

namespace PCProps::UnitOps
{

    class PipeSegmentLiquid::impl {

        double m_length {};
        double m_diameter {};
        double m_inclination {};
        double m_roughness {};
        double m_flowArea {};

        double computeElevationGradient(const PCPhases& fluid) const {
            return -(1.0 / fluid[0][PCMolarVolume]) * (fluid[0][PCMolarWeight] / 1000.0) * Globals::G_ACCL * sin(m_inclination * 0.01745329252);
        }

        double computeVolumeFlow(const PCPhases& fluid, double molarFlow) const {
            return molarFlow * fluid[0][PCMolarVolume];
        }

        double computeReynoldsNumber(const PCPhases& fluid, double velocity) const {
           return (1.0 / fluid[0][PCMolarVolume]) * (fluid[0][PCMolarWeight] / 1000.0) * velocity * m_diameter / fluid[0][PCViscosity];
        }

        double computeFrictionFactor(double reynoldsNumber) const {

            if (m_roughness == 0.0) return 0.184 * pow(reynoldsNumber, -0.2);

            using std::pow;
            using std::log;
            double a = 1.0 / (1.0 + pow(reynoldsNumber/2712, 8.4));
            double b = 1.0 / (1.0 + pow(reynoldsNumber/(150 * m_diameter / m_roughness), 1.8));

            return pow(64/reynoldsNumber, a) *
                   pow(0.75 * log(reynoldsNumber/5.37), 2 * (a - 1) * b) *
                   pow(0.88 * log(3.41 * m_diameter / m_roughness), 2 * (a-1) * (1-b));
        }


    public:

        impl(double length, double diameter, double inclination, double roughness)
            : m_length{length},
              m_diameter{diameter},
              m_inclination{inclination},
              m_roughness{roughness},
              m_flowArea{Globals::PI * pow(m_diameter, 2) / 4.0}
        {}

        double computeOutletPressure(const IFluid& fluid, double molarFlow) const {
            using std::sin;
            using std::pow;
            if (fluid.properties().size() > 1) throw std::invalid_argument("Invalid fluid");
            double elevationGradient = computeElevationGradient(fluid.properties());
            double volumeFlow        = computeVolumeFlow(fluid.properties(), molarFlow);
            double velocity          = volumeFlow / m_flowArea;
            double reynoldsNumber    = computeReynoldsNumber(fluid.properties(), velocity);
            double dpdl = -computeFrictionFactor(reynoldsNumber) *
                          (1.0 / fluid.properties()[0][PCMolarVolume]) *
                          (fluid.properties()[0][PCMolarWeight] / 1000.0) * pow(velocity, 2) / (2 * m_diameter);

            return fluid.properties()[0][PCPressure] + (dpdl + elevationGradient) * m_length;
        }

        double computeInletPressure(const IFluid& fluid, double molarFlow) {
            using std::sin;
            using std::pow;
            if (fluid.properties().size() > 1) throw std::invalid_argument("Invalid fluid");
            double elevationGradient = computeElevationGradient(fluid.properties());
            double volumeFlow        = computeVolumeFlow(fluid.properties(), molarFlow);
            double velocity          = volumeFlow / m_flowArea;
            double reynoldsNumber    = computeReynoldsNumber(fluid.properties(), velocity);
            double dpdl = -computeFrictionFactor(reynoldsNumber) *
                (1.0 / fluid.properties()[0][PCMolarVolume]) *
                (fluid.properties()[0][PCMolarWeight] / 1000.0) * pow(velocity, 2) / (2 * m_diameter);

            return fluid.properties()[0][PCPressure] - (dpdl + elevationGradient) * m_length;
        }

    };

    PipeSegmentLiquid::PipeSegmentLiquid() = default;

    PipeSegmentLiquid::PipeSegmentLiquid(double length, double diameter, double inclination, double roughness)
        : m_impl(std::make_unique<impl>(length, diameter, inclination, roughness))
    {}

    PipeSegmentLiquid::PipeSegmentLiquid(const PipeSegmentLiquid& other) : m_impl(std::make_unique<impl>(*other.m_impl)) {};

    PipeSegmentLiquid::PipeSegmentLiquid(PipeSegmentLiquid&& other) noexcept = default;

    PipeSegmentLiquid::~PipeSegmentLiquid() = default;

    PipeSegmentLiquid& PipeSegmentLiquid::operator=(const PipeSegmentLiquid& other)
    {
        PipeSegmentLiquid copy = other;
        *this                  = std::move(copy);
        return *this;
    }

    PipeSegmentLiquid& PipeSegmentLiquid::operator=(PipeSegmentLiquid&& other) noexcept = default;

    double PipeSegmentLiquid::computeOutletPressure(const IFluid& inletFluid, double molarFlow) const
    {
        return m_impl->computeOutletPressure(inletFluid, molarFlow);
    }

    double PipeSegmentLiquid::computeInletPressure(const IFluid& outletFluid, double molarFlow) const
    {
        return m_impl->computeInletPressure(outletFluid, molarFlow);
    }

    double PipeSegmentLiquid::computeMolarFlow(const IFluid& inletFluid, const IFluid& outletFluid) const
    {
        return 0;
    }

}  // PCProps::UnitOps