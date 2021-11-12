//
// Created by Kenneth Balslev on 19/01/2021.
//

#ifndef PCPROPS_PIPESEGMENTLIQUID_HPP
#define PCPROPS_PIPESEGMENTLIQUID_HPP

#include <common/PropertyData.hpp>
#include <IFluid.hpp>

#include <memory>

namespace PCProps::UnitOps
{
    class PipeSegmentLiquid
    {
    public:

        // =====================================================================
        // CONSTRUCTORS & ASSIGNMENT OPERATORS
        // =====================================================================

        PipeSegmentLiquid();

        PipeSegmentLiquid(double length, double diameter, double inclination, double roughness);

        PipeSegmentLiquid(const PipeSegmentLiquid& other);

        PipeSegmentLiquid(PipeSegmentLiquid&& other) noexcept;

        ~PipeSegmentLiquid();

        PipeSegmentLiquid& operator=(const PipeSegmentLiquid& other);

        PipeSegmentLiquid& operator=(PipeSegmentLiquid&& other) noexcept;


        // =====================================================================
        // FLOW CALCULATIONS
        // =====================================================================

        double computeOutletPressure(const IFluid& inletFluid, double molarFlow) const;

        double computeInletPressure(const IFluid& outletFluid, double molarFlow) const;

        double computeMolarFlow(const IFluid& inletFluid, const IFluid& outletFluid) const;

    private:
        class impl;
        std::unique_ptr<impl> m_impl;

    };
}    // namespace PCProps::UnitOps

#endif    // PCPROPS_PIPESEGMENTLIQUID_HPP
