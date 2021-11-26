//
// Created by Kenneth Balslev on 19/01/2021.
//

#ifndef PCPROPS_PROPERTYLIB_HPP
#define PCPROPS_PROPERTYLIB_HPP

#include "../PureComponent/PureComponent.hpp"
#include "ConstantData/Joback.hpp"

#include "CompressedLiquidVolume/Aalto.hpp"
#include "CompressedLiquidVolume/Thomson.hpp"

#include "HeatCapacity/AlyLee.hpp"
#include "HeatCapacity/Polynomial.hpp"
#include "HeatCapacity/PPDSLiquid.hpp"

#include "SaturatedLiquidVolume/Elbro.hpp"
#include "SaturatedLiquidVolume/HankinsonThomson.hpp"
#include "SaturatedLiquidVolume/Rackett.hpp"
#include "SaturatedLiquidVolume/YenWoods.hpp"

#include "VaporPressure/AmbroseWalton.hpp"
#include "VaporPressure/AntoineExtended.hpp"
#include "VaporPressure/Riedel.hpp"
#include "VaporPressure/Wagner.hpp"

#include "Viscosity/DIPPR102.hpp"
#include "Viscosity/KirchhoffExtended.hpp"
#include "Viscosity/Lucas.hpp"

#include "CompressedLiquidViscosity/Lucas.hpp"
#include "CompressedVaporViscosity/Lucas.hpp"


#endif    // PCPROPS_PROPERTYLIB_HPP
