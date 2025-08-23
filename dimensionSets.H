#ifndef dimensionSets_H
#define dimensionSets_H

#include "dimensionSet.H"

#include <dimensionSets.H>

namespace Foam
{
    static const dimensionSet dimTempRate(0, 0, -1, 1, 0, 0, 0);
    static const dimensionSet dimSpecificHeat(0, 2, -2, -1, 0, 0, 0);
    static const dimensionSet dimThermalCond(1, 1, -3, -1, 0, 0, 0);
    static const dimensionSet dimHeatSource(1, -1, -3, 0, 0, 0, 0);
    static const dimensionSet dimLatentHeat(0, 2, -2, 0, 0, 0, 0);
    static const dimensionSet dimMomentum(1, 1, -2, 0, 0, 0, 0);
    static const dimensionSet dimSpecificMomentum(0, 1, -2, 0, 0, 0, 0);
    static const dimensionSet dimInverseA(-1, 3, 1, 0, 0, 0, 0);
}

#endif
