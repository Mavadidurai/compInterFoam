#include "DimensionValidator.H"

namespace Foam
{

// Initialise static dimension sets
const dimensionSet DimensionValidator::dimVelocity(0, 1, -1, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimArea(0, 2, 0, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimPressure(1, -1, -2, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimTemperature(0, 0, 0, 1, 0, 0, 0);
const dimensionSet DimensionValidator::dimEnergy(1, 2, -2, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimTempRate(0, 0, -1, 1, 0, 0, 0);
const dimensionSet DimensionValidator::dimDensity(1, -3, 0, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimSpecificHeat(0, 2, -2, -1, 0, 0, 0);
const dimensionSet DimensionValidator::dimThermalCond(1, 1, -3, -1, 0, 0, 0);
const dimensionSet DimensionValidator::dimHeatSource(1, -1, -3, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimLatentHeat(0, 2, -2, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimAcceleration(0, 1, -2, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimMomentum(1, 1, -2, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimSpecificMomentum(0, 1, -2, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimInverseA(-1, 3, 1, 0, 0, 0, 0);

} // End namespace Foam
