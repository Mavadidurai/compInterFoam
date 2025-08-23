#ifndef DimensionValidator_H
#define DimensionValidator_H
#include "dimensionSet.H"

namespace Foam
{

//- Holder for common dimension sets used across the solver.
class DimensionValidator
{
public:
    // Standard dimension sets
    static const dimensionSet dimVelocity;          // [m/s]
    static const dimensionSet dimArea;              // [m²]
    static const dimensionSet dimPressure;          // [kg⋅m⁻¹⋅s⁻²]
    static const dimensionSet dimTemperature;       // [K]
    static const dimensionSet dimEnergy;            // [kg⋅m²⋅s⁻²]
    static const dimensionSet dimTempRate;          // [K/s]
    static const dimensionSet dimDensity;           // [kg/m³]
    static const dimensionSet dimSpecificHeat;      // [J/kg·K]
    static const dimensionSet dimThermalCond;       // [W/m·K]
    static const dimensionSet dimHeatSource;        // [W/m³]
    static const dimensionSet dimLatentHeat;        // [J/kg]
    static const dimensionSet dimAcceleration;      // [m/s²]
    static const dimensionSet dimMomentum;          // [kg⋅m/s²]
    static const dimensionSet dimSpecificMomentum;  // [m/s²]
    static const dimensionSet dimInverseA;          // [1/A] (inverse of UEqn.A())
};

} // End namespace Foam

#endif
