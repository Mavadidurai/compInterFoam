#ifndef DimensionValidator_H
#define DimensionValidator_H

#include "dimensionSet.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrix.H"
#include "GeometricField.H"
#include "fvPatchField.H"
#include "fvsPatchField.H"
#include "volMesh.H"

namespace Foam
{

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

    //- Check volume field dimensions
    template<class Type>
    static void checkVolFieldDimensions
    (
        const GeometricField<Type, fvPatchField, volMesh>& field,
        const dimensionSet& expectedDimensions,
        const word& fieldName
    )
    {
        if (field.dimensions() != expectedDimensions)
        {
            FatalErrorInFunction
                << "Dimension mismatch for volume field \"" << fieldName << "\"" << nl
                << "Expected: " << expectedDimensions << nl
                << "Actual: " << field.dimensions()
                << abort(FatalError);
        }
    }

    //- Check surface field dimensions
    template<class Type>
    static void checkSurfaceFieldDimensions
    (
        const GeometricField<Type, fvsPatchField, surfaceMesh>& field,
        const dimensionSet& expectedDimensions,
        const word& fieldName
    )
    {
        if (field.dimensions() != expectedDimensions)
        {
            FatalErrorInFunction
                << "Dimension mismatch for surface field \"" << fieldName << "\"" << nl
                << "Expected: " << expectedDimensions << nl
                << "Actual: " << field.dimensions()
                << abort(FatalError);
        }
    }

    //- Validate temperature equation terms
    static void validateTEqnTerms
    (
        const volScalarField& T,
        const volScalarField& rho,
        const surfaceScalarField& phi,
        const volScalarField& alphaEff,
        const volScalarField& Q_laser,
        const volScalarField& phaseChangeSource
    );

    //- Validate source terms
    static void validateSourceTerms
    (
        const volScalarField& Q_laser,
        const volScalarField& phaseChangeSource,
        const dimensionedScalar& Cp,
        const dimensionedScalar& rho,
        const dimensionedScalar& latentHeat
    );

    //- Check basic physical constraints
    static void checkPhysicalConstraints
    (
    volScalarField& T,
    volScalarField& p,
    volScalarField& alpha1
    );

    //- Validate material properties
    static void validateMaterialProperties
    (
        const dimensionedScalar& Cp,
        const dimensionedScalar& k,
        const dimensionedScalar& rho,
        const dimensionedScalar& latentHeat
    );

    //- Validate momentum equation terms
    static bool validateMomentumEquation
    (
        const fvVectorMatrix& UEqn,
        const volVectorField& HbyA,
        const volScalarField& rAU,
        const surfaceScalarField& surfaceTensionForce
    );

    //- Validate pressure equation terms
    static bool validatePressureEquation
    (
        const volScalarField& p_rgh,
        const surfaceScalarField& phi,
        const surfaceScalarField& phiHbyA,
        const surfaceScalarField& rAUf
    );
};

} // End namespace Foam

#endif
