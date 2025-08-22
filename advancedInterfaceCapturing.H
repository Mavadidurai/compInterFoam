#ifndef advancedInterfaceCapturing_H
#define advancedInterfaceCapturing_H

#include "fvMesh.H"
#include "volFields.H"
#include "fvc.H"
#include "surfaceFields.H"
#include "twoPhaseMixtureThermo.H"
#include "DimensionValidator.H"

namespace Foam
{

class advancedInterfaceCapturing
{
private:
    // References
    const fvMesh& mesh_;
    volScalarField& alpha1_;         // Phase fraction
    const surfaceScalarField& phi_;  // Flux field
    const twoPhaseMixtureThermo& mixture_;
    const volScalarField& T_;        // Temperature field for phase change

    // Physical parameters - minimal set
    scalar meltingTemp_;             // Melting temperature
    scalar vaporTemp_;               // Vaporization temperature
    dimensionedScalar latentHeat_;   // Latent heat of fusion
        // User-configurable recoil pressure parameters
    scalar pressureScale_;           // Scaling factor for evaporation rate
    scalar recoilMax_;               // Maximum allowed recoil pressure
    label  recoilUpdateInterval_;    // Steps between recoil pressure updates
    scalar recoilTempOffset_;         // Temperature offset for recoil activation
    // Valid alpha range (controlDict → advancedInterfaceCapturing → alphaBounds)
    scalar alphaMin_;               // Minimum alpha for recoil activation
    scalar alphaMax_;               // Maximum alpha for recoil activation
    // Field for recoil pressure
    volScalarField recoilPressure_;  // Recoil pressure from vaporization
    
    //- Track number of correct() calls
    label callCount_;


public:
    // Constructor
    advancedInterfaceCapturing
    (
        const fvMesh& mesh,
        volScalarField& alpha1,
        const surfaceScalarField& phi,
        const twoPhaseMixtureThermo& mixture,
        const volScalarField& T
    );

    //- Delete default copy constructor and assignment operator
    advancedInterfaceCapturing(const advancedInterfaceCapturing&) = delete;
    void operator=(const advancedInterfaceCapturing&) = delete;

    // Member functions
    void correct();  // Main solver routine
    void calculateRecoilPressure();  // Calculate recoil pressure

    // Access functions
    const volScalarField& recoilPressure() const { return recoilPressure_; }
    
    // Access functions for validation
    dimensionedScalar latentHeat() const { return latentHeat_; }

    scalar meltingTemperature() const { return meltingTemp_; }
    scalar vaporTemperature() const { return vaporTemp_; }

    // Write functions
    void write() const;
    void clearInternalFields()
{
    // Reset internal fields to prevent memory issues
    recoilPressure_ = dimensionedScalar("zero", dimPressure, 0.0);
    callCount_ = 0;
}
};

} // End namespace Foam

#endif
