#ifndef advancedInterfaceCapturing_H
#define advancedInterfaceCapturing_H
#include "fvMesh.H"
#include "volFields.H"
#include "fvc.H"
#include "surfaceFields.H"
#include "twoPhaseMixtureThermo.H"
namespace Foam
{
class advancedInterfaceCapturing
{
private:
    // References
    const fvMesh& mesh_;
    volScalarField& alpha1_;         // Phase fraction
    const twoPhaseMixtureThermo& mixture_;
    const volScalarField& T_;        // Temperature field for phase change
    // Phase-change parameters provided by twoPhaseMixtureThermo
    dimensionedScalar meltingTemp_;             // Melting temperature (Tsol)
    dimensionedScalar vaporTemp_;               // Vaporization temperature (Tvap)
        // User-configurable recoil pressure parameters
    scalar pressureScale_;           // Scaling factor for evaporation rate
    scalar recoilMax_;               // Base recoil pressure cap
    dimensionedScalar recoilTempOffset_;        // Temperature offset for recoil activation
    dimensionedScalar phaseChangeTempOffset_;   // Offset for phase-change threshold
    Switch clampRecoil_;             // Enable/disable recoil pressure clamping
    Switch scaleRecoilMax_;          // Scale recoilMax with interface weighting
    // Valid alpha range (controlDict → advancedInterfaceCapturing)
    scalar alphaMin_;               // Minimum alpha for recoil activation
    scalar alphaMax_;               // Maximum alpha for recoil activation
    // Field for recoil pressure
    volScalarField recoilPressure_;  // Recoil pressure from vaporization

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
    dimensionedScalar meltingTemperature() const { return meltingTemp_; }
    dimensionedScalar vaporTemperature() const { return vaporTemp_; }
    // Write functions
    void write() const;
    void clearInternalFields()
{
    // Reset internal fields to prevent memory issues
    recoilPressure_ = dimensionedScalar("zero", dimPressure, 0.0);
}
};
} // End namespace Foam
#endif
