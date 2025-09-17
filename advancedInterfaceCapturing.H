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
    const surfaceScalarField& phi_;  // Flux field
    const twoPhaseMixtureThermo& mixture_;
    const volScalarField& T_;        // Temperature field for phase change
    // Phase-change parameters provided by twoPhaseMixtureThermo
    dimensionedScalar meltingTemp_;             // Melting temperature (Tsol)
    dimensionedScalar vaporTemp_;               // Vaporization temperature (Tvap)
    dimensionedScalar latentHeat_;              // Latent heat of fusion (hf)
        // User-configurable recoil pressure parameters
    scalar pressureScale_;           // Scaling factor for evaporation rate
    scalar recoilMax_;               // Base recoil pressure cap
    label  recoilUpdateInterval_;    // Steps between recoil pressure updates
    dimensionedScalar recoilTempOffset_;        // Temperature offset for recoil activation
    dimensionedScalar phaseChangeTempOffset_;   // Offset for phase-change threshold
    Switch clampRecoil_;             // Enable/disable recoil pressure clamping
    Switch scaleRecoilMax_;          // Scale recoilMax with phase-change rate
    scalar relaxFactor_;             // Relaxation factor for alpha transport
    // Valid alpha range (controlDict → advancedInterfaceCapturing)
    scalar alphaMin_;               // Minimum alpha for recoil activation
    scalar alphaMax_;               // Maximum alpha for recoil activation
   // Ejection velocity coefficients
    scalar C0_;
    scalar C1_;
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
    label recoilUpdateInterval() const { return throttleRecoilUpdates_ ? recoilUpdateInterval_ : 1; }
    Switch throttleRecoilUpdates() const { return throttleRecoilUpdates_; }
    scalar relaxFactor() const { return relaxFactor_; }
    scalar C0() const { return C0_; }
    scalar C1() const { return C1_; }
    // Access functions for validation
    dimensionedScalar latentHeat() const { return latentHeat_; }
    dimensionedScalar meltingTemperature() const { return meltingTemp_; }
    dimensionedScalar vaporTemperature() const { return vaporTemp_; }
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
