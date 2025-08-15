#include "advancedInterfaceCapturing.H"
#include "fvc.H"
#include "fvm.H" 
#include "wallFvPatch.H"
#include "fvPatchField.H"

namespace Foam
{

advancedInterfaceCapturing::advancedInterfaceCapturing
(
    const fvMesh& mesh,
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const twoPhaseMixtureThermo& mixture,
    const volScalarField& T
)
:
    mesh_(mesh),
    alpha1_(alpha1),
    phi_(phi),
    mixture_(mixture),
    T_(T),
    meltingTemp_
    (
        mesh.lookupObject<dictionary>("transportProperties")
        .lookupOrDefault<scalar>("meltingTemperature", 1941.0)
    ),
    vaporTemp_
    (
        mesh.lookupObject<dictionary>("transportProperties")
        .lookupOrDefault<scalar>("vaporTemperature", 3560.0)
    ),
    latentHeat_
    (
        mesh.lookupObject<dictionary>("transportProperties")
        .lookupOrDefault<scalar>("latentHeat", 
            mesh.lookupObject<dictionary>("thermophysicalProperties")
            .subDict("metal").lookupOrDefault<scalar>("hf", 435e3))
    ),
    recoilPressure_
    (
        IOobject
        (
            "recoilPressure",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0)
    )
{
    // Simple initialization, no calculations in constructor to avoid MPI issues
    Info<< "Advanced interface capturing initialized" << endl;
}

void advancedInterfaceCapturing::calculateRecoilPressure()
{
    // OPTIMIZED: Calculate once, reuse multiple times
    const scalar currentTime = mesh_.time().value();
    const scalar maxTemp = max(T_).value();
    const scalar minTempThreshold = meltingTemp_ - 100.0;
    
    Info<< "Calculating recoil pressure at t = " << currentTime 
        << "s, max T = " << maxTemp 
        << "K, melting temp = " << meltingTemp_ << "K" << endl;
    
    // OPTIMIZED: Early exit using cached values
    if (maxTemp < minTempThreshold)
    {
        // OPTIMIZED: Use cached zero value instead of creating new one
        recoilPressure_ = dimensionedScalar("zero", dimPressure, 0.0);
        Info<< "Temperature too low for recoil pressure" << endl;
        return;
    }
    
// Initialize recoil pressure field
    recoilPressure_ = dimensionedScalar("zero", dimPressure, 0.0);

    // Constants for recoil pressure model
    const scalar pressureScale = 20000.0;
    const scalar recoilMax = 5e6;

    // Access fields only once for efficiency
    const scalarField& TField = T_.primitiveField();
    const scalarField& alpha1Field = alpha1_.primitiveField();
    const volScalarField& evap = mixture_.phaseChangeSource();
    const scalarField& evapField = evap.primitiveField();
    scalarField& recoilField = recoilPressure_.primitiveFieldRef();

    // Compute recoil pressure based on evaporation rate
    forAll(TField, cellI)
    {
        if (TField[cellI] < meltingTemp_ - 100.0) continue;

        const scalar alpha = alpha1Field[cellI];
        scalar alphaDamp = 4.0 * alpha * (1.0 - alpha);
        if (alpha < 0.01 || alpha > 0.99)
        {
            alphaDamp = 0.0;
        }

        const scalar pressureValue = pressureScale * max(evapField[cellI], 0.0);
        recoilField[cellI] = min(pressureValue * alphaDamp, recoilMax);
    }
    // Ensure boundary conditions are correct
    recoilPressure_.correctBoundaryConditions();
}

void Foam::advancedInterfaceCapturing::correct()
{
    Info<< "Performing simplified interface capturing" << endl;
    
    // First update the recoil pressure field
    // Only do this every few steps to reduce MPI communication
    static int callCount = 0;
    if (callCount++ % 5 == 0) // Only update every 5 steps
    {
        calculateRecoilPressure();
    }
    
    // Create a simple robust transport equation for alpha1
    fvScalarMatrix alpha1Eqn
    (
        fvm::ddt(alpha1_)
      + fvm::div(phi_, alpha1_)
    );
    
    // Apply relaxation for stability
    scalar relaxFactor = 0.5;
    alpha1Eqn.relax(relaxFactor);
    
    // Solve with standard settings
    alpha1Eqn.solve();
    
    // Simple boundedness check without per-cell limiting
    // This approach reduces MPI communication
    alpha1_ = max(min(alpha1_, scalar(1)), scalar(0));
    alpha1_.correctBoundaryConditions();

    Info<< "Phase-1 volume fraction = "
        << alpha1_.weightedAverage(mesh_.V()).value()
        << "  Min(alpha1) = " << min(alpha1_).value()
        << "  Max(alpha1) = " << max(alpha1_).value()
        << "  Max recoil pressure = " << max(recoilPressure_).value()
        << endl;
}

void advancedInterfaceCapturing::write() const
{
    recoilPressure_.write();
}

} // End namespace Foam
