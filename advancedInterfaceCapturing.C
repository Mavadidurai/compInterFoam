#include "advancedInterfaceCapturing.H"
#include "fvc.H"
#include "fvm.H"
#include "wallFvPatch.H"
#include "fvPatchField.H"
#include "DimensionValidator.H"

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
        .lookupOrDefault<dimensionedScalar>
        (
            "latentHeat",
            dimensionedScalar
            (
                "latentHeat",
                DimensionValidator::dimLatentHeat,
                mesh.lookupObject<dictionary>("thermophysicalProperties")
                .subDict("metal").lookupOrDefault<scalar>("hf", 435e3)
            )
        )
    ),
    pressureScale_(20000.0),
    recoilMax_(5e6),
    recoilUpdateInterval_(5),
    recoilTempOffset_(100.0),
    alphaMin_(0.01),
    alphaMax_(0.99),
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
    ),
    callCount_(0)
{
        const dictionary& aicDict =
        mesh.time().controlDict().subOrEmptyDict("advancedInterfaceCapturing");

    pressureScale_ = aicDict.lookupOrDefault<scalar>
    (
        "pressureScale",
        pressureScale_
    );
    recoilMax_ = aicDict.lookupOrDefault<scalar>("recoilMax", recoilMax_);
    recoilUpdateInterval_ = aicDict.lookupOrDefault<label>
    (
        "recoilUpdateInterval",
        recoilUpdateInterval_
    );
    recoilTempOffset_ = aicDict.lookupOrDefault<scalar>
    (
        "recoilTempOffset",
        recoilTempOffset_
    );

    const dictionary& alphaBounds = aicDict.subOrEmptyDict("alphaBounds");
    alphaMin_ = alphaBounds.lookupOrDefault<scalar>("alphaMin", alphaMin_);
    alphaMax_ = alphaBounds.lookupOrDefault<scalar>("alphaMax", alphaMax_);
    // Simple initialization, no calculations in constructor to avoid MPI issues
    const bool verbose = mesh.time().controlDict().lookupOrDefault<Switch>("verbose", false);
    if (verbose)
    {
        Info<< "Advanced interface capturing initialized" << endl;
    }
}

void advancedInterfaceCapturing::calculateRecoilPressure()
{
    // OPTIMIZED: Calculate once, reuse multiple times
    const scalar currentTime = mesh_.time().value();
    const scalar maxTemp = max(T_).value();
    const scalar minTempThreshold = meltingTemp_ - recoilTempOffset_;

    
    const bool verbose = mesh_.time().controlDict().lookupOrDefault<Switch>("verbose", false);

    if (verbose)
    {
        Info<< "Calculating recoil pressure at t = " << currentTime
            << "s, max T = " << maxTemp
            << "K, melting temp = " << meltingTemp_ << "K" << endl;
    }

    
    // OPTIMIZED: Early exit using cached values
    if (maxTemp < minTempThreshold)
    {
        // OPTIMIZED: Use cached zero value instead of creating new one
        recoilPressure_ = dimensionedScalar("zero", dimPressure, 0.0);
        if (verbose)
        {
            Info<< "Temperature too low for recoil pressure" << endl;
        }
        return;
    }
    
// Initialize recoil pressure field
    recoilPressure_ = dimensionedScalar("zero", dimPressure, 0.0);

    // Constants for recoil pressure model sourced from dictionaries
    const scalar pressureScale = pressureScale_;
    const scalar recoilMax = recoilMax_;

    // Access fields only once for efficiency
    const scalarField& TField = T_.primitiveField();
    const scalarField& alpha1Field = alpha1_.primitiveField();
    const volScalarField& phaseChange = mixture_.phaseChangeSource();
    const scalarField& phaseChangeField = phaseChange.primitiveField();
    scalarField& recoilField = recoilPressure_.primitiveFieldRef();

    // Compute recoil pressure based on evaporation rate
    forAll(TField, cellI)
    {
        if (TField[cellI] < minTempThreshold) continue;

        const scalar alpha = alpha1Field[cellI];
        scalar alphaDamp = 4.0 * alpha * (1.0 - alpha);
        if (alpha < alphaMin_ || alpha > alphaMax_)
        {
            alphaDamp = 0.0;
        }

        const scalar pressureValue = pressureScale * max(phaseChangeField[cellI], 0.0);
        recoilField[cellI] = min(pressureValue * alphaDamp, recoilMax);
    }
    // Ensure boundary conditions are correct
    recoilPressure_.correctBoundaryConditions();
}

void Foam::advancedInterfaceCapturing::correct()
{
    const bool verbose = mesh_.time().controlDict().lookupOrDefault<Switch>("verbose", false);
    if (verbose)
    {
        Info<< "Performing simplified interface capturing" << endl;
    }
    
    // First update the recoil pressure field
    // Only do this every few steps to reduce MPI communication
    if (recoilUpdateInterval_ <= 1 || callCount_++ % recoilUpdateInterval_ == 0)
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

    if (verbose)
    {
        Info<< "Phase-1 volume fraction = "
            << alpha1_.weightedAverage(mesh_.V()).value()
            << "  Min(alpha1) = " << min(alpha1_).value()
            << "  Max(alpha1) = " << max(alpha1_).value()
            << "  Max recoil pressure = " << max(recoilPressure_).value()
            << endl;
    }
}

void advancedInterfaceCapturing::write() const
{
    recoilPressure_.write();
}

} // End namespace Foam
