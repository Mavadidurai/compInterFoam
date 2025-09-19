#include "advancedInterfaceCapturing.H"
#include "fvc.H"
#include "fvm.H"
#include "wallFvPatch.H"
#include "fvPatchField.H"
#include <dimensionSets.H>

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
        dimensionedScalar
        (
            "meltingTemperature",
            dimTemperature,
            mixture.T_melt()
        )
    ),
    vaporTemp_
    (
        dimensionedScalar
        (
            "vaporTemperature",
            dimTemperature,
            mixture.T_vapor()
        )
    ),
    latentHeat_(mixture.latentHeat()),
    pressureScale_(20000.0),
    recoilMax_(5e6),
    recoilUpdateInterval_(1),
    throttleRecoilUpdates_(false),
    recoilTempOffset_
    (
        dimensionedScalar("recoilTempOffset", dimTemperature, 0.0)
    ),
    phaseChangeTempOffset_
    (
        dimensionedScalar("phaseChangeTempOffset", dimTemperature, 0.0)
    ),
    clampRecoil_(true),
    scaleRecoilMax_(false),
    relaxFactor_(0.5),
    // Relaxed default bounds to reduce clipping of interface values
    alphaMin_(0.001),
    alphaMax_(0.999),
    C0_(0.0),
    C1_(0.0),
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
    const dictionary& aicDict =
        mesh.time().controlDict().subOrEmptyDict("advancedInterfaceCapturing");

        const bool verbose =
        mesh.time().controlDict().lookupOrDefault<Switch>("verbose", false);

    meltingTemp_ = aicDict.lookupOrDefault<dimensionedScalar>
    (
        "meltingTemperature",
        meltingTemp_
    );
    vaporTemp_ = aicDict.lookupOrDefault<dimensionedScalar>
    (
        "vaporTemperature",
        vaporTemp_
    );

    if (verbose)
    {
        Info<< "advancedInterfaceCapturing temperature bounds: melting = "
            << meltingTemp_.value() << " K, vapor = " << vaporTemp_.value()
            << " K" << endl;
    }

    if (meltingTemp_.value() >= vaporTemp_.value())
    {
        FatalErrorInFunction
            << "meltingTemperature (" << meltingTemp_.value()
            << " K) must be less than vaporTemperature ("
            << vaporTemp_.value() << " K). Values read from the"
            << " 'advancedInterfaceCapturing' dictionary."
            << abort(FatalError);
    }
    phaseChangeTempOffset_ = aicDict.lookupOrDefault<dimensionedScalar>
    (
        "phaseChangeTempOffset",
        phaseChangeTempOffset_
    );

    pressureScale_ = aicDict.lookupOrDefault<scalar>
    (
        "pressureScale",
        pressureScale_
    );
    recoilMax_ = aicDict.lookupOrDefault<scalar>("recoilMax", recoilMax_);

    recoilTempOffset_ = aicDict.lookupOrDefault<dimensionedScalar>
    (
        "recoilTempOffset",
        recoilTempOffset_
    );
    if (recoilTempOffset_.value() < 0)
    {
        FatalErrorInFunction
            << "recoilTempOffset (" << recoilTempOffset_.value()
            << ") must be non-negative"
            << abort(FatalError);
    }

    if (recoilTempOffset_.value() >= vaporTemp_.value())
    {
        FatalErrorInFunction
            << "recoilTempOffset (" << recoilTempOffset_.value()
            << ") must be less than vaporTemperature ("
            << vaporTemp_.value() << ")"
            << abort(FatalError);
    }
    clampRecoil_ = aicDict.lookupOrDefault<Switch>
    (
        "clampRecoil",
        clampRecoil_
    );
    scaleRecoilMax_ = aicDict.lookupOrDefault<Switch>
    (
        "scaleRecoilMax",
        scaleRecoilMax_
    );
    relaxFactor_ = aicDict.lookupOrDefault<scalar>
    (
        "relaxFactor",
        relaxFactor_
    );
    alphaMin_ = aicDict.lookupOrDefault<scalar>("alphaMin", alphaMin_);
    alphaMax_ = aicDict.lookupOrDefault<scalar>("alphaMax", alphaMax_);
    C0_ = aicDict.lookupOrDefault<scalar>("C0", C0_);
    C1_ = aicDict.lookupOrDefault<scalar>("C1", C1_);
    // Simple initialization, no calculations in constructor to avoid MPI issues
    if (verbose)
    {
        Info<< "Advanced interface capturing initialized" << endl;
    }
}
void advancedInterfaceCapturing::calculateRecoilPressure()
{
    // Check field validity using size and dimensions
    if (T_.size() != mesh_.nCells() || alpha1_.size() != mesh_.nCells())
    {
        FatalErrorIn("advancedInterfaceCapturing::calculateRecoilPressure()")
            << "Field size mismatch detected: T:" << T_.size() 
            << " alpha1:" << alpha1_.size()
            << " mesh:" << mesh_.nCells()
            << abort(FatalError);
    }

    // OPTIMIZED: Calculate once, reuse multiple times
    const scalar currentTime = mesh_.time().value();
    // Filter the maximum temperature to metal cells using the alpha threshold
    tmp<volScalarField> maskedT = T_*pos0(alpha1_ - alphaMin_);
    const scalar maxTemp = gMax(maskedT());
    const dimensionedScalar minTempThreshold = vaporTemp_ - recoilTempOffset_;
    const bool verbose = mesh_.time().controlDict().lookupOrDefault<Switch>("verbose", false);
    if (verbose)
    {
        Info<< "Calculating recoil pressure at t = " << currentTime
            << "s, max T = " << maxTemp
            << "K, vapor temp = " << vaporTemp_.value() << "K" << endl;
    }
    // OPTIMIZED: Early exit using cached values
    if (maxTemp < minTempThreshold.value())
    {
        // OPTIMIZED: Use cached zero value instead of creating new one
        recoilPressure_ = dimensionedScalar("zero", dimPressure, 0.0);
        if (verbose)
        {
            Info<< "Temperature too low for recoil pressure" << endl;
        }
        recoilPressure_.correctBoundaryConditions();
        return;
    }
// Initialize recoil pressure field
    recoilPressure_ = dimensionedScalar("zero", dimPressure, 0.0);
    // Constants for recoil pressure model sourced from dictionaries
    const scalar pressureScale = pressureScale_;
    const bool clampRecoil = clampRecoil_;
    const bool scaleRecoilMax = scaleRecoilMax_;
    // Access fields only once for efficiency and validate sizes
    if (T_.size() != alpha1_.size() || T_.size() != mesh_.nCells())
    {
        FatalError << "Field size mismatch in calculateRecoilPressure()" << abort(FatalError);
    }
    
    const scalarField& TField = T_.primitiveField();
    const scalarField& alpha1Field = alpha1_.primitiveField();
    
    // Check for non-finite values in fields
    forAll(TField, cellI)
    {
        if (!std::isfinite(TField[cellI]))
        {
            FatalErrorIn("advancedInterfaceCapturing::calculateRecoilPressure()")
                << "Non-finite temperature value detected at cell " << cellI
                << ". Value: " << TField[cellI]
                << abort(FatalError);
        }
        if (!std::isfinite(alpha1Field[cellI]))
        {
            FatalErrorIn("advancedInterfaceCapturing::calculateRecoilPressure()")
                << "Non-finite alpha1 value detected at cell " << cellI
                << ". Value: " << alpha1Field[cellI]
                << abort(FatalError);
        }
    }

    const volScalarField& massRate = mixture_.dgdt();        // [1/s]
    
    if (massRate.size() != mesh_.nCells())
    {
        FatalErrorIn("advancedInterfaceCapturing::calculateRecoilPressure()")
            << "Invalid mass transfer rate field size: " << massRate.size()
            << " expected: " << mesh_.nCells()
            << abort(FatalError);
    }
    
    const scalarField& massRateField = massRate.primitiveField();
    scalarField& recoilField = recoilPressure_.primitiveFieldRef();

    // Utility: smooth step for alpha ramping
    auto smoothStep = [](const scalar x, const scalar x1, const scalar x2)
    {
        if (x <= x1) return scalar(0);
        if (x >= x2) return scalar(1);
        const scalar y = (x - x1)/(x2 - x1);
        return y*y*(3.0 - 2.0*y);
    };

    // Compute recoil pressure based on evaporation rate
    forAll(TField, cellI)
    {
        if (cellI >= TField.size() || cellI >= alpha1Field.size() || cellI >= massRateField.size())
        {
            FatalError << "Cell index out of bounds in calculateRecoilPressure()" << abort(FatalError);
        }
        
        if (TField[cellI] < minTempThreshold.value()) continue;

        const scalar alpha = alpha1Field[cellI];
        scalar ramp = 1.0;
        if (alpha < alphaMin_)
        {
            ramp = smoothStep(alpha, 0.0, alphaMin_);
        }
        else if (alpha > alphaMax_)
        {
            ramp = 1.0 - smoothStep(alpha, alphaMax_, 1.0);
        }
        const scalar alphaDamp = 4.0 * alpha * (1.0 - alpha) * ramp;

        scalar phaseChangeVal =
            TField[cellI] >= (vaporTemp_.value() + phaseChangeTempOffset_.value())
            ? mag(massRateField[cellI])   // vapor rate drives recoil
            : 0.0;
        const scalar pressureValue = pressureScale * phaseChangeVal;
        const scalar localRecoilMax = scaleRecoilMax
            ? recoilMax_ * alphaDamp
            : recoilMax_;
        const scalar unclamped = pressureValue * alphaDamp;
        recoilField[cellI] = clampRecoil ? min(unclamped, localRecoilMax) : unclamped;
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
    // Update the recoil pressure field every invocation so it follows fast transients.
    calculateRecoilPressure();

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

