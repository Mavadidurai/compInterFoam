fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  + transportModel.divDevRhoReff(U)
 ==
    fvOptions(rho, U)
);

UEqn.relax();
fvOptions.constrain(UEqn);
const dictionary& compCoeffs = compInterFoamCoeffsDict(mesh);
// CRITICAL FIX: Only validate, don't clamp unless catastrophic
const scalar maxReasonableVel =
    compCoeffs.lookupOrDefault<scalar>("maxReasonableVelocity", 5000.0);

if (pimple.momentumPredictor())
{
    if (recoilPressurePtr)
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(mixture.sigmaK())
                  * fvc::snGrad(alpha1)
                )*mesh.magSf()
            )
          - fvc::grad(p_rgh)
          - fvc::grad(alpha1*(*recoilPressurePtr))
          + rho*g
        );
    }
    else
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(mixture.sigmaK())
                  * fvc::snGrad(alpha1)
                )*mesh.magSf()
            )
          - fvc::grad(p_rgh)
          + rho*g
        );
    }

    fvOptions.correct(U);

    // CRITICAL FIX: Only catch catastrophic non-physical values
    vectorField& Uinternal = U.primitiveFieldRef();
    label nCatastrophic = 0;
    
    forAll(Uinternal, cellI)
    {
        vector& Ui = Uinternal[cellI];
        const scalar magU = mag(Ui);
        
        if (!std::isfinite(magU))
        {
            Ui = vector::zero;
            nCatastrophic++;
        }
        else if (maxReasonableVel > 0 && magU > maxReasonableVel)
        {
            // WARNING only - don't clamp unless truly catastrophic
            if (verbose && master && nCatastrophic == 0)
            {
                WarningInFunction
                    << "High velocity detected: " << magU << " m/s at cell " 
                    << cellI << " (exceeds " << maxReasonableVel << " m/s)" 
                    << endl;
            }
            nCatastrophic++;
        }
    }
    
    if (nCatastrophic > 0 && verbose && master)
    {
        Info<< "Velocity validation: " << nCatastrophic 
            << " cells with extreme values" << endl;
    }

    U.correctBoundaryConditions();
    K = 0.5*magSqr(U);
}

const scalar diagFloorInput =
    compCoeffs.lookupOrDefault<scalar>("minUEqnDiag", 1e-9);
const scalar diagFloor = Foam::max(diagFloorInput, Foam::SMALL);

tmp<volScalarField> tA = UEqn.A();
volScalarField& A = tA.ref();

bool adjustedDiag = false;
auto enforceDiagFloor = [&](scalar& coeff)
{
    if (!std::isfinite(coeff) || coeff < diagFloor)
    {
        coeff = diagFloor;
        adjustedDiag = true;
    }
};

scalarField& AInternal = A.primitiveFieldRef();
forAll(AInternal, cellI)
{
    enforceDiagFloor(AInternal[cellI]);
}

volScalarField::Boundary& Aboundary = A.boundaryFieldRef();
forAll(Aboundary, patchI)
{
    scalarField& patchField = Aboundary[patchI];
    forAll(patchField, faceI)
    {
        enforceDiagFloor(patchField[faceI]);
    }
}

if (adjustedDiag && verbose && master)
{
    WarningInFunction
        << "Detected non-physical UEqn diagonal coefficients; enforcing "
        << "a minimum magnitude of " << diagFloor
        << " before computing rAU" << endl;
}

volScalarField rAU("rAU", 1.0/A);
{
    const Switch enableClamp = compCoeffs.lookupOrDefault<Switch>("enableRAUClamp", false);

    if (enableClamp)
    {
        const scalar minRAU = compCoeffs.lookupOrDefault<scalar>("minRAU", 1e-10);
        const scalar maxRAU = compCoeffs.lookupOrDefault<scalar>("maxRAU", GREAT);

        auto clampValue = [&](scalar& value)
        {
            if (!std::isfinite(value))
            {
                value = minRAU;
            }
            else
            {
                value = Foam::max(minRAU, Foam::min(maxRAU, value));
            }
        };

        scalarField& rAUInternal = rAU.primitiveFieldRef();
        forAll(rAUInternal, cellI)
        {
            clampValue(rAUInternal[cellI]);
        }

        volScalarField::Boundary& rAUb = rAU.boundaryFieldRef();
        forAll(rAUb, patchI)
        {
            scalarField& patchField = rAUb[patchI];
            forAll(patchField, faceI)
            {
                clampValue(patchField[faceI]);
            }
        }
    }
}
phi = fvc::flux(U);
const scalar safeRAU = 1.0/diagFloor;

scalarField& rAUInternal = rAU.primitiveFieldRef();
bool adjustedRAU = false;
forAll(rAUInternal, cellI)
{
    scalar& val = rAUInternal[cellI];
    if (!std::isfinite(val) || val <= 0)
    {
        val = safeRAU;
        adjustedRAU = true;
    }
}

volScalarField::Boundary& rAUb = rAU.boundaryFieldRef();
forAll(rAUb, patchI)
{
    scalarField& patchField = rAUb[patchI];
    forAll(patchField, faceI)
    {
        scalar& val = patchField[faceI];
        if (!std::isfinite(val) || val <= 0)
        {
            val = safeRAU;
            adjustedRAU = true;
        }
    }
}

if (adjustedRAU && verbose && master)
{
    WarningInFunction << "Clamped non-physical rAU entries" << endl;
}

const scalarField& phiInternal = phi.internalField();
forAll(phiInternal, faceI)
{
    if (!std::isfinite(phiInternal[faceI]))
    {
        FatalErrorInFunction
            << "Detected non-finite phi at face " << faceI
            << " during time " << runTime.timeName()
            << exit(FatalError);
    }
}

//const fvMesh& mesh_ = phi.mesh();
const surfaceScalarField::Boundary& phiBf = phi.boundaryField();
forAll(phiBf, patchI)
{
    const auto& patchField = phiBf[patchI];
    const label patchStart = patchField.patch().start();

    forAll(patchField, faceI)
    {
        if (!std::isfinite(patchField[faceI]))
        {
            FatalErrorInFunction
                << "Detected non-finite phi at face "
                << patchStart + faceI
                << " (patch '" << patchField.patch().name() << "')"
                << " during time " << runTime.timeName()
                << exit(FatalError);
        }
    }
}
rhoPhi = fvc::interpolate(rho)*phi;
fvc::makeRelative(rhoPhi, U);

U.correctBoundaryConditions();
fvOptions.correct(U);
