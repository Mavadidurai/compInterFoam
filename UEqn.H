fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  + transportModel.divDevRhoReff(U)
 ==
    fvOptions(rho, U)
);

{
    volScalarField::Internal contErr
    (
        (
            fvc::ddt(rho)
          + fvc::div(rhoPhi)
        )()
    );

    UEqn -= fvm::Sp(contErr, U);
}

UEqn.relax();
fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
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

    fvOptions.correct(U);

    // Limit unphysical velocities that can destabilise the pressure solve
    const scalar maxRealisticVel = 2000.0;
    vectorField& Uinternal = U.primitiveFieldRef();
    forAll(Uinternal, cellI)
    {
        vector& Ui = Uinternal[cellI];
        const scalar magU = mag(Ui);
        if (std::isfinite(magU))
        {
            if (magU > maxRealisticVel)
            {
                Ui *= maxRealisticVel/magU;
            }
        }
        else
        {
            Ui = vector::zero;
        }
    }

    K = 0.5*magSqr(U);
}


const dictionary& compCoeffs = compInterFoamCoeffsDict(mesh);

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
// Validate the computed flux to ensure all values are finite
{
    const scalar minRAU = 1e-12;

    scalarField& rAUInternal = rAU.primitiveFieldRef();
    bool adjustedRAU = false;
    forAll(rAUInternal, cellI)
    {
        scalar& val = rAUInternal[cellI];
        if (!std::isfinite(val) || val <= minRAU)
        {
            val = minRAU;
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
            if (!std::isfinite(val) || val <= minRAU)
            {
                val = minRAU;
                adjustedRAU = true;
            }
        }
    }

    if (adjustedRAU && verbose && master)
    {
        WarningInFunction << "Clamped non-physical rAU entries" << endl;
    }
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
