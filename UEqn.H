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

if (recoilPressurePtr)
{
    tmp<volVectorField> tRecoilGrad = fvc::grad(*recoilPressurePtr);
    volVectorField& recoilGrad = tRecoilGrad.ref();

    const scalar maxAllowedGrad = 1e7; // Pa/m
    forAll(recoilGrad, cellI)
    {
        const scalar gradMag = mag(recoilGrad[cellI]);
        if (gradMag > maxAllowedGrad)
        {
            recoilGrad[cellI] *= maxAllowedGrad/gradMag;
        }
    }

    UEqn -= recoilGrad;
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
    forAll(U, cellI)
    {
        const scalar magU = mag(U[cellI]);
        if (magU > maxRealisticVel)
        {
            U[cellI] *= maxRealisticVel/magU;
        }
    }

    K = 0.5*magSqr(U);
}


volScalarField rAU("rAU", 1.0/UEqn.A());

phi = fvc::flux(U);
// Validate the computed flux to ensure all values are finite
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
