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
    UEqn -= fvc::grad(*recoilPressurePtr);
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
