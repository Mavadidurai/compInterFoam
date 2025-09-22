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
          - (fvOptions(rho, U) & rho)
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
rhoPhi = fvc::interpolate(rho)*phi;
fvc::makeRelative(rhoPhi, U);

U.correctBoundaryConditions();
fvOptions.correct(U);
