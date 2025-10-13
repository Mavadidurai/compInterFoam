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

if (pimple.momentumPredictor())
{
    tmp<volVectorField> tCapillaryForce
    (
        fvc::reconstruct
        (
            (
                fvc::interpolate(mixture.sigmaK())
              * fvc::snGrad(alpha1)
            )*mesh.magSf()
        )
    );
    volVectorField momentumSource
    (
        tCapillaryForce()
      - fvc::grad(p_rgh)
      + rho*g
    );

    if (recoilPressurePtr)
    {
        momentumSource -= fvc::grad(alpha1*(*recoilPressurePtr));
    }

    solve(UEqn == momentumSource);

    fvOptions.correct(U);
    U.correctBoundaryConditions();
    K = 0.5*magSqr(U);
}

volScalarField rAU("rAU", 1.0/UEqn.A());

phi = fvc::flux(U);
rhoPhi = fvc::interpolate(rho)*phi;
fvc::makeRelative(rhoPhi, U);

U.correctBoundaryConditions();
fvOptions.correct(U);
