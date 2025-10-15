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
    // The capillary force already appears in phig within pEqn.H via
    // mixture.surfaceTensionForce().  Including it here as well was
    // double-counting the surface-tension contribution, which the run
    // logs revealed as steadily increasing velocities even when the
    // recoil pressure field remained identically zero.  Only assemble
    // the forces that are not already handled through the pressure
    // correction stage.

    // Start with the forces that are *not* treated in phig
    volVectorField momentumSource
    (
        IOobject
        (
            "momentumSource",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        - fvc::grad(p_rgh)
      + rho*g
    );

    if (recoilPressurePtr)
    {
        const volScalarField& recoilP = *recoilPressurePtr;
        
        // Validate recoil pressure before using
        const scalar maxRecoil = gMax(recoilP);
        const scalar minRecoil = gMin(recoilP);
        
        if (maxRecoil > 1e8 || minRecoil < -1e8)
        {
            WarningInFunction
                << "Recoil pressure out of physical range: ["
                << minRecoil << ", " << maxRecoil << "] Pa"
                << nl << "Clamping to [-100 MPa, 100 MPa]"
                << endl;
                
            volScalarField recoilClamped
            (
                Foam::min
                (
                    Foam::max(recoilP, dimensionedScalar("minR", dimPressure, -1e8)),
                    dimensionedScalar("maxR", dimPressure, 1e8)
                )
            );
            
            momentumSource -= fvc::grad(alpha1*recoilClamped);
        }
        else
        {
            momentumSource -= fvc::grad(alpha1*recoilP);
        }
        
        if (verbose && master && runTime.timeIndex() % 10 == 0)
        {
            const dimensionedScalar recoilForce = 
                fvc::domainIntegrate(mag(fvc::grad(alpha1*recoilP)));
                
            Info<< "    Recoil force magnitude: " << recoilForce.value() 
                << " N" << nl
                << "    Max recoil pressure: " << maxRecoil/1e6 
                << " MPa" << endl;
        }
    }

    // Solve with complete momentum source
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
