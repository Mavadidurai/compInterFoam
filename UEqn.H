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
    // ═══════════════════════════════════════════════════════════
    // CRITICAL FIX: Assemble ALL momentum sources BEFORE solve
    // ═══════════════════════════════════════════════════════════
    
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
    
    // Start with basic forces
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
        tCapillaryForce()
      - fvc::grad(p_rgh)
      + rho*g
    );

    // ✅ ADD RECOIL PRESSURE *BEFORE* SOLVE
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
    
    // ═══════════════════════════════════════════════════════════
    // DIVERGENCE DETECTION: Stop if velocity explodes
    // ═══════════════════════════════════════════════════════════
    const scalar maxPhysicalVelocity = 5000.0;  // 5 km/s (generous for LIFT)
    const scalar maxU = gMax(mag(U)());
    
    if (maxU > maxPhysicalVelocity)
    {
        FatalErrorInFunction
            << "Velocity field has DIVERGED!" << nl
            << "Max velocity: " << maxU << " m/s" << nl
            << "Physical limit: " << maxPhysicalVelocity << " m/s" << nl
            << nl
            << "This indicates numerical instability. Check:" << nl
            << "  1. Time step (current: " << runTime.deltaTValue() << " s)" << nl
            << "  2. Max recoil pressure: " 
            << (recoilPressurePtr ? gMax(*recoilPressurePtr)/1e6 : 0) << " MPa" << nl
            << "  3. Courant number" << nl
            << "  4. Temperature gradients" << nl
            << abort(FatalError);
    }
    
    K = 0.5*magSqr(U);
}

volScalarField rAU("rAU", 1.0/UEqn.A());

phi = fvc::flux(U);
rhoPhi = fvc::interpolate(rho)*phi;
fvc::makeRelative(rhoPhi, U);

U.correctBoundaryConditions();
fvOptions.correct(U);
