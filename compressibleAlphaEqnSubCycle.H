// Check if using standard subcycling or advanced interface capturing
// After: guard against null pointer and use the correct pointer name
if (pInterfaceCapturing.valid())
{
    // perform advanced interface capturing only if the object exists
    pInterfaceCapturing->correct();

    // Store the alphaPhi10 flux for further calculations
    tmp<surfaceScalarField> talphaPhi1(alphaPhi10);
    const surfaceScalarField& alphaPhi1 = talphaPhi1();

    // Update density field based on new alpha1 distribution
    rho = alpha1*rho1 + alpha2*rho2;

    // Create alphaPhi2 field based on conservation
    surfaceScalarField alphaPhi2("alphaPhi2", phi - alphaPhi1);

    // Calculate continuity error for momentum equation
    volScalarField::Internal contErr
    (
        (
            fvc::ddt(rho) + fvc::div(rhoPhi)
          - (fvOptions(alpha1, mixture.thermo1().rho())&rho1)
          - (fvOptions(alpha2, mixture.thermo2().rho())&rho2)
        )()
    );

    // Ensure alpha1 field is bounded between 0 and 1
    alpha1.min(0.0);
    alpha1.max(1.0);

    // Check for phase fraction correctness
    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}

else
{
    // Use standard alpha subcycling
    tmp<surfaceScalarField> talphaPhi1(alphaPhi10);
    
    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
    
        talphaPhi1.reset
        (
            new surfaceScalarField
            (
                mesh.newIOobject("alphaPhi1"),
                mesh,
                dimensionedScalar(alphaPhi10.dimensions(), Zero)
            )
        );
    
        surfaceScalarField rhoPhiSum
        (
            mesh.newIOobject("rhoPhiSum"),
            mesh,
            dimensionedScalar(rhoPhi.dimensions(), Zero)
        );
    
        tmp<volScalarField> trSubDeltaT;
    
        if (LTS)
        {
            trSubDeltaT =
                fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
        }
    
        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            #include "alphaEqn.H"
            talphaPhi1.ref() += (runTime.deltaT()/totalDeltaT)*alphaPhi10;
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        }
    
        rhoPhi = rhoPhiSum;
    }
    else
    {
        #include "alphaEqn.H"
    }
    
    rho = alpha1*rho1 + alpha2*rho2;
    
    const surfaceScalarField& alphaPhi1 = talphaPhi1();
    surfaceScalarField alphaPhi2("alphaPhi2", phi - alphaPhi1);
    
    volScalarField::Internal contErr
    (
        (
            fvc::ddt(rho) + fvc::div(rhoPhi)
          - (fvOptions(alpha1, mixture.thermo1().rho())&rho1)
          - (fvOptions(alpha2, mixture.thermo2().rho())&rho2)
        )()
    );
}
