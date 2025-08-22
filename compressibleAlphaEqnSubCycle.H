// Check if using standard subcycling or advanced interface capturing
// After: guard against null pointer and use the correct pointer name
if (pInterfaceCapturing.valid())
{
    // Flag that the alpha subcycle has been executed
    alphaSubCycleExecuted = true;

    // perform advanced interface capturing only if the object exists
    pInterfaceCapturing->correct();

    // Ensure alpha1 field is bounded between 0 and 1
    alpha1 = max(min(alpha1, scalar(1)), scalar(0));
    alpha1.correctBoundaryConditions();

    // Update the complementary phase fraction
    alpha2 = scalar(1) - alpha1;
    alpha2.correctBoundaryConditions();

    // Recompute density and fluxes based on bounded alpha values
    alphaPhi10 = fvc::interpolate(alpha1)*phi;
    rho = alpha1*rho1 + alpha2*rho2;
    rhoPhi = alphaPhi10*fvc::interpolate(rho1) + (phi - alphaPhi10)*fvc::interpolate(rho2);

          if (verbose)
    {
        // Calculate continuity error for diagnostics
        tmp<volScalarField> tContErr = fvc::ddt(rho) + fvc::div(rhoPhi);
        const volScalarField& contErr = tContErr();

        scalar sumLocalContErr = runTime.deltaT().value()
            * mag(contErr)().weightedAverage(mesh.V()).value();
        scalar globalContErr = runTime.deltaT().value()
            * contErr.weightedAverage(mesh.V()).value();

        if (!std::isfinite(sumLocalContErr)) sumLocalContErr = 0.0;
        if (!std::isfinite(globalContErr)) globalContErr = 0.0;
        Info<< "Continuity errors:" << nl
            << "  Sum local: " << sumLocalContErr << nl
            << "  Global: " << globalContErr << endl;

        // Check for phase fraction correctness
        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.V()).value()
            << "  Min(alpha1) = " << min(alpha1).value()
            << "  Max(alpha1) = " << max(alpha1).value()
            << endl;
    }
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
                IOobject
                (
                    "alphaPhi1",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(alphaPhi10.dimensions(), Zero)
            )
        );
    
        surfaceScalarField rhoPhiSum
        (
            IOobject
            (
                "rhoPhiSum",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
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
        
    if (verbose)
    {
        // Calculate continuity error for diagnostics
        tmp<volScalarField> tContErr = fvc::ddt(rho) + fvc::div(rhoPhi);
        const volScalarField& contErr = tContErr();

        scalar sumLocalContErr = runTime.deltaT().value()
            * mag(contErr)().weightedAverage(mesh.V()).value();
        scalar globalContErr = runTime.deltaT().value()
            * contErr.weightedAverage(mesh.V()).value();

        if (!std::isfinite(sumLocalContErr)) sumLocalContErr = 0.0;
        if (!std::isfinite(globalContErr)) globalContErr = 0.0;

        Info<< "Continuity errors:" << nl
            << "  Sum local: " << sumLocalContErr << nl
            << "  Global: " << globalContErr << endl;
    }
}
