// Execute interface capturing if the model is available
auto checkAlphaFinite = [&](const char* context)
{
    const scalarField& alpha1Internal = alpha1.internalField();
    const scalarField& alpha2Internal = alpha2.internalField();

    forAll(alpha1Internal, cellI)
    {
        if
        (
            !std::isfinite(alpha1Internal[cellI])
         || !std::isfinite(alpha2Internal[cellI])
        )
        {
            FatalErrorInFunction
                << "Non-finite alpha values detected (" << context << ")"
                << " at cell " << cellI
                << " alpha1:" << alpha1Internal[cellI]
                << " alpha2:" << alpha2Internal[cellI]
                << abort(FatalError);
        }
    }
};

auto enforceAlphaBounds = [&]()
{
    alpha1 = max(min(alpha1, scalar(1)), scalar(0));
    alpha1.correctBoundaryConditions();
    alpha2 = scalar(1) - alpha1;
    alpha2.correctBoundaryConditions();
};

auto stabilizeAlphaFields = [&](const char* context)
{
    checkAlphaFinite(context);
    enforceAlphaBounds();
    checkAlphaFinite(context);
};
if (pInterfaceCapturing.valid())
{
    // Flag that the alpha subcycle has been executed
    //alphaSubCycleExecuted = true;
    
    // Validate fields before operations
    if (alpha1.size() != mesh.nCells() || alpha2.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Alpha field size mismatch: alpha1:" << alpha1.size()
            << " alpha2:" << alpha2.size()
            << " mesh:" << mesh.nCells()
            << abort(FatalError);
    }

    stabilizeAlphaFields("initial alpha fields before interface capturing");
    
    const bool useSubCycling = nAlphaSubCycles > 1;

    if (useSubCycling)
    {
        tmp<surfaceScalarField> talphaPhi1(alphaPhi10);
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
        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            pInterfaceCapturing->correct();
            interfaceCorrectionAppliedInSubCycle = true;
            stabilizeAlphaFields("after interface correction (sub-cycle)");

            #include "alphaEqn.H"

            stabilizeAlphaFields("post alphaEqn sub-cycle (interface capturing)");

            dimensionedScalar subCycleWeight = runTime.deltaT()/totalDeltaT;
            talphaPhi1.ref() += subCycleWeight*alphaPhi10;
            rhoPhiSum += subCycleWeight*rhoPhi;
        }

        alphaPhi10 = talphaPhi1();
        rhoPhi = rhoPhiSum;
    }
    else
    {
        pInterfaceCapturing->correct();
        interfaceCorrectionAppliedInSubCycle = true;
        stabilizeAlphaFields("after interface correction");

        #include "alphaEqn.H"

        stabilizeAlphaFields("post alphaEqn interface capturing");
    }

    stabilizeAlphaFields("post interface-capturing alpha solve");

    // Recompute density, fluxes and provide diagnostics
    rho = alpha1*rho1 + alpha2*rho2;

    if (!useSubCycling)
    {
        // Store intermediate fields to avoid redundant temporaries
        tmp<surfaceScalarField> trho1 = fvc::interpolate(rho1);
        tmp<surfaceScalarField> trho2 = fvc::interpolate(rho2);
        surfaceScalarField phiDiff = phi - alphaPhi10;
        rhoPhi = alphaPhi10*trho1() + phiDiff*trho2();
    }

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
        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            stabilizeAlphaFields("alpha subcycle (no interface capturing)");
            #include "alphaEqn.H"
            stabilizeAlphaFields("post alpha subcycle (no interface capturing)");
            talphaPhi1.ref() += (runTime.deltaT()/totalDeltaT)*alphaPhi10;
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        }
        rhoPhi = rhoPhiSum;
    }
    else
    {
        stabilizeAlphaFields("alpha solve (no interface capturing)");
        #include "alphaEqn.H"
        stabilizeAlphaFields("post alpha solve (no interface capturing)");
    }
    stabilizeAlphaFields("final alpha state (no interface capturing)");
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

