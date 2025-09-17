// Execute interface capturing if the model is available
if (pInterfaceCapturing.valid())
{
    // Flag that the alpha subcycle has been executed
    alphaSubCycleExecuted = true;
    
    // Validate fields before operations
    if (alpha1.size() != mesh.nCells() || alpha2.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Alpha field size mismatch: alpha1:" << alpha1.size()
            << " alpha2:" << alpha2.size()
            << " mesh:" << mesh.nCells()
            << abort(FatalError);
    }

    // Check for non-finite values
    const scalarField& alpha1Internal = alpha1.internalField();
    const scalarField& alpha2Internal = alpha2.internalField();
    forAll(alpha1Internal, cellI)
    {
        if (!std::isfinite(alpha1Internal[cellI]) || !std::isfinite(alpha2Internal[cellI]))
        {
            FatalErrorInFunction
                << "Non-finite alpha values detected at cell " << cellI
                << " alpha1:" << alpha1Internal[cellI]
                << " alpha2:" << alpha2Internal[cellI]
                << abort(FatalError);
        }
    }
    
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
            pInterfaceCapturing->correct();
            interfaceCorrectionAppliedInSubCycle = true;

            const scalarField& alphaField = alpha1.internalField();
            forAll(alphaField, i)
            {
                if (!std::isfinite(alphaField[i]))
                {
                    FatalErrorInFunction
                        << "Non-finite value detected in alpha1 at cell " << i
                        << ". Value: " << alphaField[i]
                        << abort(FatalError);
                }
            }

            alpha1 = max(min(alpha1, scalar(1)), scalar(0));
            alpha1.correctBoundaryConditions();
            alpha2 = scalar(1) - alpha1;
            alpha2.correctBoundaryConditions();

            #include "alphaEqn.H"

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

        const scalarField& alphaField = alpha1.internalField();
        forAll(alphaField, i)
        {
            if (!std::isfinite(alphaField[i]))
            {
                FatalErrorInFunction
                    << "Non-finite value detected in alpha1 at cell " << i
                    << ". Value: " << alphaField[i]
                    << abort(FatalError);
            }
        }

        alpha1 = max(min(alpha1, scalar(1)), scalar(0));
        alpha1.correctBoundaryConditions();
        alpha2 = scalar(1) - alpha1;
        alpha2.correctBoundaryConditions();

        #include "alphaEqn.H"
    }

    // Final safety check on alpha1 after sub-cycling
    {
        const scalarField& alphaField = alpha1.internalField();
        forAll(alphaField, i)
        {
            if (!std::isfinite(alphaField[i]))
            {
                FatalErrorInFunction
                    << "Non-finite value detected in alpha1 at cell " << i
                    << ". Value: " << alphaField[i]
                    << abort(FatalError);
            }
        }
    }

    alpha1 = max(min(alpha1, scalar(1)), scalar(0));
    alpha1.correctBoundaryConditions();
    alpha2 = scalar(1) - alpha1;
    alpha2.correctBoundaryConditions();

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

