/*---------------------------------------------------------------------------*\
  Temperature equation for LIFT simulation with two-temperature model
\*---------------------------------------------------------------------------*/
extern Foam::Switch verbose;
{
    if (verbose)
    {
        Info<< "Solving temperature equation" << endl;
    }
    // Solve the compressible temperature equation to evolve the gas phase
    tmp<volScalarField> tContErr = fvc::div(rhoPhi);
    const volScalarField& contErr = tContErr();

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(rhoPhi, T)
      - fvm::Sp(contErr, T)
      - fvm::laplacian(transportModel.alphaEff(), T)
      + (
            fvc::div(fvc::absolute(phi, U), p)()()
          + (fvc::ddt(rho, K) + fvc::div(rhoPhi, K))()()
          - (contErr*K).dimensionedInternalField()
        )
       *
        (
            alpha1/mixture.thermo1().Cv()()
          + alpha2/mixture.thermo2().Cv()()
        )
     ==
        fvOptions(rho, T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    tContErr.clear();

    // Store the solved temperature field for gas-dominated regions
    const scalarField solvedTInternal(T.internalField());

    // Since ttm.solve(laserSource, mixture.phaseChangeSource()) is called
    // before this in compInterFoam.C, the two-temperature model has already
    // solved for Te and Tl.
    // Selectively synchronize the mixture temperature with the lattice
    // temperature within metal-rich cells while retaining the gas solution.
    const dictionary& ctrl = runTime.controlDict();
    const scalar alphaMetalThreshold =
        ctrl.lookupOrDefault<scalar>("alphaMetalThreshold", 0.7);

    const volScalarField& TlField = ttm.Tl();

    label metalSyncedCells = 0;
    label gasCells = 0;
    scalar metalMaxAdjustment = 0.0;
    scalar gasTempSum = 0.0;
    scalar gasMinTemp = GREAT;
    scalar gasMaxTemp = -GREAT;
    scalar maxGasTlDiff = 0.0;

    forAll(T, cellI)
    {
        const scalar alphaVal = alpha1[cellI];
        const scalar TlVal = TlField[cellI];
        const scalar solvedT = solvedTInternal[cellI];

        if (alphaVal >= alphaMetalThreshold)
        {
            T[cellI] = TlVal;
            metalMaxAdjustment = Foam::max
            (
                metalMaxAdjustment,
                Foam::mag(solvedT - TlVal)
            );
            ++metalSyncedCells;
        }
        else
        {
            T[cellI] = solvedT;
            gasMinTemp = Foam::min(gasMinTemp, solvedT);
            gasMaxTemp = Foam::max(gasMaxTemp, solvedT);
            gasTempSum += solvedT;
            maxGasTlDiff = Foam::max
            (
                maxGasTlDiff,
                Foam::mag(solvedT - TlVal)
            );
            ++gasCells;
        }
    }

    reduce(metalSyncedCells, sumOp<label>());
    reduce(gasCells, sumOp<label>());
    reduce(metalMaxAdjustment, maxOp<scalar>());
    reduce(gasMinTemp, minOp<scalar>());
    reduce(gasMaxTemp, maxOp<scalar>());
    reduce(gasTempSum, sumOp<scalar>());
    reduce(maxGasTlDiff, maxOp<scalar>());

    scalar meanGasTemp = 0.0;
    if (gasCells > 0)
    {
        meanGasTemp = gasTempSum/Foam::max(static_cast<scalar>(gasCells), VSMALL);
    }

    if (verbose)
    {
        Info<< "Selective temperature synchronisation:" << nl
            << "  alpha1 threshold: " << alphaMetalThreshold << nl
            << "  Cells synced to Tl: " << metalSyncedCells << nl
            << "  Max |Tl - solved T| correction: "
            << metalMaxAdjustment << " K" << endl;

        if (gasCells > 0)
        {
            Info<< "  Gas-dominated cells: " << gasCells << nl
                << "    Gas T range: " << gasMinTemp << " - "
                << gasMaxTemp << " K" << nl
                << "    Gas T mean: " << meanGasTemp << " K" << nl
                << "    Max |Tgas - Tl|: " << maxGasTlDiff << " K" << endl;
        }
        else
        {
            Info<< "  Gas-dominated cells: none detected" << endl;
        }
    }

    // Evaluate global temperature statistics after selective synchronisation
    scalar minT = min(T).value();
    scalar maxT = max(T).value();
    scalar meanT = gAverage(T);

    if (verbose)
    {
        Info<< "Temperature field after selective synchronisation:" << nl
            << "  T range: " << minT << " - " << maxT << " K" << nl
            << "  Mean T: " << meanT << " K" << endl;
    }

    // Apply physical bounds to prevent numerical issues
    const scalar minPhysicalTemp = ctrl.lookupOrDefault<scalar>("minT", 250.0);
    const scalar maxPhysicalTemp = ctrl.lookupOrDefault<scalar>("maxT", 15000.0);
    const Switch enableTClamp = ctrl.lookupOrDefault<Switch>("enableTClamp", true);
    label minTCells = 0;
    label maxTCells = 0;
    if (enableTClamp)
    {
        forAll(T, cellI)
        {
            if (T[cellI] < minPhysicalTemp)
            {
                T[cellI] = minPhysicalTemp;
                minTCells++;
            }
            else if (T[cellI] > maxPhysicalTemp)
            {
                T[cellI] = maxPhysicalTemp;
                maxTCells++;
            }
        }
        reduce(minTCells, sumOp<label>());
        reduce(maxTCells, sumOp<label>());
        if (minTCells > 0 || maxTCells > 0)
        {
            if (verbose)
            {
                Info<< "Temperature clamping stats: " << minTCells
                    << " cells to minT (" << minPhysicalTemp << "), "
                    << maxTCells << " cells to maxT (" << maxPhysicalTemp << ")" << endl;
            }
        }
    }
    else if (minT < minPhysicalTemp || maxT > maxPhysicalTemp)
    {
        if (verbose)
        {
            Info<< "Temperature out of bounds: range " << minT << " - " << maxT
                << " K; allowed range " << minPhysicalTemp << " - "
                << maxPhysicalTemp << " K" << endl;
        }
    }
    // Update statistics after any clamping
    minT = min(T).value();
    maxT = max(T).value();
    meanT = gAverage(T);
    // Ensure boundary conditions are properly applied
    T.correctBoundaryConditions();
    // Update mixture thermodynamic properties based on new temperature
    mixture.correctThermo();
    mixture.correct();
    // Validate thermodynamic consistency
    const volScalarField& rho1 = mixture.thermo1().rho();
    const volScalarField& rho2 = mixture.thermo2().rho();
    // Check for any thermodynamic inconsistencies
    scalar minRho1 = min(rho1).value();
    scalar minRho2 = min(rho2).value();
    if (minRho1 <= 0 || minRho2 <= 0)
    {
        WarningInFunction
            << "Non-positive density detected after temperature update:" << nl
            << "  min(rho1) = " << minRho1 << " kg/m³" << nl
            << "  min(rho2) = " << minRho2 << " kg/m³" << nl
            << "Check thermodynamic properties" << endl;
    }
    // Report temperature statistics specific to LIFT process
    if (maxT > mixture.T_melt())
    {
        // Count cells above melting point
        label meltedCells = 0;
        scalar meltedVolume = 0.0;  
        forAll(T, cellI)
        {
            if (T[cellI] > mixture.T_melt())
            {
                meltedCells++;
                meltedVolume += mesh.V()[cellI] * alpha1[cellI];
            }
        }
        reduce(meltedCells, sumOp<label>());
        reduce(meltedVolume, sumOp<scalar>());
        if (verbose)
        {
            Info<< "LIFT Process Statistics:" << nl
                << "  Cells above melting point: " << meltedCells << nl
                << "  Melted metal volume: " << meltedVolume*1e18 << " μm³" << nl
                << "  Max temperature: " << maxT << " K" << nl
                << "  Melting point: " << mixture.T_melt() << " K" << endl;

            if (maxT > mixture.T_melt() + 1000.0)
            {
                Info<< "Note: Significant overheating detected - typical for laser heating" << endl;
            }
        }
    }
    // Check electron-lattice temperature difference for TTM validation
    volScalarField tempDiff = mag(ttm.Te() - ttm.Tl());
    scalar maxTempDiff = max(tempDiff).value();
    if (maxTempDiff > 500.0)  // Significant electron-lattice temperature difference
    {
        if (verbose)
        {
            Info<< "Two-temperature model status:" << nl
                << "  Max Te-Tl difference: " << maxTempDiff << " K" << nl
                << "  This indicates strong electron-phonon coupling effects" << endl;
        }
    }
    if (verbose)
    {
        Info<< "Temperature equation completed successfully" << endl;
    }
}
