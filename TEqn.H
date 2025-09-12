/*---------------------------------------------------------------------------*\
  Temperature equation for LIFT simulation with two-temperature model
\*---------------------------------------------------------------------------*/
extern Foam::Switch verbose;
{
    if (verbose)
    {
        Info<< "Solving temperature equation" << endl;
    }
   // Since ttm.solve(laserSource, mixture.phaseChangeSource()) is called
    // before this in compInterFoam.C, the two-temperature model has already
    // solved for Te and Tl.
    // We now synchronize the mixture temperature with the lattice temperature
    // and ensure thermodynamic consistency.
    // Synchronize mixture temperature with lattice temperature from TTM
    T = ttm.Tl();
    // Validate temperature field for physical bounds
    scalar minT = min(T).value();
    scalar maxT = max(T).value();
    scalar meanT = gAverage(T);
    if (verbose)
    {
        Info<< "Temperature field after TTM synchronization:" << nl
            << "  T range: " << minT << " - " << maxT << " K" << nl
            << "  Mean T: " << meanT << " K" << endl;
    }
    // Apply physical bounds to prevent numerical issues
    const dictionary& ctrl = runTime.controlDict();
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
    // Ensure boundary conditions are properly applied
    T.correctBoundaryConditions();
    // Update mixture thermodynamic properties based on new temperature
    mixture.correctThermo();
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
