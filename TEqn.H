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
    tmp<volScalarField> tContErrK(contErr*K);
    const volScalarField& contErrK = tContErrK();
    tmp<volScalarField> tG = ttm.electronPhononCoupling();
    const volScalarField& GField = tG();
    const volScalarField& TeField = ttm.Te();
    const volScalarField& TlField = ttm.Tl();
    volScalarField& gasMetalHeatFlux =
        mesh.lookupObjectRef<volScalarField>("gasMetalHeatFlux");
    tmp<volScalarField> tInterfaceExchange = ttm.gasMetalExchangeCoeffField();
    const volScalarField& interfaceExchange = tInterfaceExchange();
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(rhoPhi, T)
      - fvm::Sp(contErr, T)
      - fvm::laplacian(transportModel.alphaEff(), T)
      + fvm::Sp(alpha1*GField, T)
     ==
        fvOptions(rho, T)
      + alpha1*GField*TeField
    );

    tmp<volScalarField> tDivAbsolutePhiP
    (
        fvc::div(fvc::absolute(phi, U), p)
    );

    tmp<volScalarField> tKTransport
    (
        fvc::ddt(rho, K) + fvc::div(rhoPhi, K)
    );

    tmp<volScalarField> tCv1 = mixture.thermo1().Cv();
    tmp<volScalarField> tCv2 = mixture.thermo2().Cv();

    tmp<volScalarField> tHeatCapacityWeighting
    (
        alpha1/tCv1()
      + alpha2/tCv2()
    );

    tCv1.clear();
    tCv2.clear();

    tmp<volScalarField> tPressureKSource
    (
        tHeatCapacityWeighting()
       *
        (
            tDivAbsolutePhiP()
          + tKTransport()
          - contErrK
        )
    );

    TEqn += tPressureKSource();

    tPressureKSource.clear();
    tHeatCapacityWeighting.clear();
    tKTransport.clear();
    tDivAbsolutePhiP.clear();

    TEqn += fvm::Sp(interfaceExchange, T);
    TEqn.source() += interfaceExchange*TlField;

    TEqn.relax();

    fvOptions.constrain(TEqn);
    const dictionary& ctrl = runTime.controlDict();

    TEqn.solve();

    fvOptions.correct(T);

    scalar fluxRelax = Foam::max(scalar(0), Foam::min(scalar(1), thermalFluxRelax));
    tmp<volScalarField> tRawFlux = interfaceExchange*(T - TlField);
    const volScalarField& rawFlux = tRawFlux();

    if (fluxRelax < 1.0 - SMALL)
    {
        gasMetalHeatFlux.internalField() =
            fluxRelax*rawFlux.internalField()
          + (1.0 - fluxRelax)*gasMetalHeatFlux.internalField();

        fvBoundaryScalarField& gasMetalHeatFluxBf =
            gasMetalHeatFlux.boundaryFieldRef();
        const fvBoundaryScalarField& rawFluxBf = rawFlux.boundaryField();

        forAll(gasMetalHeatFluxBf, patchI)
        {
            gasMetalHeatFluxBf[patchI] =
                fluxRelax*rawFluxBf[patchI]
              + (1.0 - fluxRelax)*gasMetalHeatFluxBf[patchI];
        }
    }
    else
    {
        gasMetalHeatFlux = rawFlux;
    }

    gasMetalHeatFlux.correctBoundaryConditions();

    if (verbose)
    {
        Info<< "  Updated gas-metal heat flux range: ["
            << gMin(gasMetalHeatFlux) << ", "
            << gMax(gasMetalHeatFlux) << "] W/m³" << endl;
    }

    tRawFlux.clear();
    tContErr.clear();
    tContErrK.clear();

    scalar couplingMin = gMin(interfaceExchange);
    scalar couplingMax = gMax(interfaceExchange);

    tmp<volScalarField> tCouplingDelta = mag(T - TlField);
    const volScalarField& couplingDelta = tCouplingDelta();
    scalar meanCouplingDelta = gAverage(couplingDelta);
    scalar maxCouplingDelta = gMax(couplingDelta);
    tCouplingDelta.clear();

    if (verbose)
    {
        Info<< "Interface coupling diagnostics:" << nl
            << "  Exchange coeff range: " << couplingMin << " - "
            << couplingMax << " W/m³/K" << nl
            << "  Heat flux range: " << gMin(gasMetalHeatFlux)
            << " - " << gMax(gasMetalHeatFlux) << " W/m³" << nl
            << "  Mean |T - Tl|: " << meanCouplingDelta << " K" << nl
            << "  Max |T - Tl|: " << maxCouplingDelta << " K" << endl;
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
    volScalarField tempDiff = mag(TeField - TlField);
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
    tInterfaceExchange.clear();    
    if (verbose)
    {
        Info<< "Temperature equation completed successfully" << endl;
    }
}
