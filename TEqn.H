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
     ==
        fvOptions(rho, T)
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
    const volScalarField& Cv1 = tCv1();
    const volScalarField& Cv2 = tCv2();

    tmp<volScalarField> tRho1 = mixture.thermo1().rho();
    tmp<volScalarField> tRho2 = mixture.thermo2().rho();
    const volScalarField& rho1 = tRho1();
    const volScalarField& rho2 = tRho2();

    const dimensionedScalar cvFloor
    (
        "cvFloor",
        dimSpecificHeat,
        1e-12
    );

    const dimensionedScalar rhoFloor
    (
        "rhoFloor",
        dimDensity,
        1e-12
    );

    tmp<volScalarField> tCvMetal(max(Cv1, cvFloor));
    tmp<volScalarField> tCvGas(max(Cv2, cvFloor));
    const volScalarField& CvMetal = tCvMetal();
    const volScalarField& CvGas = tCvGas();

    tmp<volScalarField> tYMetal
    (
        alpha1*max(rho1, rhoFloor)/max(rho, rhoFloor)
    );
    tmp<volScalarField> tYGas
    (
        alpha2*max(rho2, rhoFloor)/max(rho, rhoFloor)
    );
    const volScalarField& YMetal = tYMetal();
    const volScalarField& YGas = tYGas();

    tmp<volScalarField> tCvMix
    (
        max
        (
            YMetal*CvMetal + YGas*CvGas,
            cvFloor
        )
    );
    const volScalarField& CvMix = tCvMix();
    tmp<volScalarField> tRhoEffective(max(rho, rhoFloor));
    const volScalarField& rhoEffective = tRhoEffective();

    tmp<volScalarField> tCvol(rhoEffective*CvMix);
    const volScalarField& Cvol = tCvol();

    tmp<volScalarField> tHeatCapacityWeighting
    (
        dimensionedScalar("one", dimless, 1.0)/Cvol
    );
    const volScalarField& heatCapacityWeighting = tHeatCapacityWeighting();

    dimensionedScalar oldGasEnergy("oldGasEnergy", dimEnergy, 0.0);

    if (verbose)
    {
        tmp<volScalarField> tHeGasOld = mixture.thermo2().he();
        const volScalarField& heGasOld = tHeGasOld();
        oldGasEnergy = fvc::domainIntegrate
        (
            alpha2*rho2*(heGasOld - p/rho2)
        );
        tHeGasOld.clear();
    }

    tYMetal.clear();
    tYGas.clear();
    tCvMix.clear();
    tCvMetal.clear();
    tCvGas.clear();
    tCv1.clear();
    tCv2.clear();
    tRho1.clear();
    tRho2.clear();
    tCvol.clear();
    tmp<volScalarField> tPressureKSource
    (
        rhoEffective
       *heatCapacityWeighting
       *(
            tDivAbsolutePhiP()
          + tKTransport()
          - contErrK
        )
    );

    TEqn += tPressureKSource();

    tPressureKSource.clear();
    tRhoEffective.clear();

    tKTransport.clear();
    tDivAbsolutePhiP.clear();

    TEqn += fvm::Sp(interfaceExchange*heatCapacityWeighting, T);
    TEqn.source() +=
        interfaceExchange*heatCapacityWeighting*TlField;

    tHeatCapacityWeighting.clear();

    TEqn.relax();

    fvOptions.constrain(TEqn);
    const dictionary& ctrl = runTime.controlDict();

    const dictionary solvers(mesh.solutionDict().subOrEmptyDict("solvers"));
    const word finalName(T.name() + "Final");

    if (pimple.finalIter() && !solvers.found(finalName))
    {
        dictionary fallback(mesh.solverDict(T.name()));
        TEqn.solve(fallback);
    }
    else
    {
        TEqn.solve();
    }

    fvOptions.correct(T);

    tContErr.clear();
    tContErrK.clear();

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

    // Recompute the gas-metal heat flux using the clamped temperature field
    scalar fluxRelax = Foam::max(scalar(0), Foam::min(scalar(1), thermalFluxRelax));
    tmp<volScalarField> tRawFlux = interfaceExchange*(T - TlField);
    const volScalarField& rawFlux = tRawFlux();

    if (fluxRelax < 1.0 - SMALL)
    {
        scalarField& gasMetalHeatFluxPf =
            gasMetalHeatFlux.primitiveFieldRef();
        const scalarField& rawFluxPf = rawFlux.primitiveField();

        gasMetalHeatFluxPf =
            fluxRelax*rawFluxPf
          + (1.0 - fluxRelax)*gasMetalHeatFluxPf;

        volScalarField::Boundary& gasMetalHeatFluxBf =
            gasMetalHeatFlux.boundaryFieldRef();
        const volScalarField::Boundary& rawFluxBf = rawFlux.boundaryField();

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
        tmp<volScalarField> tRho2Diag = mixture.thermo2().rho();
        tmp<volScalarField> tHeGasNew = mixture.thermo2().he();
        const volScalarField& rho2Diag = tRho2Diag();
        const volScalarField& heGasNew = tHeGasNew();

        dimensionedScalar newGasEnergy = fvc::domainIntegrate
        (
            alpha2*rho2Diag*(heGasNew - p/rho2Diag)
        );
        dimensionedScalar deltaGasEnergy = newGasEnergy - oldGasEnergy;

        dimensionedScalar gasMetalPower = fvc::domainIntegrate(gasMetalHeatFlux);
        dimensionedScalar exchangeEnergy = -gasMetalPower*runTime.deltaT();
        dimensionedScalar energyResidual = deltaGasEnergy - exchangeEnergy;

        scalar denom = Foam::max
        (
            Foam::max
            (
                Foam::mag(deltaGasEnergy.value()),
                Foam::mag(exchangeEnergy.value())
            ),
            VSMALL
        );
        scalar relResidual = Foam::mag(energyResidual.value())/denom;

        Info<< "Gas energy audit:" << nl
            << "  E_gas (old) = " << oldGasEnergy.value() << " J" << nl
            << "  E_gas (new) = " << newGasEnergy.value() << " J" << nl
            << "  ΔE_gas      = " << deltaGasEnergy.value() << " J" << nl
            << "  ∫Q_gm dt    = " << exchangeEnergy.value() << " J" << nl
            << "  Residual    = " << energyResidual.value() << " J" << nl
            << "  Relative residual = " << relResidual << nl
            << "  Gas-metal power   = " << gasMetalPower.value() << " W" << endl;

        tHeGasNew.clear();
        tRho2Diag.clear();
    }

    scalar couplingMin = gMin(interfaceExchange);
    scalar couplingMax = gMax(interfaceExchange);

    const dimensionedScalar interfaceExchangeFloor
    (
        "interfaceExchangeFloor",
        interfaceExchange.dimensions(),
        SMALL
    );

    tmp<volScalarField> tCouplingDelta =
        Foam::mag(gasMetalHeatFlux)
       /Foam::max(interfaceExchange, interfaceExchangeFloor);
    const volScalarField& couplingDelta = tCouplingDelta();
    scalar meanCouplingDelta = gAverage(couplingDelta);
    scalar maxCouplingDelta = gMax(couplingDelta);

    if (verbose)
    {
        Info<< "  Updated gas-metal heat flux range: ["
            << gMin(gasMetalHeatFlux) << ", "
            << gMax(gasMetalHeatFlux) << "] W/m³" << endl;

        Info<< "Interface coupling diagnostics:" << nl
            << "  Exchange coeff range: " << couplingMin << " - "
            << couplingMax << " W/m³/K" << nl
            << "  Heat flux range: " << gMin(gasMetalHeatFlux)
            << " - " << gMax(gasMetalHeatFlux) << " W/m³" << nl
            << "  Mean |T - Tl| (flux-derived): " << meanCouplingDelta << " K" << nl
            << "  Max |T - Tl| (flux-derived): " << maxCouplingDelta << " K" << endl;
    }

    tCouplingDelta.clear();
    tRawFlux.clear();
    // Validate thermodynamic consistency
    {
        const volScalarField& rho1Diag = mixture.thermo1().rho();
        const volScalarField& rho2Diag = mixture.thermo2().rho();
        // Check for any thermodynamic inconsistencies
        scalar minRho1 = min(rho1Diag).value();
        scalar minRho2 = min(rho2Diag).value();
        if (minRho1 <= 0 || minRho2 <= 0)
        {
            WarningInFunction
                << "Non-positive density detected after temperature update:" << nl
                << "  min(rho1) = " << minRho1 << " kg/m³" << nl
                << "  min(rho2) = " << minRho2 << " kg/m³" << nl
                << "Check thermodynamic properties" << endl;
        }
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
