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
    const dictionary& ctrl = runTime.controlDict();
    const scalar minPhysicalTemp = ctrl.lookupOrDefault<scalar>("minT", 250.0);
    const scalar maxPhysicalTemp = ctrl.lookupOrDefault<scalar>("maxT", 15000.0);
    const Switch enableTClamp = ctrl.lookupOrDefault<Switch>("enableTClamp", true);

    if (enableTClamp)
    {
        const dimensionedScalar minTempDim
        (
            "minT",
            dimTemperature,
            minPhysicalTemp
        );
        const dimensionedScalar maxTempDim
        (
            "maxT",
            dimTemperature,
            maxPhysicalTemp
        );

        T = Foam::max(T, minTempDim);
        T = Foam::min(T, maxTempDim);
        T.correctBoundaryConditions();
    }

    const volScalarField& TeField = ttm.Te();
    const volScalarField& TlField = ttm.Tl();
    volScalarField& gasMetalHeatFlux =
        mesh.lookupObjectRef<volScalarField>("gasMetalHeatFlux");
    tmp<volScalarField> tInterfaceExchange = ttm.gasMetalExchangeCoeffField();
    volScalarField& interfaceExchange = tInterfaceExchange.ref();
    const volScalarField& gasFraction = alpha2;
    const volScalarField& Qlaser = mixture.Q_laser();

    const scalar maxInterfaceExchange =
        ctrl.lookupOrDefault<scalar>("maxGasMetalExchangeCoeff", 1e18);
    const scalar minInterfaceExchange =
        ctrl.lookupOrDefault<scalar>("minGasMetalExchangeCoeff", 0.0);

    label interfaceClampedCells = 0;
    label interfaceClampedFaces = 0;

    {
        scalarField& coeffInternal = interfaceExchange.primitiveFieldRef();

        forAll(coeffInternal, cellI)
        {
            scalar& value = coeffInternal[cellI];

            if (value < minInterfaceExchange)
            {
                value = minInterfaceExchange;
                ++interfaceClampedCells;
            }
            else if (value > maxInterfaceExchange)
            {
                value = maxInterfaceExchange;
                ++interfaceClampedCells;
            }
        }

        volScalarField::Boundary& coeffBoundary =
            interfaceExchange.boundaryFieldRef();

        forAll(coeffBoundary, patchI)
        {
            scalarField& patchField = coeffBoundary[patchI];

            forAll(patchField, faceI)
            {
                scalar& value = patchField[faceI];

                if (value < minInterfaceExchange)
                {
                    value = minInterfaceExchange;
                    ++interfaceClampedFaces;
                }
                else if (value > maxInterfaceExchange)
                {
                    value = maxInterfaceExchange;
                    ++interfaceClampedFaces;
                }
            }
        }
    }

    reduce(interfaceClampedCells, sumOp<label>());
    reduce(interfaceClampedFaces, sumOp<label>());

    if (verbose && (interfaceClampedCells > 0 || interfaceClampedFaces > 0))
    {
        Info<< "Interface exchange coefficient limited in "
            << interfaceClampedCells << " cells and "
            << interfaceClampedFaces << " faces (range "
            << minInterfaceExchange << " - "
            << maxInterfaceExchange << " W/m³/K)" << endl;
    }


    tmp<volScalarField> tCv2 = mixture.thermo2().Cv();
    const volScalarField& Cv2 = tCv2();

    tmp<volScalarField> tRho2 = mixture.thermo2().rho();
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
    tmp<volScalarField> tCvGas(max(Cv2, cvFloor));
    const volScalarField& CvGas = tCvGas();

    tmp<volScalarField> tGasDensity(gasFraction*max(rho2, rhoFloor));
    const volScalarField& gasDensity = tGasDensity();

    tmp<volScalarField> tCvolGas
    (
        max(gasDensity, rhoFloor)*CvGas
    );
    const volScalarField& CvolGas = tCvolGas();

    tmp<surfaceScalarField> tGasCvolPhi
    (
        fvc::interpolate(CvolGas)*phi
    );
    const surfaceScalarField& gasCvolPhi = tGasCvolPhi();

    tmp<volScalarField> tContErr = fvc::div(gasCvolPhi);
    const volScalarField& contErr = tContErr();

    fvScalarMatrix TEqn
    (
        fvm::ddt(CvolGas, T)
      + fvm::div(gasCvolPhi, T)
      - fvm::Sp(contErr, T)
      - fvm::laplacian(transportModel.kappaEff(), T)
     ==
        fvOptions(CvolGas, T)
    );

    tmp<volScalarField> tCvolGasClamp
    (
        max(CvolGas, rhoFloor*cvFloor)
    );
    const volScalarField& CvolGasClamp = tCvolGasClamp();

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
    tCvGas.clear();
    tCv2.clear();
    tRho2.clear();

    tmp<volScalarField> tInterfaceRelaxCoeff
    (
        interfaceExchange/CvolGasClamp
    );
    const volScalarField& interfaceRelaxCoeff = tInterfaceRelaxCoeff();

    TEqn += fvm::Sp
    (
        interfaceRelaxCoeff,
        T
    );
    TEqn.source() += interfaceRelaxCoeff*TlField;
    TEqn.source() += alpha2*Qlaser/CvolGasClamp;

    tGasCvolPhi.clear();
    tGasDensity.clear();
    tInterfaceRelaxCoeff.clear();
    tCvolGas.clear();
    tCvolGasClamp.clear();
    TEqn.relax();
    fvOptions.constrain(TEqn);
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
    const scalar maxPhysicalFlux =
        ctrl.lookupOrDefault<scalar>("maxGasMetalHeatFlux", 1e15);

    label fluxLimitedCells = 0;
    label fluxLimitedFaces = 0;

    if (maxPhysicalFlux > SMALL)
    {
        scalarField& fluxInternal = gasMetalHeatFlux.primitiveFieldRef();

        forAll(fluxInternal, cellI)
        {
            scalar& value = fluxInternal[cellI];

            if (Foam::mag(value) > maxPhysicalFlux)
            {
                value = Foam::sign(value)*maxPhysicalFlux;
                ++fluxLimitedCells;
            }
        }
    }

    gasMetalHeatFlux.correctBoundaryConditions();

    if (maxPhysicalFlux > SMALL)
    {
        volScalarField::Boundary& fluxBoundary =
            gasMetalHeatFlux.boundaryFieldRef();

        forAll(fluxBoundary, patchI)
        {
            scalarField& patchField = fluxBoundary[patchI];

            forAll(patchField, faceI)
            {
                scalar& value = patchField[faceI];

                if (Foam::mag(value) > maxPhysicalFlux)
                {
                    value = Foam::sign(value)*maxPhysicalFlux;
                    ++fluxLimitedFaces;
                }               
            }
        }
        // Final safeguard using the controlDict-configured limit
        scalarField& fluxInternal = gasMetalHeatFlux.primitiveFieldRef();

        forAll(fluxInternal, cellI)
        {
            scalar& value = fluxInternal[cellI];

            if (Foam::mag(value) > maxPhysicalFlux)
            {
                value = Foam::sign(value)*maxPhysicalFlux;
            }
        }        
    }

    reduce(fluxLimitedCells, sumOp<label>());
    reduce(fluxLimitedFaces, sumOp<label>());

    if (verbose && (fluxLimitedCells > 0 || fluxLimitedFaces > 0))
    {
        Info<< "Gas-metal heat flux limited in " << fluxLimitedCells
            << " cells and " << fluxLimitedFaces << " faces (|q| ≤ "
            << maxPhysicalFlux << " W/m³)" << endl;
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
