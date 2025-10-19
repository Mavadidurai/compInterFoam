/*---------------------------------------------------------------------------*\
  Temperature equation for LIFT simulation with the two-temperature model
  Fixed version with proper gas-metal thermal coupling
\*---------------------------------------------------------------------------*/

extern Foam::Switch verbose;
extern const bool master;
{
    if (verbose && master)
    {
        Info<< "Solving temperature equation" << endl;
    }

    const dictionary& ctrl = runTime.controlDict();
    const Switch enableTClamp = ctrl.lookupOrDefault<Switch>("enableTClamp", false);
    const scalar minPhysicalTemp = Foam::max
    (
        ctrl.lookupOrDefault<scalar>("minT", 0.0),
        Foam::SMALL
    );
    const scalar maxPhysicalTemp = ctrl.lookupOrDefault<scalar>
    (
        "maxT",
        Foam::GREAT
    );
    const dimensionedScalar minTDim
    (
        "minT",
        dimTemperature,
        minPhysicalTemp
    );

    const dimensionedScalar maxTDim
    (
        "maxT",
        dimTemperature,
        maxPhysicalTemp
    );
    if (enableTClamp && maxPhysicalTemp <= minPhysicalTemp)
    {
        FatalErrorInFunction
            << "Temperature clamp bounds invalid: minT=" << minPhysicalTemp
            << " must be less than maxT=" << maxPhysicalTemp
            << abort(FatalError);
    }
    
    const scalar thermalFluxRelax =
        ctrl.lookupOrDefault<scalar>("thermalFluxRelax", 1.0);
    
    scalar totalCellsClamped = 0.0;
    scalar totalLowClamped = 0.0;
    scalar totalHighClamped = 0.0;
    
    auto clampTemperature =
        [&](const word& stage)
        {
            const scalarField TintPreClamp(T.internalField());
            label nLow = 0;
            label nHigh = 0;
            forAll(TintPreClamp, cellI)
            {
                const scalar val = TintPreClamp[cellI];
                if (val < minPhysicalTemp)
                {
                    ++nLow;
                }
                else if (val > maxPhysicalTemp)
                {
                    ++nHigh;
                }
            }

            scalar localLow = scalar(nLow);
            scalar localHigh = scalar(nHigh);
            scalar localCells = scalar(TintPreClamp.size());
            Foam::reduce(localLow, sumOp<scalar>());
            Foam::reduce(localHigh, sumOp<scalar>());
            Foam::reduce(localCells, sumOp<scalar>());
            totalLowClamped = localLow;
            totalHighClamped = localHigh;
            totalCellsClamped = localCells;

            T = Foam::max(Foam::min(T, maxTDim), minTDim);
            T.correctBoundaryConditions();

            if (Pstream::master())
            {
                if (localLow > 0 || localHigh > 0)
                {
                    WarningInFunction
                        << "Temperature clamp activated during " << stage
                        << ": enforcing [" << minPhysicalTemp << ", "
                        << maxPhysicalTemp << "] K" << endl;
                }
                else
                {
                    Info<< "T within physical bounds during " << stage << endl;
                }
            }
        };

    if (enableTClamp)
    {
        clampTemperature("pre-solve");
    }

    volScalarField& gasMetalHeatFlux =
        mesh.lookupObjectRef<volScalarField>("gasMetalHeatFlux");

    tmp<volScalarField> tInterfaceExchange = ttm.gasMetalExchangeCoeffField();
    volScalarField& interfaceExchange = tInterfaceExchange.ref();
    tmp<volScalarField> tCappedInterfaceExchange(interfaceExchange.clone());
    volScalarField& cappedInterfaceExchange = tCappedInterfaceExchange.ref();
    
    const dictionary twoTemperatureDict =
        ctrl.subOrEmptyDict("twoTemperatureProperties");

    const Switch userGasLaserHeating
    (
        twoTemperatureDict.lookupOrDefault<Switch>
        (
            "applyGasLaserHeating",
            false
        )
    );

    const bool laserSupportsGasHeating =
        laser.gasAbsorptionCoeff().value() > VSMALL;

    const bool applyGasLaserHeating =
        userGasLaserHeating && laserSupportsGasHeating;

    if (userGasLaserHeating && !laserSupportsGasHeating && verbose && master)
    {
        WarningInFunction
            << "applyGasLaserHeating requested but gasAbsorptionCoeff ≈ 0;"
            << " skipping direct gas heating."
            << endl;
    }


    const dimensionedScalar minCoeffDefault
    (
        "minGasMetalExchangeCoeff",
        interfaceExchange.dimensions(),
        ctrl.lookupOrDefault<scalar>("minGasMetalExchangeCoeff", 0.0)
    );

    const dimensionedScalar maxCoeffDefault
    (
        "maxGasMetalExchangeCoeff",
        interfaceExchange.dimensions(),
        ctrl.lookupOrDefault<scalar>("maxGasMetalExchangeCoeff", Foam::GREAT)
    );

    const dimensionedScalar minCoeff =
        twoTemperatureDict.lookupOrDefault<dimensionedScalar>
        (
            "minGasMetalExchangeCoeff",
            minCoeffDefault
        );

    const dimensionedScalar maxCoeff =
        twoTemperatureDict.lookupOrDefault<dimensionedScalar>
        (
            "maxGasMetalExchangeCoeff",
            maxCoeffDefault
        );
        
    cappedInterfaceExchange =
        Foam::min
        (
            Foam::max(cappedInterfaceExchange, minCoeff),
            maxCoeff
        );
    cappedInterfaceExchange.correctBoundaryConditions();
    
    const volScalarField& Tl = ttm.Tl();

    const scalar qCap =
        twoTemperatureDict.lookupOrDefault<scalar>
        (
            "maxGasMetalHeatFlux",
            -1.0
        );

    const dimensionedScalar qCapDim
    (
        "maxGasMetalHeatFlux",
        dimPower/dimArea,
        qCap
    );
    
    const label timeIndex = runTime.timeIndex();
    const dimensionedScalar maxTlDim = gMax(Tl);
    const scalar maxTl = maxTlDim.value();
    const dimensionedScalar dt = runTime.deltaT();
    
    const scalar metalFractionFloor =
        twoTemperatureDict.lookupOrDefault<scalar>
        (
            "metalFractionFloor",
            1e-6
        );
    const scalar blendWidth =
        twoTemperatureDict.lookupOrDefault<scalar>
        (
            "metalAmbientBlendWidth",
            1e-3
        );

    const scalar activeCutoff = Foam::max(metalFractionFloor, VSMALL);
    const dimensionedScalar cutoffDim("metalCutoff", dimless, activeCutoff);

    const volScalarField& metalFraction = mixture.alpha1();

    tmp<volScalarField> tActiveMask;
    if (blendWidth <= SMALL)
    {
        tActiveMask = Foam::pos(metalFraction - cutoffDim);
    }
    else
    {
        const scalar lowerBoundValue =
            Foam::max(metalFractionFloor - blendWidth, scalar(0));
        const dimensionedScalar lowerBound
        (
            "metalLowerBound",
            dimless,
            lowerBoundValue
        );
        const dimensionedScalar blendDim
        (
            "metalBlendWidth",
            dimless,
            Foam::max(blendWidth, SMALL)
        );

        tActiveMask =
            Foam::min
            (
                dimensionedScalar("one", dimless, 1.0),
                Foam::max
                (
                    dimensionedScalar("zero", dimless, 0.0),
                    (metalFraction - lowerBound)/blendDim
                )
            );
    }
if (verbose && master && timeIndex % 10 == 0)
{
    // Calculate cooling rate
    static scalar lastMaxTl = maxTl;
    const scalar coolingRate = (lastMaxTl - maxTl) / dt.value();

    Info<< "Thermal diagnostics:" << nl
        << "  Cooling rate: " << coolingRate << " K/s" << nl
        << "  Expected equilibrium time: "
        << (maxTl - 300.0) / max(coolingRate, SMALL) << " s" << endl;

    lastMaxTl = maxTl;
}
    volScalarField& activeMask = tActiveMask.ref();
    activeMask = Foam::min(Foam::max(activeMask, scalar(0)), scalar(1));

    const volScalarField& gasFraction = mixture.alpha2();

    tmp<volScalarField> tGasMask
    (
        Foam::min(Foam::max(gasFraction, scalar(0)), scalar(1))
    );
    const volScalarField& gasMask = tGasMask();

    tmp<volScalarField> tNeighbourGas(fvc::average(gasFraction));
    volScalarField& neighbourGas = tNeighbourGas.ref();
    neighbourGas = Foam::min(Foam::max(neighbourGas, scalar(0)), scalar(1));

    tmp<volScalarField> tContactMask(Foam::max(gasMask, neighbourGas));
    volScalarField& contactMask = tContactMask.ref();
    contactMask = Foam::min(Foam::max(contactMask, scalar(0)), scalar(1));
     contactMask.correctBoundaryConditions();

    tmp<volScalarField> tNeighbourMetal(fvc::average(metalFraction));
    volScalarField& neighbourMetal = tNeighbourMetal.ref();
    neighbourMetal = Foam::min
    (
        Foam::max(neighbourMetal, scalar(0)),
        scalar(1)
    );

    tmp<volScalarField> tInterfaceMask
    (
        Foam::max(activeMask, neighbourMetal)
    );
    volScalarField& interfaceMask = tInterfaceMask.ref();
    interfaceMask = Foam::min(Foam::max(interfaceMask, scalar(0)), scalar(1));
    interfaceMask = Foam::min
    (
        Foam::max(interfaceMask*contactMask, scalar(0)),
        scalar(1)
    );
    interfaceMask.correctBoundaryConditions();

    // Limit the gas-metal exchange to cells where gas contacts metal. Include
    // a neighbour-weighted metal mask so pure-gas cells adjacent to metal
    // retain a non-zero coupling weight.
    cappedInterfaceExchange *= interfaceMask;
    cappedInterfaceExchange.correctBoundaryConditions();

    tmp<volVectorField> tGradMetal(fvc::grad(metalFraction));
    const volVectorField& gradMetal = tGradMetal();

    tmp<volScalarField> tInterfaceDelta
    (
        new volScalarField
        (
            IOobject
            (
                "gasMetalInterfaceThickness",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar
            (
                "gasMetalInterfaceThickness",
                dimLength,
                1e-9
            )
        )
    );
    volScalarField& interfaceDelta = tInterfaceDelta.ref();

    scalarField& deltaInternal = interfaceDelta.primitiveFieldRef();

    forAll(deltaInternal, cellI)
    {
        const scalar cellVolume = mesh.V()[cellI];
        scalar interfaceArea = mag(gradMetal[cellI])*cellVolume;
        const scalar geometricArea = Foam::pow(cellVolume, scalar(2.0/3.0));

        if (interfaceArea <= SMALL)
        {
            interfaceArea = Foam::max(geometricArea, SMALL);
        }

        deltaInternal[cellI] = Foam::max(cellVolume/interfaceArea, SMALL);
    }

    interfaceDelta.correctBoundaryConditions();

    tmp<volScalarField> tSurfaceExchange(cappedInterfaceExchange*interfaceDelta);
    volScalarField& surfaceExchange = tSurfaceExchange.ref();

    if (qCap > 0)
    {
        const scalar tempDiffFloorValue =
            twoTemperatureDict.lookupOrDefault<scalar>
            (
                "temperatureDifferenceFloor",
                1e-3
            );

        const dimensionedScalar tempDiffFloor
        (
            "temperatureDifferenceFloor",
            dimTemperature,
            Foam::max(tempDiffFloorValue, SMALL)
        );

        tmp<volScalarField> tTempDiff(mag(T - Tl));
        volScalarField& tempDiff = tTempDiff.ref();
        tempDiff = Foam::max(tempDiff, tempDiffFloor);

        tmp<volScalarField> tCoeffCap(qCapDim/tempDiff);
        const volScalarField& coeffCap = tCoeffCap();

        surfaceExchange =
            Foam::min
            (
                Foam::max(surfaceExchange, -coeffCap),
                coeffCap
            );
        surfaceExchange.correctBoundaryConditions();
    }

    cappedInterfaceExchange = surfaceExchange/interfaceDelta;
    cappedInterfaceExchange.correctBoundaryConditions();

    tmp<volScalarField> tCv2 = mixture.thermo2().Cv();
    tmp<volScalarField> tRho2 = mixture.thermo2().rho();

    const scalar cvFloorValue =
        twoTemperatureDict.lookupOrDefault<scalar>("CvGasFloor", 100.0);

    const dimensionedScalar cvFloor
    (
        "cvFloor",
        dimEnergy/dimMass/dimTemperature,
        cvFloorValue
    );
    const dimensionedScalar rhoFloor("rhoFloor", dimDensity, 1e-12);

    tmp<volScalarField> tCvGas(Foam::max(tCv2(), cvFloor));
    const volScalarField& CvGas = tCvGas();

    tmp<volScalarField> tRhoGasCore(Foam::max(tRho2(), rhoFloor));
    const volScalarField& rhoGasCore = tRhoGasCore();

    tmp<volScalarField> tCvolGas(alpha2*rhoGasCore*CvGas);

    const volScalarField& CvolGas = tCvolGas();
    const scalar cvolFloorValue =
        twoTemperatureDict.lookupOrDefault<scalar>
        (
            "CvolGasFloor",
            1e-12
        );

    const dimensionedScalar cvolFloor
    (
        "CvolGasFloor",
        dimEnergy/dimVolume/dimTemperature,
        cvolFloorValue
    );

    tmp<volScalarField> tCvolLimited(Foam::max(CvolGas, cvolFloor));
    const volScalarField& CvolLimited = tCvolLimited();

    const scalar maxThermalCourant =
        Foam::max
        (
            scalar(0),
            Foam::min
            (
                ctrl.lookupOrDefault<scalar>("maxThermalCourant", 1.0),
                scalar(1)
            )
        );

    if (maxThermalCourant > SMALL && maxThermalCourant < 1.0 - SMALL)
    {
        const scalar limiterAlphaMin =
            Foam::max
            (
                ctrl.lookupOrDefault<scalar>("thermalLimiterAlphaMin", 1e-6),
                0.0
            );

        const dimensionedScalar alphaLimiterFloor
        (
            "thermalLimiterAlphaMin",
            dimless,
            limiterAlphaMin
        );

        tmp<volScalarField> tThermalCo
        (
            (cappedInterfaceExchange*runTime.deltaT())/CvolLimited
        );
        volScalarField& thermalCo = tThermalCo.ref();
        thermalCo = Foam::max
        (
            thermalCo,
            dimensionedScalar("zeroThermalCo", dimless, 0.0)
        );

        const dimensionedScalar maxThermalCourantDim
        (
            "maxThermalCourant",
            dimless,
            maxThermalCourant
        );
        const dimensionedScalar oneDim("one", dimless, 1.0);
        const dimensionedScalar epsCo
        (
            "epsThermalCourant",
            dimless,
            SMALL
        );

        tmp<volScalarField> tLimiter
        (
            maxThermalCourantDim/(thermalCo + epsCo)
        );
        volScalarField& limiter = tLimiter.ref();
        limiter = Foam::min(limiter, oneDim);
        limiter = Foam::max
        (
            limiter,
            dimensionedScalar("zeroLimiter", dimless, 0.0)
        );

        tmp<volScalarField> tLimiterActive
        (
            Foam::pos(alpha2 - alphaLimiterFloor)
        );
        const volScalarField& limiterActive = tLimiterActive();

        limiter = limiterActive*limiter + (oneDim - limiterActive)*oneDim;

        cappedInterfaceExchange *= limiter;
        cappedInterfaceExchange.correctBoundaryConditions();
    }

    interfaceExchange = cappedInterfaceExchange;
    interfaceExchange.correctBoundaryConditions();

    tmp<surfaceScalarField> tCvolPhi(fvc::interpolate(CvolLimited)*phi);
    const surfaceScalarField& CvolPhi = tCvolPhi();
    tmp<volScalarField> tContErr
    (
        fvc::ddt(CvolLimited)
      + fvc::div(CvolPhi)
    );
    volScalarField& contErr = tContErr.ref();

    const dimensionedScalar zeroContErr
    (
        "zeroContErr",
        contErr.dimensions(),
        0.0
    );

    contErr = Foam::min(contErr, zeroContErr);

    T = Foam::max(Foam::min(T, maxTDim), minTDim);
    T.correctBoundaryConditions();

    fvScalarMatrix TEqn
    (
        fvm::ddt(CvolLimited, T)
      + fvm::div(CvolPhi, T)
      - fvm::Sp(contErr, T)
      - fvm::laplacian(transportModel.kappaEff(), T)
      ==
        fvOptions(CvolLimited, T)
    );


    TEqn += fvm::Sp(cappedInterfaceExchange, T);
    TEqn.source() += cappedInterfaceExchange*Tl;


    if (applyGasLaserHeating)
    {
        const scalar metalFractionCutoff =
            twoTemperatureDict.lookupOrDefault<scalar>
            (
                "metalFractionCutoff",
                twoTemperatureDict.lookupOrDefault<scalar>
                (
                    "metalFractionFloor",
                    1e-6
                )
            );

        const dimensionedScalar metalCutoffDim
        (
            "metalFractionCutoff",
            dimless,
            Foam::max(metalFractionCutoff, Foam::SMALL)
        );

        const volScalarField& metalFraction = mixture.alpha1();

        tmp<volScalarField> tGasLaserMask
        (
            Foam::pos(metalCutoffDim - metalFraction)
        );

        const volScalarField& gasLaserMask = tGasLaserMask();

        TEqn.source() += contactMask*gasLaserMask*mixture.Q_laser();
    }

    TEqn.relax();
    fvOptions.constrain(TEqn);

    const dictionary solvers(mesh.solutionDict().subOrEmptyDict("solvers"));
    const word finalName(T.name() + "Final");

    if (pimple.finalIter() && !solvers.found(finalName))
    {
        dictionary fallbackSolver(mesh.solverDict(T.name()));
        TEqn.solve(fallbackSolver);
    }
    else
    {
        TEqn.solve();
    }

    fvOptions.correct(T);

    if (enableTClamp)
    {
        clampTemperature("post-solve");
    }
    else
    {
        T = Foam::max(T, minTDim);
        T.correctBoundaryConditions();
    }
    // Automatically enforce physical gas-temperature bounds when the
    // solution produces extreme or non-finite values.  Leaving these
    // unchecked causes the gas-metal heat flux to explode, which in turn
    // triggers the energy audit failure observed in the run log.
    const scalar Tmin = gMin(T);
    const scalar Tmax = gMax(T);

    const bool invalidExtremes =
        (!std::isfinite(Tmin) || !std::isfinite(Tmax));
    const bool violatesClampBounds =
        (Tmin < minPhysicalTemp || Tmax > maxPhysicalTemp);

    if (invalidExtremes || violatesClampBounds)
    {
        if (verbose && master)
        {
            WarningInFunction
                << "Gas temperature produced "
                << (invalidExtremes ? "non-finite" : "out-of-range")
                << " values: Tmin=" << Tmin
                << " K, Tmax=" << Tmax << " K.";

            if (enableTClamp)
            {
                Info<< " Enforcing clamp [" << minPhysicalTemp
                    << ", " << maxPhysicalTemp << "] K" << endl;
            }
            else
            {
                Info<< " Resetting to fallback bounds to recover stability" << endl;
            }
        }

        T = Foam::max(Foam::min(T, maxTDim), minTDim);
        T.correctBoundaryConditions();
    }

    const scalar fluxRelax =
        Foam::max(scalar(0), Foam::min(scalar(1), thermalFluxRelax));

    tmp<volScalarField> tRawFlux(cappedInterfaceExchange*(Tl - T));
    volScalarField& rawFlux = tRawFlux.ref();

    tmp<volScalarField> tSurfaceFlux(rawFlux*interfaceDelta);
    volScalarField& surfaceFlux = tSurfaceFlux.ref();

    if (qCap > 0)
    {
        surfaceFlux = Foam::min(Foam::max(surfaceFlux, -qCapDim), qCapDim);
        surfaceFlux.correctBoundaryConditions();
    }

    rawFlux = surfaceFlux/interfaceDelta;

    const volScalarField& limitedFlux = tRawFlux();

    tmp<volScalarField> tMaskedFlux(new volScalarField(limitedFlux));
    volScalarField& maskedFlux = tMaskedFlux.ref();

    const volScalarField& metalToGasFlux = tMaskedFlux();

    if (fluxRelax < 1.0 - SMALL)
    {
        gasMetalHeatFlux =
            fluxRelax*metalToGasFlux + (1.0 - fluxRelax)*gasMetalHeatFlux;
    }
    else
    {
        gasMetalHeatFlux = metalToGasFlux;
    }

    if (qCap > 0)
    {
        tmp<volScalarField> tRelaxedSurfaceFlux(gasMetalHeatFlux*interfaceDelta);
        volScalarField& relaxedSurfaceFlux = tRelaxedSurfaceFlux.ref();

        relaxedSurfaceFlux =
            Foam::min
            (
                Foam::max(relaxedSurfaceFlux, -qCapDim),
                qCapDim
            );
        relaxedSurfaceFlux.correctBoundaryConditions();

        gasMetalHeatFlux = relaxedSurfaceFlux/interfaceDelta;
    }

    gasMetalHeatFlux.correctBoundaryConditions();

    if (verbose && master)
    {
        Info<< "    min(T) = " << min(T).value()
            << "  max(T) = " << max(T).value()
            << "  mean(T) = " << gAverage(T) << endl;

        tmp<volScalarField> tDiagSurfaceFlux(gasMetalHeatFlux*interfaceDelta);
        const volScalarField& diagSurfaceFlux = tDiagSurfaceFlux();

        Info<< "    mean exchange coeff (vol) = "
            << gAverage(interfaceExchange)
            << "  mean q_gm (surface) = "
            << gAverage(diagSurfaceFlux) << endl;
        Info<< "    max |q_gm| (surface) = "
            << gMax(mag(diagSurfaceFlux)().internalField())
            << "  max |q_gm| (vol) = "
            << gMax(mag(gasMetalHeatFlux)().internalField()) << endl;

        // Diagnostic: Total heat flux from metal to gas
        const dimensionedScalar totalHeatOut =
            fvc::domainIntegrate(gasMetalHeatFlux);
        Info<< "    net heat flux metal→gas: "
            << totalHeatOut.value() << " W" << endl;

        if (enableTClamp)
        {
            const scalar denom = Foam::max(scalar(1), totalCellsClamped);
            Info<< "    clamp fractions: low=" << totalLowClamped/denom
                << "  high=" << totalHighClamped/denom << endl;
        }
    }

    // Enforce lattice temperature only in cells that are overwhelmingly metal.
    const scalar strictMetalThreshold = 0.9;
    const dimensionedScalar strictMetalDim
    (
        "strictMetalThreshold",
        dimless,
        strictMetalThreshold
    );

    tmp<volScalarField> tPureMetal
    (
        Foam::pos(metalFraction - strictMetalDim)
    );
    const volScalarField& pureMetal = tPureMetal();

    T = (scalar(1) - pureMetal)*T + pureMetal*Tl;
    if (enableTClamp)
    {
        T = Foam::max(Foam::min(T, maxTDim), minTDim);
    }
    else
    {
        T = Foam::max(T, minTDim);
    }
    T.correctBoundaryConditions();

    mixture.correctThermo();
    mixture.correct();
}

// ************************************************************************* //
