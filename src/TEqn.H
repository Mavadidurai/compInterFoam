/*---------------------------------------------------------------------------*\
  Temperature equation for LIFT simulation with the two-temperature model
  Fixed version with proper gas-metal thermal coupling
\*---------------------------------------------------------------------------*/

extern Foam::Switch verbose;
{
    if (verbose && Foam::Pstream::master())
    {
        Info<< "Solving temperature equation" << endl;
    }

    const dictionary& ctrl = runTime.controlDict();
    const Switch femtosecondShot =
        ctrl.lookupOrDefault<Switch>("femtosecondShot", false);

    // The option enableTClamp (an elaborate replacement of out-of-range gas
    // temperatures by neighbour averages, about 400 lines) was removed: no
    // case used it and it had never been tested. The gas temperature is
    // bounded to [minT, maxT] as before.
    if (ctrl.lookupOrDefault<Switch>("enableTClamp", false))
    {
        FatalErrorInFunction
            << "enableTClamp true is no longer supported: the temperature"
            << " projection was removed. The gas temperature is bounded to"
            << " [minT, maxT]; remove the entry or set it to false."
            << exit(FatalError);
    }

    const scalar minPhysicalTemp = Foam::max
    (
        ctrl.lookupOrDefault<scalar>("minT", 0.0),
        Foam::SMALL
    );
    scalar maxPhysicalTemp = ctrl.lookupOrDefault<scalar>
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

    const dictionary twoTemperatureDict =
        ctrl.subOrEmptyDict("twoTemperatureProperties");

    const dimensionedScalar maxTDim
    (
        "maxT",
        dimTemperature,
        maxPhysicalTemp
    );
    
    const scalar thermalFluxRelax =
        ctrl.lookupOrDefault<scalar>("thermalFluxRelax", 1.0);
    const scalar fluxRelax =
        Foam::max(scalar(0), Foam::min(scalar(1), thermalFluxRelax));
    const volScalarField& Tl = ttm.Tl();

    volScalarField& gasMetalHeatFlux =
        mesh.lookupObjectRef<volScalarField>("gasMetalHeatFlux");

    tmp<volScalarField> tInterfaceExchange = ttm.gasMetalExchangeCoeffField();
    volScalarField& interfaceExchange = tInterfaceExchange.ref();
    tmp<volScalarField> tCappedInterfaceExchange(interfaceExchange.clone());
    volScalarField& cappedInterfaceExchange = tCappedInterfaceExchange.ref();
    static bool gasMetalExchangeLimiterReported = false;
    
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

    if
    (
        userGasLaserHeating
     && !laserSupportsGasHeating
     && verbose
     && Foam::Pstream::master()
    )
    {
        WarningInFunction
            << "applyGasLaserHeating requested but gasAbsorptionCoeff ~ 0;"
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

    const bool limitInterfaceExchange =
        maxCoeff.value() > SMALL && maxCoeff.value() < Foam::GREAT;

    if
    (
        !gasMetalExchangeLimiterReported
     && limitInterfaceExchange
    )
    {
        const dimensionedScalar unboundedMax = gMax(interfaceExchange);
        if
        (
            unboundedMax.value() > maxCoeff.value()*(1.0 + SMALL)
        )
        {
            Info<< "Limiting gas-metal exchange coefficient to "
                << maxCoeff.value() << " W/m^3/K at time "
                << runTime.timeName() << " (peak="
                << unboundedMax.value() << ")." << endl;
            gasMetalExchangeLimiterReported = true;
        }
    }
    if (limitInterfaceExchange)
    {
        cappedInterfaceExchange =
            Foam::min
            (
                Foam::max(cappedInterfaceExchange, minCoeff),
                maxCoeff
            );
    }
    else
    {
        cappedInterfaceExchange = Foam::max(cappedInterfaceExchange, minCoeff);
    }
    cappedInterfaceExchange.correctBoundaryConditions();
    
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
    const bool applyHeatFluxCap = qCap > 0 && !femtosecondShot;
    
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
    const volScalarField& gasFraction = mixture.alpha2();

    tmp<volScalarField> tGasMask
    (
        Foam::min(Foam::max(gasFraction, scalar(0)), scalar(1))
    );
    const volScalarField& gasMask = tGasMask();
    const scalar gasFractionFloorValue = twoTemperatureDict.lookupOrDefault<scalar>
    (
        "gasFractionFloor",
        metalFractionFloor
    );

    const dimensionedScalar gasMaskFloorDim
    (
        "gasFractionFloor",
        dimless,
        Foam::max(gasFractionFloorValue, VSMALL)
    );

    tmp<volScalarField> tGasCapacityMask
    (
        Foam::max(gasMask, gasMaskFloorDim)
    );
    const volScalarField& gasCapacityMask = tGasCapacityMask();
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
if (verbose && Foam::Pstream::master() && timeIndex % 10 == 0)
{
    // Calculate cooling rate
    static scalar lastMaxTl = maxTl;
    static scalar lastSampleTime = runTime.value();
    const scalar currentTime = runTime.value();
    const scalar elapsed = currentTime - lastSampleTime;

    Info<< "Thermal diagnostics:" << nl;

    if (elapsed > SMALL)
    {
        const scalar coolingRate = (lastMaxTl - maxTl)/elapsed;

        Info<< "  Cooling rate: " << coolingRate << " K/s" << nl;

        if (coolingRate > 0)
        {
            Info<< "  Expected equilibrium time: "
                << (maxTl - 300.0)/coolingRate << " s" << endl;
        }
        else
        {
            Info<< "  Cooling rate = 0 (heating or steady); skipping equilibrium"
                " estimate." << endl;
        }
    }
    else
    {
        Info<< "  Cooling rate unavailable (elapsed time too small)." << endl;
    }

    lastMaxTl = maxTl;
    lastSampleTime = currentTime;
}
    volScalarField& activeMask = tActiveMask.ref();
     activeMask = Foam::min(Foam::max(activeMask, scalar(0)), scalar(1));

    tmp<volVectorField> tGradMetal(fvc::grad(metalFraction));
    const volVectorField& gradMetal = tGradMetal();

    const scalar interfaceGradientCutoff =
        twoTemperatureDict.lookupOrDefault<scalar>
        (
            "interfaceGradientCutoff",
            Foam::SMALL
        );

    const scalar interfaceMetalLowerBound =
        twoTemperatureDict.lookupOrDefault<scalar>
        (
            "interfaceMetalLowerBound",
            metalFractionFloor
        );

    const scalar interfaceMetalUpperBound =
        twoTemperatureDict.lookupOrDefault<scalar>
        (
            "interfaceMetalUpperBound",
            scalar(1)
        );

    const scalar boundedMetalLower =
        Foam::min(Foam::max(interfaceMetalLowerBound, scalar(0)), scalar(1));
    const scalar boundedMetalUpper =
        Foam::max
        (
            boundedMetalLower,
            Foam::min(Foam::max(interfaceMetalUpperBound, scalar(0)), scalar(1))
        );

    tmp<volScalarField> tInterfaceMask(activeMask.clone());
    volScalarField& interfaceMask = tInterfaceMask.ref();
    scalarField& interfaceMaskInternal = interfaceMask.primitiveFieldRef();
    interfaceMaskInternal = scalar(0);

    forAll(interfaceMaskInternal, cellI)
    {
        const scalar mf = metalFraction[cellI];
        const scalar gradMag = mag(gradMetal[cellI]);

        if
        (
            mf >= boundedMetalLower
         && mf <= boundedMetalUpper + Foam::max(SMALL, metalFractionFloor)
         && gradMag >= interfaceGradientCutoff
        )
        {
            interfaceMaskInternal[cellI] = scalar(1);
        }
    }

    interfaceMask.correctBoundaryConditions();

    // Limit the gas-metal exchange strictly to the identified interface cells
    // so that pure-gas regions do not inherit a neighbour-averaged coupling.
    cappedInterfaceExchange = interfaceMask*cappedInterfaceExchange;
    cappedInterfaceExchange.correctBoundaryConditions();

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

    // (A multiply-by-delta then divide-by-delta round trip that changed
    // nothing was removed here.)
    cappedInterfaceExchange.correctBoundaryConditions();


    tmp<volScalarField> tMaxInterfaceExchange;
    const volScalarField* maxInterfaceExchangePtr = nullptr;

    if (applyHeatFluxCap)
    {
        // Use the actual (lagged, pre-solve) gas/lattice temperature
        // difference rather than a fictitious 1 K reference: qCap is a
        // physical W/m^2 flux cap, so converting it to a volumetric
        // exchange coefficient must divide by the real driving delta-T,
        // not an assumed 1 K one -- otherwise the cap only bounds the
        // coefficient for a 1 K difference and is decorative at the much
        // larger delta-T seen during a femtosecond pulse.
        const dimensionedScalar minReferenceDeltaT
        (
            "minReferenceDeltaT",
            dimTemperature,
            1.0
        );

        tmp<volScalarField> tReferenceDeltaT
        (
            Foam::max(mag(Tl - T), minReferenceDeltaT)
        );
        const volScalarField& referenceDeltaT = tReferenceDeltaT();

        tMaxInterfaceExchange = (qCapDim/referenceDeltaT)/interfaceDelta;
        const volScalarField& maxInterfaceExchange = tMaxInterfaceExchange();
        maxInterfaceExchangePtr = &maxInterfaceExchange;

        cappedInterfaceExchange =
            Foam::min
            (
                cappedInterfaceExchange,
                maxInterfaceExchange
            );
        cappedInterfaceExchange.correctBoundaryConditions();
    }

    tmp<volScalarField> tCv2 = mixture.thermo2().Cv();
    tmp<volScalarField> tRho2 = mixture.thermo2().rho();

    const scalar cvFloorValue =
        twoTemperatureDict.lookupOrDefault<scalar>("CvGasFloor", 1.0);

    const dimensionedScalar cvFloor
    (
        "cvFloor",
        dimEnergy/dimMass/dimTemperature,
        cvFloorValue
    );
    const scalar rhoFloorValue = Foam::max
    (
        twoTemperatureDict.lookupOrDefault<scalar>("rhoGasFloor", 1e-6),
        SMALL
    );

    const dimensionedScalar rhoFloor
    (
        "rhoFloor",
        dimDensity,
        rhoFloorValue
    );

    tmp<volScalarField> tCvGas(Foam::max(tCv2(), cvFloor));
    const volScalarField& CvGas = tCvGas();

    tmp<volScalarField> tRhoGasCore(Foam::max(tRho2(), rhoFloor));
    const volScalarField& rhoGasCore = tRhoGasCore();

    tmp<volScalarField> tCvolGas(gasFraction*rhoGasCore*CvGas);

    const volScalarField& CvolGas = tCvolGas();
    const scalar cvolFloorDefault = rhoFloorValue*cvFloorValue;

    const scalar cvolFloorValue = Foam::max
    (
        twoTemperatureDict.lookupOrDefault<scalar>
        (
            "CvolGasFloor",
            cvolFloorDefault
        ),
        SMALL
    );

    const dimensionedScalar cvolFloor
    (
        "CvolGasFloor",
        dimEnergy/dimVolume/dimTemperature,
        cvolFloorValue
    );

    tmp<volScalarField> tCvolFloorField(gasCapacityMask*cvolFloor);
    const volScalarField& CvolFloorField = tCvolFloorField();

    tmp<volScalarField> tCvolLimited
    (
        Foam::max
        (
            CvolGas,
            CvolFloorField
        )
    );
    const volScalarField& CvolLimited = tCvolLimited();
    const dimensionedScalar minCvol
    (
        "minCvol",
        CvolLimited.dimensions(),
        SMALL
    );

    // bound the gas temperature to [minT, maxT]
    auto projectTemperature =
        [&](const word&)
        {
            T = Foam::max(T, minTDim);
            T = Foam::min(T, maxTDim);
            T.correctBoundaryConditions();
        };

    projectTemperature("pre-solve");


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
            (cappedInterfaceExchange*runTime.deltaT())
          / Foam::max(CvolLimited, minCvol)
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
            Foam::pos(gasFraction - alphaLimiterFloor)
        );
        const volScalarField& limiterActive = tLimiterActive();

        limiter = limiterActive*limiter + (oneDim - limiterActive)*oneDim;

        cappedInterfaceExchange *= limiter;
        cappedInterfaceExchange.correctBoundaryConditions();
    }

    if (fluxRelax < 1.0 - SMALL)
    {
        tmp<volScalarField> tDeltaTmag(mag(Tl - T));
        volScalarField& deltaTmag = tDeltaTmag.ref();

        dimensionedScalar deltaTFloorInput
        (
            "temperatureDifferenceFloor",
            dimTemperature,
            10.0
        );

        if (twoTemperatureDict.found("temperatureDifferenceFloor"))
        {
            const dimensionedScalar rawFloor =
                twoTemperatureDict.lookupOrDefault<dimensionedScalar>
                (
                    "temperatureDifferenceFloor",
                    deltaTFloorInput
                );

            if (rawFloor.dimensions() == dimTemperature)
            {
                deltaTFloorInput = rawFloor;
            }
            else if (rawFloor.dimensions() == dimless)
            {
                deltaTFloorInput = dimensionedScalar
                (
                    rawFloor.name(),
                    dimTemperature,
                    rawFloor.value()
                );
            }
            else
            {
                FatalIOErrorInFunction(twoTemperatureDict)
                    << "Entry 'temperatureDifferenceFloor' has dimensions "
                    << rawFloor.dimensions()
                    << ", expected " << dimTemperature
                    << " or " << dimless
                    << exit(FatalIOError);
            }
        }

        const dimensionedScalar deltaTFloor
        (
            "deltaTFloor",
            dimTemperature,
            Foam::max(deltaTFloorInput.value(), scalar(SMALL))
        );

        deltaTmag = Foam::max(deltaTmag, deltaTFloor);

        tmp<volScalarField> tPrevCoeff(gasMetalHeatFlux/deltaTmag);
        volScalarField& prevCoeff = tPrevCoeff.ref();

        prevCoeff =
            Foam::min
            (
                Foam::max(prevCoeff, minCoeff),
                maxCoeff
            );

        if (maxInterfaceExchangePtr)
        {
            prevCoeff = Foam::min(prevCoeff, *maxInterfaceExchangePtr);
        }

        // Restrict the relaxation memory to current interface cells;
        // otherwise the minimum-coefficient floor re-introduces a
        // domain-wide gas-lattice coupling through the blending below.
        prevCoeff = interfaceMask*prevCoeff;

        prevCoeff.correctBoundaryConditions();

        cappedInterfaceExchange =
            fluxRelax*cappedInterfaceExchange
          + (scalar(1) - fluxRelax)*prevCoeff;

        if (maxInterfaceExchangePtr)
        {
            cappedInterfaceExchange =
                Foam::min
                (
                    cappedInterfaceExchange,
                    *maxInterfaceExchangePtr
                );
        }

        cappedInterfaceExchange.correctBoundaryConditions();
    }

    interfaceExchange = cappedInterfaceExchange;

    if (maxInterfaceExchangePtr)
    {
        interfaceExchange =
            Foam::min
            (
                interfaceExchange,
                *maxInterfaceExchangePtr
            );
    }

    interfaceExchange.correctBoundaryConditions();

    const volScalarField& relaxedInterfaceExchange = interfaceExchange;

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

    // Only the negative part of the continuity residual is kept, on purpose:
    // -fvm::Sp(contErr, T) with contErr < 0 adds to the matrix diagonal, while
    // a positive residual would subtract from it and can make the gas
    // temperature equation unbounded in the nearly empty cells (capacity at
    // its floor). The price is that the advective form is exact only where
    // the residual is negative. The gas temperature does not act back on the
    // flow, so this choice has been left as it is.
    contErr = Foam::min(contErr, zeroContErr);

    projectTemperature("pre-assembly");

    fvScalarMatrix TEqn
    (
        fvm::ddt(CvolLimited, T)
      + fvm::div(CvolPhi, T)
      - fvm::Sp(contErr, T)
      - fvm::laplacian(transportModel.kappaEff(), T)
      ==
        fvOptions(CvolLimited, T)
    );

    TEqn += fvm::Sp(relaxedInterfaceExchange, T);
    TEqn.source() += mesh.V()*relaxedInterfaceExchange*Tl;

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
        volScalarField& gasLaserMask = tGasLaserMask.ref();

        // contactMask (interfaceMask) requires metalFraction above the
        // interface floor, so contactMask*gasLaserMask was structurally
        // always ~zero: a cell cannot simultaneously be "gas-dominated"
        // (gasLaserMask) and "contain metal above the floor" (contactMask).
        // What this feature needs is GAS cells adjacent to the interface
        // (close enough to the metal surface for the laser to also deposit
        // energy directly into the gas), identified the same way
        // interfaceMask identifies metal-side interface cells: a nonzero
        // metal-fraction gradient.
        scalarField& gasLaserMaskInternal = gasLaserMask.primitiveFieldRef();

        forAll(gasLaserMaskInternal, cellI)
        {
            if
            (
                gasLaserMaskInternal[cellI] > SMALL
             && mag(gradMetal[cellI]) < interfaceGradientCutoff
            )
            {
                gasLaserMaskInternal[cellI] = 0.0;
            }
        }

        gasLaserMask.correctBoundaryConditions();

        TEqn.source() += mesh.V()*gasLaserMask*mixture.Q_laser();
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

    projectTemperature("post-solve");
    
    // Preserve the gas temperature prior to enforcing the metal projection so
    // that interfacial fluxes use the true post-solve gas state.
    tmp<volScalarField> tInterfaceTemperature(T.clone());
    const volScalarField& interfaceTemperature = tInterfaceTemperature();

    // Enforce lattice temperature only in cells that are overwhelmingly metal.
    const scalar strictMetalThreshold = 0.9;
    const dimensionedScalar strictMetalDim
    (
        "strictMetalThreshold",
        dimless,
        strictMetalThreshold
    );

    // Blend over the same width used for activeMask above instead of a
    // hard step at strictMetalThreshold, which otherwise creates a sharp T
    // discontinuity exactly at metalFraction = 0.9.
    tmp<volScalarField> tPureMetal;

    if (blendWidth <= SMALL)
    {
        tPureMetal = Foam::pos(metalFraction - strictMetalDim);
    }
    else
    {
        const scalar strictMetalLowerValue =
            Foam::max(strictMetalThreshold - blendWidth, scalar(0));
        const dimensionedScalar strictMetalLower
        (
            "strictMetalLowerBound",
            dimless,
            strictMetalLowerValue
        );
        const dimensionedScalar strictMetalBlendDim
        (
            "strictMetalBlendWidth",
            dimless,
            Foam::max(blendWidth, SMALL)
        );

        tPureMetal =
            Foam::min
            (
                dimensionedScalar("one", dimless, 1.0),
                Foam::max
                (
                    dimensionedScalar("zero", dimless, 0.0),
                    (metalFraction - strictMetalLower)/strictMetalBlendDim
                )
            );
    }

    const volScalarField& pureMetal = tPureMetal();

    T = (scalar(1) - pureMetal)*T + pureMetal*Tl;
    projectTemperature("post-metal");

    tmp<volScalarField> tRawFlux
    (
        relaxedInterfaceExchange*(Tl - interfaceTemperature)
    );
    volScalarField& rawFlux = tRawFlux.ref();

    rawFlux.correctBoundaryConditions();

    // Snapshot of the flux implied by the already-solved (unclamped) TEqn,
    // used below to feed any subsequent clamping back into T so that gas
    // never gains more energy than what gasMetalHeatFlux later removes
    // from the lattice equation.
    tmp<volScalarField> tOriginalRawFlux(rawFlux.clone());
    const volScalarField& originalRawFlux = tOriginalRawFlux();

    // Use VSMALL as minimum timestep to prevent division by near-zero values
    // that could produce unrealistic heating/cooling rates. VSMALL ~ 1e-300
    // provides better numerical stability than SMALL ~ 1e-15.
    const scalar dtVal = Foam::max(dt.value(), VSMALL);
    // Safety check: warn if timestep becomes dangerously small
    static bool timestepWarningIssued = false;
    if
    (
        dt.value() < 1e-15
     && !timestepWarningIssued
     && verbose
     && Foam::Pstream::master()
    )
    {
        WarningInFunction
            << "Time step (" << dt.value() << " s) is extremely small. "
            << "Flux limiting may produce very large heating/cooling rates." << nl
            << "Consider adjusting maxCo or time step controls." << endl;
        timestepWarningIssued = true;
    }
    scalarField& rawFluxInternal = rawFlux.primitiveFieldRef();
    const scalarField& capInternal = CvolLimited.internalField();
    const scalarField& TintInternal = T.internalField();

    forAll(rawFluxInternal, cellI)
    {
        const scalar cap = Foam::max(capInternal[cellI], VSMALL);
        const scalar Tcell = TintInternal[cellI];
        // Limit heat flux to prevent non-physical temperature changes
        // maxHeating: maximum energy input that won't exceed maxPhysicalTemp
        // maxCooling: maximum energy removal that won't drop below minPhysicalTemp
        const scalar maxHeating = cap*(maxPhysicalTemp - Tcell)/dtVal;
        const scalar maxCooling = cap*(Tcell - minPhysicalTemp)/dtVal;

        scalar limitedFlux = rawFluxInternal[cellI];

        if (limitedFlux > maxHeating)
        {
            limitedFlux = maxHeating;
        }
        else if (limitedFlux < -maxCooling)
        {
            limitedFlux = -maxCooling;
        }

        rawFluxInternal[cellI] = limitedFlux;
    }
    // Apply same flux limiting to boundary faces
    forAll(rawFlux.boundaryField(), patchI)
    {
        scalarField& fluxPatch = rawFlux.boundaryFieldRef()[patchI];
        const scalarField& capPatch = CvolLimited.boundaryField()[patchI];
        const scalarField& Tpatch = interfaceTemperature.boundaryField()[patchI];
        forAll(fluxPatch, faceI)
        {
            const scalar cap = Foam::max(capPatch[faceI], VSMALL);
            const scalar Tface = Tpatch[faceI];
            // Same limiting strategy as internal cells
            const scalar maxHeating = cap*(maxPhysicalTemp - Tface)/dtVal;
            const scalar maxCooling = cap*(Tface - minPhysicalTemp)/dtVal;

            scalar limitedFlux = fluxPatch[faceI];

            if (limitedFlux > maxHeating)
            {
                limitedFlux = maxHeating;
            }
            else if (limitedFlux < -maxCooling)
            {
                limitedFlux = -maxCooling;
            }

            fluxPatch[faceI] = limitedFlux;
        }
    }

    const scalar maxGasHeatingRate =
        twoTemperatureDict.lookupOrDefault<scalar>
        (
            "maxGasHeatingRate",
            -1.0
        );

    if (maxGasHeatingRate > SMALL)
    {
        // dtVal is already protected (max(dt.value(), VSMALL)), so direct division is safe
        const scalar invDtVal = 1.0/dtVal;

        const scalarField& gasMaskInternal = gasMask.internalField();
        const scalarField& cellVolumes = mesh.V();

        forAll(rawFluxInternal, cellI)
        {
            scalar& flux = rawFluxInternal[cellI];

            const scalar mask = Foam::max(gasMaskInternal[cellI], scalar(0));

            if (mask <= SMALL)
            {
                continue;
            }

            const scalar cellVolume = Foam::max(cellVolumes[cellI], SMALL);
            const scalar gasVolume = Foam::max(mask*cellVolume, SMALL);
            const scalar cap = Foam::max(capInternal[cellI], VSMALL);
            const scalar gasCapacity = cap*gasVolume;
            const scalar maxDeltaT = maxGasHeatingRate*dtVal;
            const scalar maxEnergyPerStep = maxDeltaT*gasCapacity;
            const scalar volumetricLimit =
                maxEnergyPerStep*invDtVal/Foam::max(gasVolume, SMALL);

            const scalar maxHeating = volumetricLimit;
            const scalar maxCooling = volumetricLimit;

            if (flux > maxHeating)
            {
                flux = maxHeating;
            }
            else if (flux < -maxCooling)
            {
                flux = -maxCooling;
            }
        }
    }
    
    tmp<volScalarField> tSurfaceFlux(rawFlux*interfaceDelta);
    volScalarField& surfaceFlux = tSurfaceFlux.ref();

    surfaceFlux.correctBoundaryConditions();

    if (qCap > 0)
    {
        surfaceFlux = Foam::min(Foam::max(surfaceFlux, -qCapDim), qCapDim);
        surfaceFlux.correctBoundaryConditions();
    }

    rawFlux = surfaceFlux/interfaceDelta;

    // The clamps above (per-cell rate cap, gas-heating-rate cap, qCap) may
    // have reduced |rawFlux| below what TEqn's implicit solve already used
    // (relaxedInterfaceExchange*(Tl - T)). Feed the shortfall back into T so
    // the gas only retains the energy that gasMetalHeatFlux will actually
    // remove from the lattice below -- otherwise every engaged clamp
    // creates energy from nothing.
    {
        scalarField& TInternal = T.primitiveFieldRef();
        const scalarField& finalFluxInternal = rawFlux.internalField();
        const scalarField& origFluxInternal = originalRawFlux.internalField();

        forAll(TInternal, cellI)
        {
            const scalar cap = Foam::max(capInternal[cellI], VSMALL);
            TInternal[cellI] +=
                (finalFluxInternal[cellI] - origFluxInternal[cellI])*dtVal/cap;
        }

        forAll(T.boundaryField(), patchI)
        {
            scalarField& TPatch = T.boundaryFieldRef()[patchI];
            const scalarField& capPatch = CvolLimited.boundaryField()[patchI];
            const scalarField& fluxPatch = rawFlux.boundaryField()[patchI];
            const scalarField& origFluxPatch =
                originalRawFlux.boundaryField()[patchI];

            forAll(TPatch, faceI)
            {
                const scalar cap = Foam::max(capPatch[faceI], VSMALL);
                TPatch[faceI] +=
                    (fluxPatch[faceI] - origFluxPatch[faceI])*dtVal/cap;
            }
        }

        T.correctBoundaryConditions();
    }

    const volScalarField& metalToGasFlux = tRawFlux();

    gasMetalHeatFlux = metalToGasFlux;

    gasMetalHeatFlux.correctBoundaryConditions();

    if (verbose && Foam::Pstream::master())
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
        Info<< "    net heat flux metal to gas: "
            << totalHeatOut.value() << " W" << endl;

    }

    mixture.correctThermo();
    mixture.correct();
}

// ************************************************************************* //
