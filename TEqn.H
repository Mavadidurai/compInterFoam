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

    if (enableTClamp && !ctrl.found("maxT"))
    {
        if (twoTemperatureDict.found("maxTl"))
        {
            maxPhysicalTemp = twoTemperatureDict.lookup<scalar>("maxTl");

            if (verbose && master)
            {
                WarningInFunction
                    << "enableTClamp active without controlDict maxT;"
                    << " falling back to twoTemperatureProperties.maxTl="
                    << maxPhysicalTemp
                    << endl;
            }
        }
        else
        {
            FatalErrorInFunction
                << "enableTClamp requires controlDict entry maxT or"
                << " twoTemperatureProperties.maxTl"
                << abort(FatalError);
        }
    }

    if (enableTClamp && maxPhysicalTemp <= minPhysicalTemp)
    {
        FatalErrorInFunction
            << "Temperature clamp bounds invalid: minT=" << minPhysicalTemp
            << " must be less than maxT=" << maxPhysicalTemp
            << abort(FatalError);
    }

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
    scalar totalCellsProjected = 0.0;
    scalar totalLowProjected = 0.0;
    scalar totalHighProjected = 0.0;

    volScalarField& gasMetalHeatFlux =
        mesh.lookupObjectRef<volScalarField>("gasMetalHeatFlux");

    tmp<volScalarField> tInterfaceExchange = ttm.gasMetalExchangeCoeffField();
    volScalarField& interfaceExchange = tInterfaceExchange.ref();
    tmp<volScalarField> tCappedInterfaceExchange(interfaceExchange.clone());
    volScalarField& cappedInterfaceExchange = tCappedInterfaceExchange.ref();
    
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
    const bool applyHeatFluxCap = qCap > 0;
    
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

    // Limit the gas-metal exchange strictly to the identified interface cells
    // so that pure-gas regions do not inherit a neighbour-averaged coupling.
    cappedInterfaceExchange = interfaceMask*cappedInterfaceExchange;
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

    surfaceExchange.correctBoundaryConditions();
    
    cappedInterfaceExchange = surfaceExchange/interfaceDelta;
    cappedInterfaceExchange.correctBoundaryConditions();


    tmp<volScalarField> tMaxInterfaceExchange;
    const volScalarField* maxInterfaceExchangePtr = nullptr;

    if (applyHeatFluxCap)
    {
        const dimensionedScalar unitTemperature
        (
            "unitTemperature",
            dimTemperature,
            1.0
        );

        tMaxInterfaceExchange = (qCapDim/unitTemperature)/interfaceDelta;
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

    const volScalarField& Told = T.oldTime();

    auto projectTemperature =
        [&](const word& stage)
        {
            if (!enableTClamp)
            {
                return;
            }

            const scalarField TintPre(T.internalField());
            label nLow = 0;
            label nHigh = 0;
            label nNonFinite = 0;
            forAll(TintPre, cellI)
            {
                const scalar val = TintPre[cellI];
                if (!std::isfinite(val))
                {
                    ++nNonFinite;
                }
                else if (val < minPhysicalTemp)
                {
                    ++nLow;
                }
                else if (val > maxPhysicalTemp)
                {
                    ++nHigh;
                }
            }

            label localCells = TintPre.size();
            Foam::reduce(nLow, sumOp<label>());
            Foam::reduce(nHigh, sumOp<label>());
            Foam::reduce(nNonFinite, sumOp<label>());
            Foam::reduce(localCells, sumOp<label>());

            totalLowProjected = scalar(nLow);
            totalHighProjected = scalar(nHigh);
            totalCellsProjected = scalar(localCells);

            tmp<volScalarField> tEnergy(CvolLimited*T);
            volScalarField& energy = tEnergy.ref();

            tmp<volScalarField> tNeighbourT(fvc::average(T));
            const volScalarField& neighbourT = tNeighbourT();

            tmp<volScalarField> tNeighbourCvol(fvc::average(CvolLimited));
            const volScalarField& neighbourCvol = tNeighbourCvol();

            const scalarField& capInternal = CvolLimited.internalField();
            const scalarField& neighbourTInternal = neighbourT.internalField();
            const scalarField& neighbourCapInternal = neighbourCvol.internalField();
            const scalarField& TlInternal = Tl.internalField();
            const scalarField& ToldInternal = Told.internalField();
            scalarField& energyInternal = energy.primitiveFieldRef();
            scalarField& Tint = T.primitiveFieldRef();

            label correctedInternal = 0;

            forAll(Tint, cellI)
            {
                scalar& Tcell = Tint[cellI];
                scalar& Ecell = energyInternal[cellI];
                const scalar cap = capInternal[cellI];

                const bool badT =
                    (!std::isfinite(Tcell)
                  || Tcell < minPhysicalTemp
                  || Tcell > maxPhysicalTemp);

                bool corrected = false;

                if (!std::isfinite(Ecell))
                {
                    corrected = true;
                }

                if (badT || !std::isfinite(Ecell))
                {
                    scalar candidateT = neighbourTInternal[cellI];
                    if (!std::isfinite(candidateT))
                    {
                        candidateT = ToldInternal[cellI];
                    }
                    if (!std::isfinite(candidateT))
                    {
                        candidateT = TlInternal[cellI];
                    }
                    if (!std::isfinite(candidateT))
                    {
                        candidateT = minPhysicalTemp;
                    }
                    candidateT = Foam::min
                    (
                        Foam::max(candidateT, minPhysicalTemp),
                        maxPhysicalTemp
                    );

                    scalar candidateEnergy = neighbourCapInternal[cellI]*candidateT;
                    if (!std::isfinite(candidateEnergy))
                    {
                        candidateEnergy = cap*candidateT;
                    }

                    const scalar energyFloor = cap*minPhysicalTemp;
                    const scalar energyCeil = cap*maxPhysicalTemp;

                    candidateEnergy = Foam::min
                    (
                        Foam::max(candidateEnergy, energyFloor),
                        energyCeil
                    );

                    Ecell = candidateEnergy;
                    Tcell = candidateEnergy/Foam::max(cap, VSMALL);
                    corrected = true;
                }
                else
                {
                    const scalar energyFloor = cap*minPhysicalTemp;
                    const scalar energyCeil = cap*maxPhysicalTemp;
                    const scalar limitedEnergy = Foam::min
                    (
                        Foam::max(Ecell, energyFloor),
                        energyCeil
                    );

                    if (limitedEnergy != Ecell)
                    {
                        Ecell = limitedEnergy;
                        corrected = true;
                    }

                    Tcell = Ecell/Foam::max(cap, VSMALL);
                }

                if (!std::isfinite(Tcell))
                {
                    scalar fallbackT = neighbourTInternal[cellI];
                    if (!std::isfinite(fallbackT))
                    {
                        fallbackT = ToldInternal[cellI];
                    }
                    if (!std::isfinite(fallbackT))
                    {
                        fallbackT = TlInternal[cellI];
                    }
                    if (!std::isfinite(fallbackT))
                    {
                        fallbackT = minPhysicalTemp;
                    }

                    fallbackT = Foam::min
                    (
                        Foam::max(fallbackT, minPhysicalTemp),
                        maxPhysicalTemp
                    );

                    Ecell = cap*fallbackT;
                    Tcell = fallbackT;
                    corrected = true;
                }

                if (corrected)
                {
                    ++correctedInternal;
                }
            }

            label correctedBoundary = 0;

            forAll(T.boundaryField(), patchI)
            {
                scalarField& Tp = T.boundaryFieldRef()[patchI];
                scalarField& Ep = energy.boundaryFieldRef()[patchI];
                const scalarField& CapP = CvolLimited.boundaryField()[patchI];
                const scalarField& TlP = Tl.boundaryField()[patchI];
                const scalarField& ToldP = Told.boundaryField()[patchI];

                forAll(Tp, faceI)
                {
                    scalar& Tface = Tp[faceI];
                    scalar& Eface = Ep[faceI];
                    const scalar cap = CapP[faceI];

                    const bool badT =
                        (!std::isfinite(Tface)
                      || Tface < minPhysicalTemp
                      || Tface > maxPhysicalTemp);

                    bool corrected = false;

                    if (badT || !std::isfinite(Eface))
                    {
                        scalar candidateT = Tface;
                        if (!std::isfinite(candidateT))
                        {
                            candidateT = ToldP[faceI];
                        }
                        if (!std::isfinite(candidateT))
                        {
                            candidateT = TlP[faceI];
                        }
                        if (!std::isfinite(candidateT))
                        {
                            candidateT = minPhysicalTemp;
                        }

                        candidateT = Foam::min
                        (
                            Foam::max(candidateT, minPhysicalTemp),
                            maxPhysicalTemp
                        );

                        Eface = cap*candidateT;
                        corrected = true;
                    }

                    const scalar energyFloor = cap*minPhysicalTemp;
                    const scalar energyCeil = cap*maxPhysicalTemp;
                    const scalar limitedEnergy = Foam::min
                    (
                        Foam::max(Eface, energyFloor),
                        energyCeil
                    );

                    if (limitedEnergy != Eface)
                    {
                        Eface = limitedEnergy;
                        corrected = true;
                    }

                    Tface = Eface/Foam::max(cap, VSMALL);

                    if (!std::isfinite(Tface))
                    {
                        scalar fallbackT = ToldP[faceI];
                        if (!std::isfinite(fallbackT))
                        {
                            fallbackT = TlP[faceI];
                        }
                        if (!std::isfinite(fallbackT))
                        {
                            fallbackT = minPhysicalTemp;
                        }

                        fallbackT = Foam::min
                        (
                            Foam::max(fallbackT, minPhysicalTemp),
                            maxPhysicalTemp
                        );

                        Tface = fallbackT;
                        Eface = cap*fallbackT;
                        corrected = true;
                    }

                    if (corrected)
                    {
                        ++correctedBoundary;
                    }
                }
            }

            T.correctBoundaryConditions();

            Foam::reduce(correctedInternal, sumOp<label>());
            Foam::reduce(correctedBoundary, sumOp<label>());

            const label correctedTotal = correctedInternal + correctedBoundary;

            if ((correctedTotal > 0 || nNonFinite > 0) && verbose && master)
            {
                WarningInFunction
                    << "Energy-consistent gas-temperature projection activated during "
                    << stage << ": corrected " << correctedTotal << " cells"
                    << " (low=" << nLow << ", high=" << nHigh
                    << ", non-finite=" << nNonFinite << ")" << endl;
            }
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

    if (fluxRelax < 1.0 - SMALL)
    {
        tmp<volScalarField> tDeltaTmag(mag(Tl - T));
        volScalarField& deltaTmag = tDeltaTmag.ref();

        const dimensionedScalar deltaTFloor
        (
            "deltaTFloor",
            dimTemperature,
            SMALL
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
    TEqn.source() += relaxedInterfaceExchange*Tl;

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

    projectTemperature("post-solve");

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
    projectTemperature("post-metal");

    tmp<volScalarField> tRawFlux(relaxedInterfaceExchange*(Tl - T));
    volScalarField& rawFlux = tRawFlux.ref();

    tmp<volScalarField> tNeighbourRawFlux(fvc::average(rawFlux));
    const volScalarField& neighbourRawFlux = tNeighbourRawFlux();

    rawFlux =
        interfaceMask*rawFlux
      + (scalar(1) - interfaceMask)*neighbourRawFlux;
    rawFlux.correctBoundaryConditions();

    const scalar dtVal = Foam::max(dt.value(), SMALL);

    scalarField& rawFluxInternal = rawFlux.primitiveFieldRef();
    const scalarField& capInternal = CvolLimited.internalField();
    const scalarField& TintInternal = T.internalField();

    forAll(rawFluxInternal, cellI)
    {
        const scalar cap = Foam::max(capInternal[cellI], VSMALL);
        const scalar Tcell = TintInternal[cellI];

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

    forAll(rawFlux.boundaryField(), patchI)
    {
        scalarField& fluxPatch = rawFlux.boundaryFieldRef()[patchI];
        const scalarField& capPatch = CvolLimited.boundaryField()[patchI];
        const scalarField& Tpatch = T.boundaryField()[patchI];

        forAll(fluxPatch, faceI)
        {
            const scalar cap = Foam::max(capPatch[faceI], VSMALL);
            const scalar Tface = Tpatch[faceI];

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

    tmp<volScalarField> tSurfaceFlux(rawFlux*interfaceDelta);
    volScalarField& surfaceFlux = tSurfaceFlux.ref();

    tmp<volScalarField> tNeighbourSurfaceFlux(fvc::average(surfaceFlux));
    const volScalarField& neighbourSurfaceFlux = tNeighbourSurfaceFlux();

    surfaceFlux =
        interfaceMask*surfaceFlux
      + (scalar(1) - interfaceMask)*neighbourSurfaceFlux;
    surfaceFlux.correctBoundaryConditions();

    if (qCap > 0)
    {
        surfaceFlux = Foam::min(Foam::max(surfaceFlux, -qCapDim), qCapDim);
        surfaceFlux.correctBoundaryConditions();
    }

    rawFlux = surfaceFlux/interfaceDelta;

    const volScalarField& metalToGasFlux = tRawFlux();

    gasMetalHeatFlux = metalToGasFlux;

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
            const scalar denom = Foam::max(scalar(1), totalCellsProjected);
            Info<< "    projection fractions: low=" << totalLowProjected/denom
                << "  high=" << totalHighProjected/denom << endl;
        }
    }

    mixture.correctThermo();
    mixture.correct();
}

// ************************************************************************* //
