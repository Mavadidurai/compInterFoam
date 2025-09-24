/*---------------------------------------------------------------------------*\
  Temperature equation for LIFT simulation with the two-temperature model
  rewritten in a compact OpenFOAM style.
\*---------------------------------------------------------------------------*/

extern Foam::Switch verbose;

{
    if (verbose && master)
    {
        Info<< "Solving temperature equation" << endl;
    }

    const dictionary& ctrl = runTime.controlDict();
    const Switch enableTClamp = ctrl.lookupOrDefault<Switch>("enableTClamp", false);
    const scalar minPhysicalTemp = ctrl.lookupOrDefault<scalar>("minT", 0.0);
    const scalar maxPhysicalTemp = ctrl.lookupOrDefault<scalar>
    (
        "maxT",
        Foam::GREAT
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
    auto clampTemperature =
        [&](const word& stage)
        {
            const scalar oldMin = gMin(T.internalField());
            const scalar oldMax = gMax(T.internalField());

            const dimensionedScalar minTDim("minT", dimTemperature, minPhysicalTemp);
            const dimensionedScalar maxTDim("maxT", dimTemperature, maxPhysicalTemp);

            T = Foam::max(Foam::min(T, maxTDim), minTDim);
            T.correctBoundaryConditions();

            if
            (
                (oldMin < minPhysicalTemp - SMALL)
             || (oldMax > maxPhysicalTemp + SMALL)
            )
            {
                if (Pstream::master())
                {
                    WarningInFunction
                        << "Temperature clamp activated during " << stage
                        << ": enforcing [" << minPhysicalTemp << ", "
                        << maxPhysicalTemp << "] K" << endl;
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

    const dictionary twoTemperatureDict =
        ctrl.subOrEmptyDict("twoTemperatureProperties");

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
    interfaceExchange = Foam::min(Foam::max(interfaceExchange, minCoeff), maxCoeff);
    interfaceExchange.correctBoundaryConditions();

    tmp<volScalarField> tCv2 = mixture.thermo2().Cv();
    tmp<volScalarField> tRho2 = mixture.thermo2().rho();

   const dimensionedScalar cvFloor
    (
        "cvFloor",
        dimEnergy/dimMass/dimTemperature,
        1e-12
    );
    const dimensionedScalar rhoFloor("rhoFloor", dimDensity, 1e-12);

    tmp<volScalarField> tCvGas(Foam::max(tCv2(), cvFloor));
    const volScalarField& CvGas = tCvGas();

    tmp<volScalarField> tRhoGasCore(Foam::max(tRho2(), rhoFloor));
    const volScalarField& rhoGasCore = tRhoGasCore();

    tmp<volScalarField> tCvolGas(alpha2*rhoGasCore*CvGas);

    const volScalarField& CvolGas = tCvolGas();

    tmp<surfaceScalarField> tCvolPhi(fvc::interpolate(CvolGas)*phi);
    const surfaceScalarField& CvolPhi = tCvolPhi();
    tmp<volScalarField> tDivPhi(fvc::div(phi));
    const volScalarField& divPhi = tDivPhi();
    tmp<volScalarField> tContErr(CvolGas*divPhi);
    const volScalarField& contErr = tContErr();

    fvScalarMatrix TEqn
    (
        fvm::ddt(CvolGas, T)
      + fvm::div(CvolPhi, T)
      + fvm::Sp(contErr, T)
      - fvm::laplacian(transportModel.kappaEff(), T)
     ==
        fvOptions(CvolGas, T)
    );

    const volScalarField& Tl = ttm.Tl();

    TEqn += fvm::Sp(interfaceExchange, T);
    TEqn.source() += interfaceExchange*Tl;

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
        T.correctBoundaryConditions();
    }

    mixture.correctThermo();
    mixture.correct();

    const scalar fluxRelax =
        Foam::max(scalar(0), Foam::min(scalar(1), thermalFluxRelax));

    tmp<volScalarField> tRawFlux(interfaceExchange*(T - Tl));
    volScalarField& rawFlux = tRawFlux.ref();

    const scalar maxThermalCourant =
        ctrl.lookupOrDefault<scalar>("maxThermalCourant", 1.0);

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

        tmp<volScalarField> tCvolLimited(rhoGasCore*CvGas);
        const volScalarField& CvolLimited = tCvolLimited();

        tmp<volScalarField> tThermalCo
        (
            (interfaceExchange*runTime.deltaT())/CvolLimited
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
        tmp<volScalarField> tLimiterActive(Foam::pos(alpha2 - alphaLimiterFloor));
        const volScalarField& limiterActive = tLimiterActive();

        limiter = limiterActive*limiter + (oneDim - limiterActive)*oneDim;
        rawFlux *= limiter;
    }
    {
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

        volScalarField& activeMask = tActiveMask.ref();
        activeMask = Foam::min(Foam::max(activeMask, scalar(0)), scalar(1));

        tmp<volScalarField> tMetalClamp
        (
            Foam::min(Foam::max(metalFraction, scalar(0)), scalar(1))
        );
        const volScalarField& metalClamp = tMetalClamp();

        tmp<volScalarField> tContactMask
        (
            new volScalarField(fvc::average(metalFraction))
        );
        volScalarField& contactMask = tContactMask.ref();
        contactMask = Foam::min(Foam::max(contactMask, scalar(0)), scalar(1));
        contactMask = Foam::max(contactMask, metalClamp);

        activeMask *= contactMask;
        activeMask.correctBoundaryConditions();

        rawFlux *= activeMask;
    }


    const volScalarField& limitedFlux = tRawFlux();

    if (fluxRelax < 1.0 - SMALL)
    {
        gasMetalHeatFlux =
            fluxRelax*limitedFlux + (1.0 - fluxRelax)*gasMetalHeatFlux;
    }
    else
    {
        gasMetalHeatFlux = limitedFlux;
    }

    gasMetalHeatFlux.correctBoundaryConditions();

    if (verbose && master)
    {
        Info<< "    min(T) = " << min(T).value()
            << "  max(T) = " << max(T).value()
            << "  mean(T) = " << gAverage(T) << endl;
        Info<< "    mean exchange coeff = "
            << gAverage(interfaceExchange)
            << "  max |q_gm| = "
            << gMax(mag(gasMetalHeatFlux)().internalField()) << endl;
    }
}

// ************************************************************************* //
