/*---------------------------------------------------------------------------*\
  Temperature equation for LIFT simulation with the two-temperature model
  rewritten in a compact OpenFOAM style.
\*---------------------------------------------------------------------------*/

extern Foam::Switch verbose;

{
    if (verbose)
    {
        Info<< "Solving temperature equation" << endl;
    }

    const dictionary& ctrl = runTime.controlDict();
    const Switch enableTClamp = ctrl.lookupOrDefault<Switch>("enableTClamp", true);
    const scalar minPhysicalTemp = ctrl.lookupOrDefault<scalar>("minT", 250.0);
    const scalar maxPhysicalTemp = ctrl.lookupOrDefault<scalar>("maxT", 15000.0);

    if (enableTClamp)
    {
        const dimensionedScalar minTDim("minT", dimTemperature, minPhysicalTemp);
        const dimensionedScalar maxTDim("maxT", dimTemperature, maxPhysicalTemp);

        T = Foam::max(Foam::min(T, maxTDim), minTDim);
        T.correctBoundaryConditions();
    }

    volScalarField& gasMetalHeatFlux =
        mesh.lookupObjectRef<volScalarField>("gasMetalHeatFlux");

    tmp<volScalarField> tInterfaceExchange = ttm.gasMetalExchangeCoeffField();
    volScalarField& interfaceExchange = tInterfaceExchange.ref();

    const scalar minCoeffScalar =
        ctrl.lookupOrDefault<scalar>("minGasMetalExchangeCoeff", 0.0);
    const scalar maxCoeffScalar =
        ctrl.lookupOrDefault<scalar>("maxGasMetalExchangeCoeff", 1e18);

    const dimensionedScalar minCoeff
    (
        "minGasMetalExchangeCoeff",
        interfaceExchange.dimensions(),
        minCoeffScalar
    );
    const dimensionedScalar maxCoeff
    (
        "maxGasMetalExchangeCoeff",
        interfaceExchange.dimensions(),
        maxCoeffScalar
    );

    interfaceExchange = Foam::min(Foam::max(interfaceExchange, minCoeff), maxCoeff);
    interfaceExchange.correctBoundaryConditions();

    tmp<volScalarField> tCv2 = mixture.thermo2().Cv();
    tmp<volScalarField> tRho2 = mixture.thermo2().rho();

    const dimensionedScalar cvFloor("cvFloor", dimSpecificHeat, 1e-12);
    const dimensionedScalar rhoFloor("rhoFloor", dimDensity, 1e-12);

    tmp<volScalarField> tCvGas(Foam::max(tCv2(), cvFloor));
    tmp<volScalarField> tRhoGas(alpha2*Foam::max(tRho2(), rhoFloor));
    tmp<volScalarField> tCvolGas(tRhoGas()*tCvGas());

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
      - fvm::Sp(contErr, T)
      - fvm::laplacian(transportModel.kappaEff(), T)
     ==
        fvOptions(CvolGas, T)
    );

    const volScalarField& Tl = ttm.Tl();

    TEqn += fvm::Sp(interfaceExchange, T);
    TEqn.source() += interfaceExchange*Tl;
    // Laser absorption couples to the gas energy equation via the gas volume
    // fraction; metal absorption is handled by the two-temperature solver
    const volScalarField& alphaMetal = mixture.alpha1();
    const volScalarField& alphaGas = mixture.alpha2();
    TEqn.source() += alphaGas*mixture.Q_laser();
    static_cast<void>(alphaMetal);


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
        const dimensionedScalar minTDim("minT", dimTemperature, minPhysicalTemp);
        const dimensionedScalar maxTDim("maxT", dimTemperature, maxPhysicalTemp);
        T = Foam::max(Foam::min(T, maxTDim), minTDim);
    }

    T.correctBoundaryConditions();

    mixture.correctThermo();
    mixture.correct();

    const scalar fluxRelax =
        Foam::max(scalar(0), Foam::min(scalar(1), thermalFluxRelax));

    tmp<volScalarField> tRawFlux(interfaceExchange*(T - Tl));
    const volScalarField& rawFlux = tRawFlux();

    if (fluxRelax < 1.0 - SMALL)
    {
        gasMetalHeatFlux =
            fluxRelax*rawFlux + (1.0 - fluxRelax)*gasMetalHeatFlux;
    }
    else
    {
        gasMetalHeatFlux = rawFlux;
    }

    const scalar maxPhysicalFlux =
        ctrl.lookupOrDefault<scalar>("maxGasMetalHeatFlux", 1e15);

    if (maxPhysicalFlux > SMALL)
    {
        const dimensionedScalar fluxLim
        (
            "maxGasMetalHeatFlux",
            gasMetalHeatFlux.dimensions(),
            maxPhysicalFlux
        );
        const dimensionedScalar fluxLimNeg
        (
            "minGasMetalHeatFlux",
            gasMetalHeatFlux.dimensions(),
            -maxPhysicalFlux
        );

        gasMetalHeatFlux = Foam::min
        (
            Foam::max(gasMetalHeatFlux, fluxLimNeg),
            fluxLim
        );
    }

    gasMetalHeatFlux.correctBoundaryConditions();

    if (verbose)
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
