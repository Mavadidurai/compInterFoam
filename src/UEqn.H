fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  + transportModel.divDevRhoReff(U)
==
    fvOptions(rho, U)
);

const Foam::dictionary solverCoeffs = compInterFoamCoeffsDict(mesh);
const bool isMaster = Foam::Pstream::master();

//UEqn.relax();
fvOptions.constrain(UEqn);

volScalarField totalPressure(p_rgh);

// Add enhanced physics pressures
if (liftPhysics.valid())
{
    if (liftPhysics->phaseExplosionEnabled())
    {
        totalPressure += liftPhysics->explosivePressure();
    }

    if (liftPhysics->plasmaEnabled())
    {
        totalPressure += liftPhysics->plasmaPressure();
    }
}

// Read recoil pressure configuration consistently (same logic as pEqn.H)
// Priority: advancedInterfaceCapturing > compInterFoamCoeffs
const Foam::dictionary& controlDict = mesh.time().controlDict();
const Foam::dictionary& aicDict =
    controlDict.subOrEmptyDict("advancedInterfaceCapturing"); 

// The recoil is applied as an interface traction p_r*grad(alpha) below (via
// gradTotalP), like the surface tension. Its bound is maxRecoilPressure in
// advancedInterfaceCapturing. A second, optional clamp used to be computed
// and logged here, but its result was never used (dead code, removed).

// Under-relaxation of the momentum equation (URelaxationFactor), applied on
// all but the final PIMPLE iteration. Must be after the solverCoeffs declaration.
const scalar URelaxFactor = solverCoeffs.lookupOrDefault<scalar>("URelaxationFactor", 0.7);
if (URelaxFactor > 0 && URelaxFactor < 1.0 && !pimple.finalIter())
{
    // Keep the final PIMPLE sweep unrelaxed so the converged state is not
    // biased by the relaxation source term.
    UEqn.relax(URelaxFactor);
}

const Foam::scalar kineticEnergyCeiling =
    solverCoeffs.lookupOrDefault<Foam::scalar>
    (
        "maxKineticEnergyDensity",
        1e12
    );
const Foam::scalar minUEqnDiagVal =
    solverCoeffs.lookupOrDefault<Foam::scalar>
(
    "minUEqnDiag",
    1e-9
);

if (minUEqnDiagVal > 0)
{
    const tmp<volScalarField> tUEqnDiagTmp = UEqn.A();
    const volScalarField& UEqnDiagField = tUEqnDiagTmp();
    const dimensionedScalar diagFloor
    (
        "minUEqnDiag",
        UEqnDiagField.dimensions(),
        minUEqnDiagVal
    );

    const dimensionedScalar currentMinDiag = gMin(UEqnDiagField);
    if (currentMinDiag.value() < diagFloor.value())
    {
        UEqn += fvm::Sp(diagFloor, U);

        if (verbose && isMaster)
        {
            Info<< "UEqn: added stabilizing diagonal " << diagFloor.value()
                << " to avoid singular momentum matrix" << endl;
        }
    }
}

if (kineticEnergyCeiling > 0 && kineticEnergyCeiling < Foam::GREAT)
{
    // 0.5*rho*|U|^2 has SI units of J/m^3; it is the kinetic energy stored in
    // a unit volume of the flow.  Checking it before the predictor prevents
    // unphysical velocities from entering the momentum solve.
    tmp<volScalarField> tKineticEnergyDensity(0.5*rho*magSqr(U));
    const volScalarField& kineticEnergyDensity = tKineticEnergyDensity();
    const dimensionedScalar maxKineticEnergyDensity =
        gMax(kineticEnergyDensity);
    const scalar maxKineticEnergyDensityValue =
        maxKineticEnergyDensity.value();

    if
    (
        !std::isfinite(maxKineticEnergyDensityValue)
     || maxKineticEnergyDensityValue > kineticEnergyCeiling
    )
    {
        FatalErrorInFunction
            << "Maximum kinetic energy density "
            << maxKineticEnergyDensityValue
            << " J/m^3 exceeds limit "
            << kineticEnergyCeiling
            << " J/m^3. Velocity field has become unphysical."
            << exit(FatalError);
    }
}

if (pimple.momentumPredictor())
{
    // Surface tension and gravity enter through phig in pEqn.H. The stock
    // compressibleInterFoam also puts them on the right-hand side of the
    // momentum predictor; that is not a double count, because HbyA =
    // rAU*UEqn.H() does not contain the predictor's right-hand side. Leaving
    // them out of the predictor here is a choice that only affects the
    // predicted velocity, not the converged PIMPLE solution, which the
    // pressure correctors determine.
    tmp<volVectorField> tGradTotalP = fvc::grad(totalPressure);
    volVectorField gradTotalP
    (
        IOobject
        (
            "grad(totalPressure)",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tGradTotalP()
    );
    tGradTotalP.clear();

    // Fix: recoil traction F = p_r*grad(alpha1), pushing into the metal.
    // Face form, identical to the phig term in pEqn.H, so predictor and
    // corrector see the same force. (Not fvc::grad(alpha1): with
    // cellLimited it vanishes at a sharp interface.)
    if (recoilPressurePtr)
    {
        gradTotalP -= fvc::reconstruct
        (
            recoilFaceValue(*recoilPressurePtr)*fvc::snGrad(alpha1)
           *mesh.magSf()
        );
    }

    const Foam::scalar gradLimit = solverCoeffs.lookupOrDefault<Foam::scalar>
    (
        "maxPressureGradient",
        Foam::GREAT
    );

    bool nonFiniteGrad = false;
    bool gradClamped = false;

    vectorField& gradInternal = gradTotalP.internalFieldRef();

    forAll(gradInternal, celli)
    {
        vector& g = gradInternal[celli];

        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            scalar& val = g[cmpt];

            if (!std::isfinite(val))
            {
                val = 0.0;
                nonFiniteGrad = true;
            }
        }

        if (gradLimit < Foam::GREAT)
        {
            const scalar magG = Foam::mag(g);

            if (magG > gradLimit && magG > SMALL)
            {
                g *= gradLimit/magG;
                gradClamped = true;
            }
        }
    }

    forAll(gradTotalP.boundaryField(), patchi)
    {
        vectorField& gradPatch = gradTotalP.boundaryFieldRef()[patchi];

        forAll(gradPatch, facei)
        {
            vector& g = gradPatch[facei];

            for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
            {
                scalar& val = g[cmpt];

                if (!std::isfinite(val))
                {
                    val = 0.0;
                    nonFiniteGrad = true;
                }
            }

            if (gradLimit < Foam::GREAT)
            {
                const scalar magG = Foam::mag(g);

                if (magG > gradLimit && magG > SMALL)
                {
                    g *= gradLimit/magG;
                    gradClamped = true;
                }
            }
        }
    }

    if (nonFiniteGrad)
    {
        WarningInFunction
            << "Non-finite entries detected in grad(totalPressure); resetting to zero"
            << nl;
    }

    if (gradClamped)
    {
        Info<< "    grad(totalPressure) limited to max magnitude "
            << gradLimit << " Pa/m" << endl;
    }

    tmp<volVectorField> tPressureForce(-gradTotalP);
    volVectorField& pressureForce = tPressureForce.ref();

    tmp<volVectorField> tMomentumSource(tPressureForce);

    const scalar domainVolume = gSum(mesh.V());
    const scalar deltaT = runTime.deltaTValue();
    const bool logVelocityDiagnostics =
        verbose && isMaster && (runTime.timeIndex() % 10 == 0);

    if (logVelocityDiagnostics && domainVolume > SMALL)
    {
        tmp<volVectorField> tPressureAccel(pressureForce/rho);
        const volVectorField& pressureAccel = tPressureAccel();
        tmp<volScalarField> tPressureAccelMag = mag(pressureAccel);
        const volScalarField& pressureAccelMag = tPressureAccelMag();

        const dimensionedScalar pressureAccelIntegral =
            fvc::domainIntegrate(pressureAccelMag);
        const scalar maxPressureAccel = gMax(pressureAccelMag);
        const scalar meanPressureAccel =
            pressureAccelIntegral.value()/domainVolume;

        Info<< "    Velocity diagnostics:" << nl
            << "      Max total-pressure acceleration: "
            << maxPressureAccel << " m/s^2 (dU ~ "
            << maxPressureAccel*deltaT << " m/s per step)" << nl
            << "      Mean total-pressure acceleration: "
            << meanPressureAccel << " m/s^2" << endl;
    }

    volVectorField momentumSource
    (
        IOobject
        (
            "momentumSource",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tMomentumSource()
    );
    tPressureForce.clear();
    tMomentumSource.clear();

    // Solve with implicit total-pressure forcing
    solve(UEqn == momentumSource);
    fvOptions.correct(U);
    U.correctBoundaryConditions();

    K = 0.5*magSqr(U);
}

tmp<volScalarField> tUEqnDiag = UEqn.A();
volScalarField& UEqnDiag = tUEqnDiag.ref();

const Foam::scalar maxUEqnVelocity = solverCoeffs.lookupOrDefault<Foam::scalar>
(
    "maxUEqnVelocity",
    Foam::GREAT
);

if (minUEqnDiagVal > 0)
{
    UEqnDiag = Foam::max
    (
        UEqnDiag,
        dimensionedScalar
        (
            "minUEqnDiag",
            UEqnDiag.dimensions(),
            minUEqnDiagVal
        )
    );
}

volScalarField rAU("rAU", 1.0/tUEqnDiag());

if (maxUEqnVelocity < Foam::GREAT)
{
    bool velocityLimited = false;

    vectorField& UInternal = U.primitiveFieldRef();

    forAll(UInternal, celli)
    {
        vector& Ui = UInternal[celli];
        scalar magUi = Foam::mag(Ui);

        if (!std::isfinite(magUi))
        {
            Ui = vector::zero;
            velocityLimited = true;
            continue;
        }

        if (magUi > maxUEqnVelocity && magUi > SMALL)
        {
            Ui *= maxUEqnVelocity/magUi;
            velocityLimited = true;
        }
    }

    forAll(U.boundaryField(), patchi)
    {
        vectorField& Up = U.boundaryFieldRef()[patchi];

        forAll(Up, facei)
        {
            vector& Ui = Up[facei];
            scalar magUi = Foam::mag(Ui);

            if (!std::isfinite(magUi))
            {
                Ui = vector::zero;
                velocityLimited = true;
                continue;
            }

            if (magUi > maxUEqnVelocity && magUi > SMALL)
            {
                Ui *= maxUEqnVelocity/magUi;
                velocityLimited = true;
            }
        }
    }

    if (velocityLimited)
    {
        Info<< "    Velocity field limited to max magnitude "
            << maxUEqnVelocity << " m/s before computing flux" << endl;
    }
}

// Optional limit on the gas-phase speed (compInterFoamCoeffs.maxVaporVelocity).
// It rescales U in every cell with alpha.metal < 0.1, which includes the plain
// gas, and it does not conserve momentum. It is therefore OFF by default (it
// used to default to 2000 m/s without being set by any case), and when it is
// used it acts BEFORE the flux is built, so that phi and U stay consistent.
const scalar maxVaporVelocity =
    solverCoeffs.lookupOrDefault<scalar>("maxVaporVelocity", Foam::GREAT);

if (maxVaporVelocity < Foam::GREAT)
{
    bool vaporVelocityLimited = false;
    const scalarField& alpha1Field = alpha1.primitiveField();
    vectorField& UInternal = U.primitiveFieldRef();

    forAll(UInternal, cellI)
    {
        // Only limit in vapor regions (low alpha.metal)
        if (alpha1Field[cellI] < 0.1)  // Vapor/air phase
        {
            vector& Ui = UInternal[cellI];
            scalar magUi = Foam::mag(Ui);

            if (magUi > maxVaporVelocity && magUi > SMALL)
            {
                Ui *= maxVaporVelocity / magUi;  // Scale down to limit
                vaporVelocityLimited = true;
            }
        }
    }


    // Also limit boundary fields
    forAll(U.boundaryField(), patchi)
    {
        const scalarField& alpha1Patch = alpha1.boundaryField()[patchi];
        vectorField& Up = U.boundaryFieldRef()[patchi];

        forAll(Up, facei)
        {
            if (alpha1Patch[facei] < 0.1)  // Vapor/air phase
            {
                vector& Ui = Up[facei];
                scalar magUi = Foam::mag(Ui);

                if (magUi > maxVaporVelocity && magUi > SMALL)
                {
                    Ui *= maxVaporVelocity / magUi;
                    vaporVelocityLimited = true;
                }
            }
        }
    }

    if (vaporVelocityLimited && verbose && isMaster)
    {
        Info<< "    Vapor-phase velocity limited to " << maxVaporVelocity
            << " m/s (compInterFoamCoeffs.maxVaporVelocity)" << endl;
    }
}

U.correctBoundaryConditions();
fvOptions.correct(U);

phi = fvc::flux(U);
rhoPhi = fvc::interpolate(rho)*phi;
fvc::makeRelative(rhoPhi, U);
