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

bool recoilClampActive = true;
scalar recoilClampMax = 5e10;  // Default: 50 GPa (realistic for fs-LIFT)

if (!aicDict.empty())
{

    recoilClampActive = aicDict.lookupOrDefault<Foam::Switch>("clampRecoil", true);
    // Try maxRecoilPressure first, then recoilMax (legacy)
    scalar maxRecoil = aicDict.lookupOrDefault<scalar>("maxRecoilPressure", 0.0);
    if (maxRecoil <= SMALL)
    {
        maxRecoil = aicDict.lookupOrDefault<scalar>("recoilMax", recoilClampMax);
    }
    if (std::isfinite(maxRecoil) && maxRecoil > SMALL)
    {
        recoilClampMax = Foam::mag(maxRecoil);
    }
}
// Fallback to solver coeffs if not in aicDict
if (recoilClampMax <= SMALL || recoilClampMax >= 0.99*Foam::GREAT)
{
    recoilClampMax = solverCoeffs.lookupOrDefault<scalar>("maxRecoilPressure", 1e8);
}

const bool clampRecoilContribution =
    recoilClampActive
 && recoilClampMax > SMALL
 && recoilClampMax < 0.99*Foam::GREAT;

if (recoilPressurePtr)
{
    const volScalarField& recoilP = *recoilPressurePtr;
    const volScalarField* recoilContribution = &recoilP;
    tmp<volScalarField> tRecoilClamped;

    if (clampRecoilContribution)
    {
        const dimensionedScalar maxRecoil
        (
            "maxRecoilPressure",
            dimPressure,
            recoilClampMax
        );
        const dimensionedScalar minRecoil
        (
            "minRecoilPressure",
            dimPressure,
            -recoilClampMax
        );

        const scalar observedMaxRecoil = gMax(recoilP);
        const scalar observedMinRecoil = gMin(recoilP);

        tRecoilClamped = Foam::min
        (
            Foam::max(recoilP, minRecoil),
            maxRecoil
        );
        recoilContribution = &tRecoilClamped();

        if
        (
            verbose
         && isMaster
         && (observedMaxRecoil > recoilClampMax
          || observedMinRecoil < -recoilClampMax)
        )
        {
            Info<< "    recoilPressureClamp: limiting recoil contribution to ±"
                << recoilClampMax/1e6 << " MPa" << endl;
        }
    }

    totalPressure += *recoilContribution;
}

// Less aggressive relaxation to allow recoil-driven acceleration
// For explicit coupling with high recoil pressures (38 GPa), 30% damping is too aggressive
const scalar URelaxFactor = solverCoeffs.lookupOrDefault<scalar>("URelaxationFactor", 0.7);
if (URelaxFactor > 0 && URelaxFactor < 1.0)
{
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
    // The capillary force already appears in phig within pEqn.H via
    // mixture.surfaceTensionForce().  Including it here as well was
    // double-counting the surface-tension contribution, which the run
    // logs revealed as steadily increasing velocities even when the
    // recoil pressure field remained identically zero.  Only assemble
    // the forces that are not already handled through the pressure
    // correction stage.

    // Start with the forces that are *not* treated in phig.  Gravity is
    // already accounted for through the -ghf*fvc::snGrad(rho) contribution
    // inside phig in pEqn.H, so including rho*g here would double-count it
    // and steadily accelerate the flow even in the absence of recoil.
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

// Limit vapor velocities to physical sound speed
// Sound speed in hot metal vapor is ~1000-2000 m/s
// Supersonic expansion without shock physics leads to numerical instability
const scalar maxVaporVelocity = solverCoeffs.lookupOrDefault<scalar>("maxVaporVelocity", 2000.0);
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
            << " m/s (sound speed in hot metal vapor)" << endl;
    }
}

phi = fvc::flux(U);
rhoPhi = fvc::interpolate(rho)*phi;
fvc::makeRelative(rhoPhi, U);

U.correctBoundaryConditions();
fvOptions.correct(U);
