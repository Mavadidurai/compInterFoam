fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  + transportModel.divDevRhoReff(U)
==
    fvOptions(rho, U)
);

const Foam::dictionary solverCoeffs = compInterFoamCoeffsDict(mesh);

//UEqn.relax();
fvOptions.constrain(UEqn);

volScalarField totalPressure(p_rgh);
//totalPressure += pInterfaceCapturing->recoilPressure();

// Add enhanced physics pressures
if (liftPhysics.valid())
{
    if (liftPhysics->phaseExplosionEnabled())
        totalPressure += liftPhysics->explosivePressure();
    
    if (liftPhysics->plasmaEnabled())
        totalPressure += liftPhysics->plasmaPressure();
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

// CRITICAL: Aggressive under-relaxation for extreme pressure gradients
// Must be AFTER solverCoeffs declaration. Reference: Tryggvason et al. (2011)
const scalar URelaxFactor = solverCoeffs.lookupOrDefault<scalar>("URelaxationFactor", 0.3);
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

const bool isMaster = Foam::Pstream::master();
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

// Recoil force calculation removed - recoil pressure now applied implicitly
// through p_rgh field to maintain proper PIMPLE pressure-velocity coupling.
// Explicit recoil forces caused triple-counting and broke momentum conservation.
//
// Reference: Brackbill et al. (1992) - VOF method couples pressure implicitly
// Reference: Knight (1979) - Recoil pressure physics (applied via p_rgh now)

tmp<volVectorField> tRecoilForce;  // Keep declaration for compatibility
bool recoilForceReady = false;     // Always false now

recoilTractionf = dimensionedVector
(
    "zeroRecoilTraction",
    dimForce/dimVolume,
    vector::zero
);

// ========================================================================
// REMOVED: Explicit recoil force calculation (formerly ~250 lines)
// ========================================================================
// Recoil pressure is now applied ONLY through the p_rgh pressure field
// in advancedInterfaceCapturing.C (Knight 1979 formula).
//
// WHY REMOVED:
// 1. Triple-counting bug: Recoil pressure was added to p_rgh, then converted
//    to explicit force here, then re-applied in pEqn.H via HbyA and phig
// 2. Breaks PIMPLE coupling: Explicit forces bypass the pressure-velocity
//    coupling that maintains momentum conservation in the PIMPLE algorithm
// 3. Unphysical results: Caused ±40 GPa pressure oscillations instead of
//    the physically realistic ~500 MPa recoil pressures
//
// PHYSICS REFERENCES:
// - Knight (1979): Recoil pressure formula (implemented in advancedInterfaceCapturing.C)
// - Brackbill et al. (1992): VOF pressure forces should be implicit, not explicit
// - Feinaeugle et al. (2017): Experimental validation showing ~100-500 MPa recoil
//
// The recoil pressure field is still calculated and added to p_rgh in
// advancedInterfaceCapturing.C. This implicit treatment through the pressure
// field maintains proper momentum conservation and PIMPLE algorithm stability.
// ========================================================================

// Add phase explosion and plasma pressure forces
tmp<volVectorField> tExplosiveForce;
tmp<volVectorField> tPlasmaForce;
bool explosiveForceReady = false;
bool plasmaForceReady = false;

if (recoilPressurePtr)  // Only if interface capturing is active
{
    const tmp<volVectorField> tInterfaceGrad = fvc::grad(alpha1);
    const volVectorField& interfaceGrad = tInterfaceGrad();
    const tmp<volScalarField> tInterfaceStrength = mag(interfaceGrad);
    const volScalarField& interfaceStrength = tInterfaceStrength();
    const tmp<volVectorField> tInterfaceNormal
    (
        -(interfaceGrad/
        (
            interfaceStrength
          + dimensionedScalar("small", dimless/dimLength, SMALL)
        ))
    );
    const volVectorField& interfaceNormal = tInterfaceNormal();

    // Add phase explosion pressure
    if (mesh.foundObject<volScalarField>("explosivePressure"))
    {
        const volScalarField& explosiveP =
            mesh.lookupObject<volScalarField>("explosivePressure");

        tExplosiveForce = tmp<volVectorField>
        (
            explosiveP*interfaceStrength*interfaceNormal
        );
        explosiveForceReady = true;

        if (verbose && isMaster && runTime.timeIndex() % 10 == 0)
        {
            const scalar maxExplosiveP = gMax(explosiveP);
            Info<< "    Max explosive pressure: " << maxExplosiveP/1e6
                << " MPa" << endl;
        }
    }

    // Add plasma pressure
    if (mesh.foundObject<volScalarField>("plasmaPressure"))
    {
        const volScalarField& plasmaP =
            mesh.lookupObject<volScalarField>("plasmaPressure");

        tPlasmaForce = tmp<volVectorField>
        (
            plasmaP*interfaceStrength*interfaceNormal
        );
        plasmaForceReady = true;

        if (verbose && isMaster && runTime.timeIndex() % 10 == 0)
        {
            const scalar maxPlasmaP = gMax(plasmaP);
            Info<< "    Max plasma pressure: " << maxPlasmaP/1e6
                << " MPa" << endl;
        }
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
    tmp<volVectorField> tGradPrgh = fvc::grad(p_rgh);
    volVectorField gradPrgh
    (
        IOobject
        (
            "grad(p_rgh)",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tGradPrgh()
    );
    tGradPrgh.clear();

    const Foam::scalar gradLimit = solverCoeffs.lookupOrDefault<Foam::scalar>
    (
        "maxPressureGradient",
        Foam::GREAT
    );

    bool nonFiniteGrad = false;
    bool gradClamped = false;

    vectorField& gradInternal = gradPrgh.internalFieldRef();

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

    forAll(gradPrgh.boundaryField(), patchi)
    {
        vectorField& gradPatch = gradPrgh.boundaryFieldRef()[patchi];

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
            << "Non-finite entries detected in grad(p_rgh); resetting to zero"
            << nl;
    }

    if (gradClamped)
    {
        Info<< "    grad(p_rgh) limited to max magnitude "
            << gradLimit << " Pa/m" << endl;
    }
    
    tmp<volVectorField> tPressureForce(-gradPrgh);
    volVectorField& pressureForce = tPressureForce.ref();

    tmp<volVectorField> tMomentumSource(tPressureForce);

    // Recoil force removed - was lines 353-358
    // recoilForceFieldPtr no longer needed since recoil handled via p_rgh

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

        // Recoil acceleration calculation removed - was lines 378-392
        // Recoil now in p_rgh, not as explicit force

        tmp<volVectorField> tTotalAccel(tMomentumSource()/rho);
        volVectorField& totalAccel = tTotalAccel.ref();

        // Recoil force removed from totalAccel - was lines 380-383
        // Recoil pressure now applied implicitly via p_rgh only

        if (explosiveForceReady)
        {
            totalAccel += tExplosiveForce()/rho;
        }

        if (plasmaForceReady)
        {
            totalAccel += tPlasmaForce()/rho;
        }
        tmp<volScalarField> tTotalAccelMag = mag(totalAccel);
        const volScalarField& totalAccelMag = tTotalAccelMag();

        const dimensionedScalar totalAccelIntegral =
            fvc::domainIntegrate(totalAccelMag);
        const scalar maxTotalAccel = gMax(totalAccelMag);
        const scalar meanTotalAccel =
            totalAccelIntegral.value()/domainVolume;

        Info<< "    Velocity diagnostics:" << nl
            << "      Max pressure acceleration: "
            << maxPressureAccel << " m/s^2 (dU ~ "
            << maxPressureAccel*deltaT << " m/s per step)" << nl
            << "      Max total acceleration:    "
            << maxTotalAccel << " m/s^2 (dU ~ "
            << maxTotalAccel*deltaT << " m/s per step)" << nl
            << "      Mean total acceleration:   "
            << meanTotalAccel << " m/s^2" << nl
            << "      Mean pressure acceleration: "
            << meanPressureAccel << " m/s^2" << endl;
        // Recoil acceleration output removed - recoil now in p_rgh field
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

    // ========== ADD ENHANCED PHYSICS FORCES ==========
    // Add phase explosion force
    if (explosiveForceReady)
    {
        momentumSource += tExplosiveForce();
    }

    // Add plasma force
    if (plasmaForceReady)
    {
        momentumSource += tPlasmaForce();
    }
    // ========== END ENHANCED PHYSICS ==========

    // Solve with complete momentum source
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

phi = fvc::flux(U);
rhoPhi = fvc::interpolate(rho)*phi;
fvc::makeRelative(rhoPhi, U);

U.correctBoundaryConditions();
fvOptions.correct(U);
