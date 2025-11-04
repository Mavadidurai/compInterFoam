fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  + transportModel.divDevRhoReff(U)
 ==
    fvOptions(rho, U)
);

UEqn.relax();
fvOptions.constrain(UEqn);

const Foam::dictionary& solverCoeffs = compInterFoamCoeffsDict(mesh);
// Read recoil pressure configuration consistently (same logic as pEqn.H)
// Priority: advancedInterfaceCapturing > compInterFoamCoeffs
const Foam::dictionary& controlDict = mesh.time().controlDict();
const Foam::dictionary& aicDict =
    controlDict.subOrEmptyDict("advancedInterfaceCapturing"); 

bool recoilClampActive = true;
scalar recoilClampMax = 1e8;  // Default: 100 MPa

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

if (!recoilClampActive)
{
    recoilClampMax = Foam::GREAT;
}

const Foam::scalar kineticEnergyCeiling =
    solverCoeffs.lookupOrDefault<Foam::scalar>
    (
        "maxKineticEnergyDensity",
        1e12
    );

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

tmp<volVectorField> tRecoilForce;
bool recoilForceReady = false;

recoilTractionf = dimensionedVector
(
    "zeroRecoilTraction",
    dimForce/dimVolume,
    vector::zero
);

if (recoilPressurePtr)
{
        const volScalarField& recoilP = *recoilPressurePtr;

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
        const bool logRecoil = verbose && master && runTime.timeIndex() % 10 == 0;

        // Validate recoil pressure before using
        const scalar maxRecoil = gMax(recoilP);
        const scalar minRecoil = gMin(recoilP);

        scalar lowerBound = -Foam::GREAT;
        scalar upperBound = Foam::GREAT;
        bool lowerConfigured = false;
        bool upperConfigured = false;

        if (solverCoeffs.found("minRecoilPressure"))
        {
            const scalar value =
                solverCoeffs.lookupOrDefault<scalar>("minRecoilPressure", lowerBound);

            if (std::isfinite(value))
            {
                lowerBound = value;
                lowerConfigured = true;
            }
        }

        if (solverCoeffs.found("maxRecoilPressure"))
        {
            const scalar value =
                solverCoeffs.lookupOrDefault<scalar>("maxRecoilPressure", upperBound);

            if (std::isfinite(value))
            {
                upperBound = value;
                upperConfigured = true;
            }
        }

        if (solverCoeffs.found("recoilMax"))
        {
            const scalar limit =
                Foam::mag(solverCoeffs.lookupOrDefault<scalar>("recoilMax", 0.0));

            if (std::isfinite(limit) && limit > SMALL)
            {
                lowerBound = -limit;
                upperBound = limit;
                lowerConfigured = true;
                upperConfigured = true;
            }
        }

        const dictionary& aicDict =
            mesh.time().controlDict().subOrEmptyDict("advancedInterfaceCapturing");

        if (!aicDict.empty())
        {
            if (aicDict.found("minRecoilPressure"))
            {
                const scalar value =
                    aicDict.lookupOrDefault<scalar>("minRecoilPressure", lowerBound);

                if (std::isfinite(value))
                {
                    lowerBound = value;
                    lowerConfigured = true;
                }
            }

            if (aicDict.found("maxRecoilPressure"))
            {
                const scalar value =
                    aicDict.lookupOrDefault<scalar>("maxRecoilPressure", upperBound);

                if (std::isfinite(value))
                {
                    upperBound = value;
                    upperConfigured = true;
                }
            }

            const Foam::Switch clampRecoil =
                aicDict.lookupOrDefault<Foam::Switch>("clampRecoil", Foam::Switch(false));

            if (clampRecoil)
            {
                const scalar limit =
                    Foam::mag(aicDict.lookupOrDefault<scalar>("recoilMax", 0.0));

                if (std::isfinite(limit) && limit > SMALL)
                {
                    lowerBound = -limit;
                    upperBound = limit;
                    lowerConfigured = true;
                    upperConfigured = true;
                }
            }
        }

        const bool invalidRecoil =
            (!std::isfinite(maxRecoil) || !std::isfinite(minRecoil));

        if (invalidRecoil && !(lowerConfigured || upperConfigured))
        {
            lowerBound = 0.0;
            upperBound = 0.0;
            lowerConfigured = true;
            upperConfigured = true;
        }

        const bool exceedsUpper = upperConfigured && maxRecoil > upperBound;
        const bool exceedsLower = lowerConfigured && minRecoil < lowerBound;
        const bool clampRequired = invalidRecoil || exceedsUpper || exceedsLower;

        tmp<volScalarField> tRecoilClamped;
        const volScalarField* recoilToUse = &recoilP;

        if (clampRequired)
        {
            Ostream& warn = WarningInFunction;
            warn << "Recoil pressure out of physical range: ["
                << minRecoil << ", " << maxRecoil << "] Pa" << nl;

            if (lowerConfigured || upperConfigured)
            {
                warn << "Clamping to [";

                if (lowerConfigured)
                {
                    warn << lowerBound;
                }
                else
                {
                    warn << "-inf";
                }

                warn << ", ";

                if (upperConfigured)
                {
                    warn << upperBound;
                }
                else
                {
                    warn << "inf";
                }

                warn << "] Pa" << endl;
            }
            else
            {
                warn << "No recoil bounds configured; clamping to zero." << endl;
            }

            if (lowerConfigured && upperConfigured)
            {
                tRecoilClamped = Foam::min
                (
                    Foam::max
                    (
                        recoilP,
                        dimensionedScalar("minRecoilPressure", dimPressure, lowerBound)
                    ),
                    dimensionedScalar("maxRecoilPressure", dimPressure, upperBound)
                );
            }
            else if (lowerConfigured)
            {
                tRecoilClamped = Foam::max
                (
                    recoilP,
                    dimensionedScalar("minRecoilPressure", dimPressure, lowerBound)
                );
            }
            else if (upperConfigured)
            {
                tRecoilClamped = Foam::min
                (
                    recoilP,
                    dimensionedScalar("maxRecoilPressure", dimPressure, upperBound)
                );
            }
            else
            {
                tRecoilClamped = recoilP;
                tRecoilClamped.ref() = dimensionedScalar("zero", dimPressure, 0.0);
            }

            recoilToUse = &tRecoilClamped();
        }

        tRecoilForce = tmp<volVectorField>
        (
            (*recoilToUse)*interfaceStrength*interfaceNormal
        );
        volVectorField& recoilForce = tRecoilForce.ref();       
        recoilTractionf = fvc::interpolate(recoilForce);
        recoilForceReady = true;

        if (logRecoil)
        {
            const dimensionedScalar recoilForceMag =
                fvc::domainIntegrate(mag(recoilForce));

            Info<< "    Recoil force magnitude: " << recoilForceMag.value()
                << " N" << nl
                << "    Max recoil pressure: " << maxRecoil/1e6
                << " MPa" << endl;
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

    const volVectorField* recoilForceFieldPtr = nullptr;

    if (recoilForceReady)
    {
        recoilForceFieldPtr = &tRecoilForce();
    }

    const scalar domainVolume = gSum(mesh.V());
    const scalar deltaT = runTime.deltaTValue();
    const bool logVelocityDiagnostics =
        verbose && master && (runTime.timeIndex() % 10 == 0);

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

        scalar maxRecoilAccel = 0.0;
        scalar meanRecoilAccel = 0.0;

        if (recoilForceFieldPtr)
        {
            tmp<volVectorField> tRecoilAccel((*recoilForceFieldPtr)/rho);
            const volVectorField& recoilAccel = tRecoilAccel();
            tmp<volScalarField> tRecoilAccelMag = mag(recoilAccel);
            const volScalarField& recoilAccelMag = tRecoilAccelMag();

            const dimensionedScalar recoilAccelIntegral =
                fvc::domainIntegrate(recoilAccelMag);
            maxRecoilAccel = gMax(recoilAccelMag);
            meanRecoilAccel = recoilAccelIntegral.value()/domainVolume;
        }

        tmp<volVectorField> tTotalAccel(tMomentumSource()/rho);
        volVectorField& totalAccel = tTotalAccel.ref();

        if (recoilForceFieldPtr)
        {
            totalAccel += (*recoilForceFieldPtr)/rho;
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
            << "      Max recoil acceleration:   "
            << maxRecoilAccel << " m/s^2 (dU ~ "
            << maxRecoilAccel*deltaT << " m/s per step)" << nl
            << "      Max total acceleration:    "
            << maxTotalAccel << " m/s^2 (dU ~ "
            << maxTotalAccel*deltaT << " m/s per step)" << nl
            << "      Mean total acceleration:   "
            << meanTotalAccel << " m/s^2" << nl
            << "      Mean pressure acceleration: "
            << meanPressureAccel << " m/s^2" << nl
            << "      Mean recoil acceleration:  "
            << meanRecoilAccel << " m/s^2" << endl;
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

const scalar minUEqnDiagVal = solverCoeffs.lookupOrDefault<scalar>
(
    "minUEqnDiag",
    1e-9
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
