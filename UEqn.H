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

// Determine the recoil-pressure clamp from user configuration.  Respect the
// advanced interface-capturing settings so we do not silently override a case
// that deliberately permits >100 MPa recoil loads.
scalar recoilClampMax = 1e8;
bool recoilClampActive = true;

const Foam::dictionary& controlDict = mesh.time().controlDict();

if (controlDict.found("advancedInterfaceCapturing"))
{
    const Foam::dictionary& aicDict = controlDict.subDict("advancedInterfaceCapturing");

    recoilClampActive = aicDict.lookupOrDefault<Foam::Switch>
    (
        "clampRecoil",
        Foam::Switch(true)
    );

    scalar userLimit = aicDict.lookupOrDefault<scalar>
    (
        "maxRecoilPressure",
        aicDict.lookupOrDefault<scalar>("recoilMax", recoilClampMax)
    );

    if (std::isfinite(userLimit) && userLimit > SMALL)
    {
        recoilClampMax = Foam::mag(userLimit);
    }
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

        tmp<volVectorField> tRecoilForce
        (
            (*recoilToUse)*interfaceStrength*interfaceNormal
        );
        volVectorField& recoilForce = tRecoilForce.ref();

        momentumSource += recoilForce;

        if (logRecoil)
        {
            vector totalFaceRecoilForce(vector::zero);

            const scalarField& recoilInternal = (*recoilToUse).internalField();
            const labelUList& owner = mesh.owner();
            const labelUList& neighbour = mesh.neighbour();

            tmp<surfaceScalarField> tInterfaceStrengthFace =
                fvc::interpolate(interfaceStrength);
            const surfaceScalarField& interfaceStrengthFace =
                tInterfaceStrengthFace();

            tmp<surfaceVectorField> tInterfaceNormalFace =
                fvc::interpolate(interfaceNormal);
            const surfaceVectorField& interfaceNormalFace =
                tInterfaceNormalFace();

            const surfaceScalarField& magSf = mesh.magSf();
            const scalarField& interfaceStrengthInternal =
                interfaceStrengthFace.internalField();
            const vectorField& interfaceNormalInternal =
                interfaceNormalFace.internalField();
            const scalarField& magSfInternal = magSf.internalField();

            forAll(owner, facei)
            {
                const label own = owner[facei];
                const label nei = neighbour[facei];

                const scalar pf =
                    0.5*(recoilInternal[own] + recoilInternal[nei]);

                totalFaceRecoilForce +=
                    pf*interfaceStrengthInternal[facei]
                     *interfaceNormalInternal[facei]
                     *magSfInternal[facei];
            }

            forAll(interfaceStrengthFace.boundaryField(), patchi)
            {
                const scalarField& strengthPatch =
                    interfaceStrengthFace.boundaryField()[patchi];
                const vectorField& normalPatch =
                    interfaceNormalFace.boundaryField()[patchi];
                const scalarField& magSfPatch =
                    magSf.boundaryField()[patchi];
                const scalarField& recoilPatch =
                    (*recoilToUse).boundaryField()[patchi];

                forAll(strengthPatch, facei)
                {
                    totalFaceRecoilForce +=
                        recoilPatch[facei]
                      * strengthPatch[facei]
                      * normalPatch[facei]
                      * magSfPatch[facei];
                }
            }

            static_cast<void>(totalFaceRecoilForce);

            const dimensionedScalar recoilForceMag =
                fvc::domainIntegrate(mag(recoilForce));

            Info<< "    Recoil force magnitude: " << recoilForceMag.value()
                << " N" << nl
                << "    Max recoil pressure: " << maxRecoil/1e6
                << " MPa" << endl;
        }

        tRecoilForce.clear();
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
        -gradPrgh

    );

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
