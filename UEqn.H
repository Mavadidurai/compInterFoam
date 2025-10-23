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

    if (recoilPressurePtr)
    {
        const volScalarField& recoilP = *recoilPressurePtr;
        
        // Validate recoil pressure before using
        const scalar maxRecoil = gMax(recoilP);
        const scalar minRecoil = gMin(recoilP);
        
        if (maxRecoil > 1e8 || minRecoil < -1e8)
        {
            WarningInFunction
                << "Recoil pressure out of physical range: ["
                << minRecoil << ", " << maxRecoil << "] Pa"
                << nl << "Clamping to [-100 MPa, 100 MPa]"
                << endl;
                
            volScalarField recoilClamped
            (
                Foam::min
                (
                    Foam::max(recoilP, dimensionedScalar("minR", dimPressure, -1e8)),
                    dimensionedScalar("maxR", dimPressure, 1e8)
                )
            );
            
            momentumSource -= fvc::grad(alpha1*recoilClamped);
        }
        else
        {
            momentumSource -= fvc::grad(alpha1*recoilP);
        }
        
        if (verbose && master && runTime.timeIndex() % 10 == 0)
        {
            const dimensionedScalar recoilForce = 
                fvc::domainIntegrate(mag(fvc::grad(alpha1*recoilP)));
                
            Info<< "    Recoil force magnitude: " << recoilForce.value() 
                << " N" << nl
                << "    Max recoil pressure: " << maxRecoil/1e6 
                << " MPa" << endl;
        }
    }

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
