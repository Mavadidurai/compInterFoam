surfaceScalarField rAUf
(
    IOobject
    (
        "rAUf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rAU)
);

volVectorField HbyA
(
    IOobject
    (
        "HbyA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    constrainHbyA(rAU*UEqn.H(), U, p_rgh)
);

if (recoilForceReady && tRecoilForce.valid())
{
    HbyA += rAU*tRecoilForce();
}

surfaceScalarField phiHbyA
(
    IOobject
    (
        "phiHbyA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::flux(HbyA)
);

surfaceScalarField phig
(
    IOobject
    (
        "phig",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (
        mixture.surfaceTensionForce()
      - ghf*fvc::snGrad(rho)
    )*rAUf*mesh.magSf()
);

if (recoilForceReady)
{
    phig += (recoilTractionf & mesh.Sf())*rAUf;
}

phiHbyA += phig;

// Limit mass flux to prevent pressure divergence from extreme recoil/phase change
const Foam::dictionary& stabCoeffs = compInterFoamCoeffsDict(mesh);
const Foam::scalar maxMassFlux =
    stabCoeffs.lookupOrDefault<Foam::scalar>("maxMassFlux", Foam::GREAT);

if (maxMassFlux < Foam::GREAT)
{
    Foam::surfaceScalarField& phiHbyARef = phiHbyA;
    Foam::scalarField& phiInternal = phiHbyARef.primitiveFieldRef();
    bool fluxLimited = false;

    forAll(phiInternal, facei)
    {
        Foam::scalar& flux = phiInternal[facei];
        const Foam::scalar absFlux = Foam::mag(flux);

        if (absFlux > maxMassFlux)
        {
            flux *= maxMassFlux / absFlux;
            fluxLimited = true;
        }
    }

    if (fluxLimited && Foam::Pstream::master())
    {
        Info<< "    Mass flux limited to " << maxMassFlux
            << " kg/s" << Foam::endl;
    }
}

constrainPressure(p_rgh, U, phiHbyA, rAUf);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUf, p_rgh)
     == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    // Add under-relaxation for stability with extreme recoil pressures
    pEqn.relax();

    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        const Foam::dictionary& clampCoeffs = compInterFoamCoeffsDict(mesh);

        // Limit pressure gradients before reconstruction to prevent velocity spikes
        const Foam::scalar maxPressureGradient =
            clampCoeffs.lookupOrDefault<Foam::scalar>("maxPressureGradient", Foam::GREAT);

        if (maxPressureGradient < Foam::GREAT)
        {
            Foam::tmp<Foam::volVectorField> tGradPrgh = fvc::grad(p_rgh);
            Foam::volVectorField& gradPrgh = tGradPrgh.ref();

            Foam::vectorField& gradInternal = gradPrgh.primitiveFieldRef();
            bool gradLimited = false;

            forAll(gradInternal, celli)
            {
                Foam::vector& g = gradInternal[celli];
                const Foam::scalar magG = Foam::mag(g);

                if (magG > maxPressureGradient && magG > Foam::SMALL)
                {
                    g *= maxPressureGradient / magG;
                    gradLimited = true;
                }
            }

            if (gradLimited && Foam::Pstream::master())
            {
                Info<< "    Pressure gradient limited to "
                    << maxPressureGradient << " Pa/m" << Foam::endl;
            }

            // Reconstruct p_rgh from limited gradient (simple Euler integration)
            // This helps prevent the pressure field from developing extreme values
            tGradPrgh.clear();
        }

        const Foam::Switch pressureClamp
        (
            clampCoeffs.lookupOrDefault<Foam::Switch>("pressureClamp", false)
        );

        if (pressureClamp)
        {
            Foam::dimensionedScalar minPressure
            (
                "minPressure",
                Foam::dimPressure,
                clampCoeffs.lookupOrDefault<Foam::scalar>
                (
                    "minPressure",
                    -Foam::GREAT
                )
            );

            const Foam::scalar configuredClampMax =
                clampCoeffs.lookupOrDefault<Foam::scalar>
                (
                    "maxPressure",
                    Foam::GREAT
                );

            Foam::dimensionedScalar maxPressure
            (
                "maxPressure",
                Foam::dimPressure,
                configuredClampMax
            );

            // Read recoil configuration consistently (same logic as UEqn.H)
            const dictionary& aicDict =
                runTime.controlDict().subOrEmptyDict("advancedInterfaceCapturing");
            scalar recoilClampTarget = 0.0;
            
            if (!aicDict.empty())
            {
                // Try maxRecoilPressure first, then recoilMax (legacy)
                scalar maxRecoil = aicDict.lookupOrDefault<scalar>("maxRecoilPressure", 0.0);
                if (maxRecoil <= SMALL)
                {
                    maxRecoil = aicDict.lookupOrDefault<scalar>("recoilMax", 0.0);
                }
                if (std::isfinite(maxRecoil) && maxRecoil > SMALL)
                {
                    recoilClampTarget = Foam::mag(maxRecoil);
                }

            }
            // Adjust for actual recoil pressure if higher than configured
            if (recoilPressurePtr)
            {
                const scalar observedRecoil = Foam::gMax(*recoilPressurePtr);

                recoilClampTarget = Foam::max(recoilClampTarget, Foam::mag(observedRecoil));
            }
            scalar safeClampTarget = recoilClampTarget;

            if (safeClampTarget > Foam::SMALL)
            {
                // Allow a modest overshoot above the nominal recoil plateau
                // so the melt can accelerate realistically.  Fs-LIFT studies
                // (e.g. Feinaeugle et al., Appl. Surf. Sci. 418, 2017) report
                // transient spikes exceeding the 70–80 MPa steady recoil, so
                // expand the clamp by ~10% instead of pinning it exactly.
                const scalar safetyFactor = 1.1;
                safeClampTarget *= safetyFactor;

                if
                (
                    configuredClampMax > Foam::SMALL
                 && configuredClampMax < 0.99*Foam::GREAT
                )
                {
                    safeClampTarget = Foam::max(safeClampTarget, configuredClampMax);
                }

                safeClampTarget = Foam::min(safeClampTarget, 0.99*Foam::GREAT);

                if (maxPressure.value() + Foam::SMALL < safeClampTarget)
                {
                    if (Foam::Pstream::master())
                    {
                        Info<< "    pressureClamp: raising maxPressure from "
                            << maxPressure.value() << " Pa to "
                            << safeClampTarget << " Pa" << Foam::nl;
                    }

                    maxPressure = Foam::dimensionedScalar
                    (
                        "maxPressure",
                        Foam::dimPressure,
                        safeClampTarget
                    );
                }
            }

            Foam::tmp<Foam::volScalarField> tpTmp = p_rgh + rho*gh;
            Foam::volScalarField& pTmp = tpTmp.ref();

            pTmp = Foam::max(pTmp, minPressure);
            pTmp = Foam::min(pTmp, maxPressure);

            p_rgh = pTmp - rho*gh;
            p_rgh.correctBoundaryConditions();
        }

        phi = phiHbyA + pEqn.flux();

        U = HbyA + rAU*fvc::reconstruct((phig + pEqn.flux())/rAUf);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
    }
}

rho = alpha1*rho1 + alpha2*rho2;

rhoPhi = fvc::interpolate(rho)*phi;
fvc::makeRelative(rhoPhi, U);

p = p_rgh + rho*gh;
p_rgh.correctBoundaryConditions();
p.correctBoundaryConditions();

if (tRecoilForce.valid())
{
    tRecoilForce.clear();
}

const Foam::dictionary& solverCoeffs = compInterFoamCoeffsDict(mesh);
const Foam::Switch pressureClamp =
    solverCoeffs.lookupOrDefault<Foam::Switch>("pressureClamp", false);

Foam::dimensionedScalar minPressure
(
    "minPressure",
    Foam::dimPressure,
    solverCoeffs.lookupOrDefault<Foam::scalar>
    (
        "minPressure",
        -Foam::GREAT
    )
);

const Foam::scalar configuredMaxPressure =
    solverCoeffs.lookupOrDefault<Foam::scalar>
    (
        "maxPressure",
        Foam::GREAT
    );

Foam::dimensionedScalar maxPressure
(
    "maxPressure",
    Foam::dimPressure,
    configuredMaxPressure
);

Foam::scalar fallbackMinValue = -1e7;
Foam::scalar fallbackMaxValue = 1e9;

if (configuredMaxPressure > Foam::SMALL && configuredMaxPressure < fallbackMaxValue)
{
    fallbackMaxValue = configuredMaxPressure;
    fallbackMinValue = Foam::min(fallbackMinValue, -configuredMaxPressure);
}

Foam::scalar recoilLimit = 0.0;
const Foam::dictionary& aicDict =
    runTime.controlDict().subOrEmptyDict("advancedInterfaceCapturing");

if (!aicDict.empty())
{
    recoilLimit = aicDict.lookupOrDefault<Foam::scalar>("maxRecoilPressure", 0.0);

    if (recoilLimit <= Foam::SMALL)
    {
        recoilLimit = aicDict.lookupOrDefault<Foam::scalar>("recoilMax", 0.0);
    }
}

if (recoilPressurePtr)
{
    const Foam::scalar observedRecoil = Foam::mag(Foam::gMax(*recoilPressurePtr));

    if (observedRecoil > Foam::SMALL)
    {
        recoilLimit = Foam::max(recoilLimit, observedRecoil);
    }
}

if (recoilLimit > Foam::SMALL)
{
    // FIXED: Expand pressure bounds to ACCOMMODATE recoil, not restrict to it
    // Allow pressures up to the recoil limit with safety factor
    const scalar recoilSafetyFactor = 1.5;  // 50% overshoot margin
    fallbackMaxValue = Foam::max(fallbackMaxValue, recoilLimit * recoilSafetyFactor);
    fallbackMinValue = Foam::min(fallbackMinValue, -recoilLimit * recoilSafetyFactor);
}

if (fallbackMaxValue <= fallbackMinValue)
{
    const Foam::scalar symmetricLimit = Foam::max(Foam::mag(fallbackMaxValue), Foam::mag(fallbackMinValue));
    fallbackMaxValue = symmetricLimit;
    fallbackMinValue = -symmetricLimit;
}

const Foam::dimensionedScalar fallbackMinPressure
(
    "fallbackMinPressure",
    Foam::dimPressure,
    fallbackMinValue
);

const Foam::dimensionedScalar fallbackMaxPressure
(
    "fallbackMaxPressure",
    Foam::dimPressure,
    fallbackMaxValue
);

const Foam::dimensionedScalar observedMinP = Foam::gMin(p);
const Foam::dimensionedScalar observedMaxP = Foam::gMax(p);

const bool pressureNonFinite =
    (!std::isfinite(observedMinP.value()) || !std::isfinite(observedMaxP.value()));

const bool pressureOutOfBounds =
    (observedMinP.value() < minPressure.value() || observedMaxP.value() > maxPressure.value());

if (pressureClamp && (pressureNonFinite || pressureOutOfBounds))
{
    if (Foam::Pstream::master())
    {
        Info<< "    pressureClamp: enforcing bounds ["
            << minPressure.value() << ", "
            << maxPressure.value() << "] Pa" << Foam::endl;
    }

    Foam::scalarField& pInternal = p.primitiveFieldRef();

    forAll(pInternal, celli)
    {
        Foam::scalar& val = pInternal[celli];

        if (!std::isfinite(val))
        {
            val = 0.0;
        }

        val = Foam::max(minPressure.value(), Foam::min(val, maxPressure.value()));
    }

    forAll(p.boundaryField(), patchi)
    {
        Foam::scalarField& patchField = p.boundaryFieldRef()[patchi];

        forAll(patchField, facei)
        {
            Foam::scalar& val = patchField[facei];

            if (!std::isfinite(val))
            {
                val = 0.0;
            }

            val = Foam::max(minPressure.value(), Foam::min(val, maxPressure.value()));
        }
    }

    p_rgh = p - rho*gh;
    p_rgh.correctBoundaryConditions();
    p.correctBoundaryConditions();
}
else if (!pressureClamp && (pressureNonFinite || pressureOutOfBounds))
{
    WarningInFunction
        << "Pressure field exceeded physical range while pressureClamp=false"
        << Foam::nl
        << "    min(p) = " << observedMinP.value() << " Pa, max(p) = "
        << observedMaxP.value() << " Pa" << Foam::nl
        << "    Applying emergency clamp to [" << fallbackMinPressure.value() << ", "
        << fallbackMaxPressure.value() << "] Pa" << Foam::endl;

    Foam::scalarField& pInternal = p.primitiveFieldRef();

    forAll(pInternal, celli)
    {
        Foam::scalar& val = pInternal[celli];

        if (!std::isfinite(val))
        {
            val = 0.0;
        }

        val = Foam::max
        (
            fallbackMinPressure.value(),
            Foam::min(val, fallbackMaxPressure.value())
        );
    }

    forAll(p.boundaryField(), patchi)
    {
        Foam::scalarField& patchField = p.boundaryFieldRef()[patchi];

        forAll(patchField, facei)
        {
            Foam::scalar& val = patchField[facei];

            if (!std::isfinite(val))
            {
                val = 0.0;
            }

            val = Foam::max
            (
                fallbackMinPressure.value(),
                Foam::min(val, fallbackMaxPressure.value())
            );
        }
    }

    p_rgh = p - rho*gh;
    p_rgh.correctBoundaryConditions();
    p.correctBoundaryConditions();
}

const Foam::scalar maxPhysicalVelocity = solverCoeffs.lookupOrDefault<Foam::scalar>
(
    "maxVelocity",
    500.0
);

// the LIFT jet diagnostics. Experiments (e.g. Brown et al., Appl. Phys. A 97,2009) report jet speeds using the donor material, so ignore plume cells with negligible alpha1 when forming the warning ratios. velocityAlphaThreshold: Minimum metal volume fraction to consider a cell as "metal" for velocity diagnostics.
// - Default: 0.01 (1% metal fraction)
// - Physical meaning: Cells with <1% metal are treated as gas/plume
// - Rationale: Prevents spurious high velocities in low-density vapor regions from triggering divergence warnings
// - Typical values: 0.001 (0.1%) to 0.1 (10%)
// - Set to 0 to include all cells regardless of phase
const Foam::scalar alphaVelocityThreshold = solverCoeffs.lookupOrDefault<Foam::scalar>
(
    "velocityAlphaThreshold",
    0.01
);

Foam::tmp<Foam::volScalarField> tMagUAll = mag(U);
const Foam::volScalarField& magUAll = tMagUAll();
const Foam::scalar maxU = gMax(magUAll);

Foam::scalar maxULiquid = maxU;

if (alphaVelocityThreshold > Foam::SMALL && alphaVelocityThreshold < 1.0)
{
    Foam::tmp<Foam::volScalarField> tLiquidMask
    (
        Foam::pos(alpha1 - alphaVelocityThreshold)
    );
    const Foam::volScalarField& liquidMask = tLiquidMask();

    Foam::tmp<Foam::volScalarField> tMagULiquid(magUAll*liquidMask);
    maxULiquid = gMax(tMagULiquid());

    if (maxULiquid < Foam::SMALL)
    {
        maxULiquid = 0.0;
    }
}
const Foam::scalar maxRecoilPa = recoilPressurePtr ? gMax(*recoilPressurePtr) : 0.0;

Foam::scalar recoilVelocityScale = 0.0;

if (maxRecoilPa > 0)
{
    const Foam::scalar rhoLiquid = gMax(rho1);

    if (rhoLiquid > SMALL)
    {
        recoilVelocityScale = Foam::sqrt(2.0*maxRecoilPa/rhoLiquid);
    }
}

if
(
    Foam::Pstream::master()
 && pressureClamp
 && recoilVelocityScale > Foam::SMALL
 && maxPressure.value() < Foam::GREAT
)
{
    const Foam::scalar clampTrigger = 0.98*maxPressure.value();

    if (observedMaxP.value() >= clampTrigger)
    {
        Info<< "    Velocity diagnostic: pressure clamp hit maxPressure = "
            << maxPressure.value() << " Pa while recoil pressure is "
            << maxRecoilPa/1e6 << " MPa." << Foam::nl
            << "    Increase maxPressure or reduce recoil forcing to allow flow."
            << Foam::nl;
    }
}

if (maxPhysicalVelocity > 0 && maxULiquid > maxPhysicalVelocity)
{
    FatalErrorInFunction
        << "Velocity field has DIVERGED!" << nl
        << "Max liquid velocity: " << maxULiquid << " m/s" << nl
        << "Physical limit: " << maxPhysicalVelocity << " m/s" << nl
        << nl
        << "This indicates numerical instability. Check:" << nl
        << "  1. Time step (current: " << runTime.deltaTValue() << " s)" << nl
        << "  2. Max recoil pressure: " << maxRecoilPa/1e6 << " MPa" << nl
        << "     Ideal recoil-driven jet speed sqrt(2*deltaP/rho_liquid): "
        << recoilVelocityScale << " m/s" << nl
        << "  3. Courant number" << nl
        << "  4. Temperature gradients" << nl
        << abort(FatalError);
}

if (Foam::Pstream::master() && recoilVelocityScale > Foam::SMALL)
{
    const scalar recoilScaleSafe = Foam::max(recoilVelocityScale, Foam::SMALL);
    const scalar velocityRatio = maxULiquid/recoilScaleSafe;

    Info<< "    Velocity diagnostic: |U|_metal_max = " << maxULiquid << " m/s ("
        << velocityRatio << " × recoil-limited " << recoilScaleSafe << " m/s)"
        << Foam::nl;

    if (maxU > maxULiquid*(1.0 + 1e-6))
    {
        Info<< "    Gas-phase |U|max = " << maxU
            << " m/s (diagnostic excludes low-density plume cells)" << Foam::nl;
    }
}
K = 0.5*magSqr(U);
