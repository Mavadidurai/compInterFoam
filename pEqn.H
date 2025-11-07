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

constrainPressure(p_rgh, U, phiHbyA, rAUf);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUf, p_rgh)
     == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);
    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        const Foam::dictionary& clampCoeffs = compInterFoamCoeffsDict(mesh);
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
            const scalar configuredClampMax =
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
            if
            (
                configuredClampMax > Foam::SMALL
             && configuredClampMax < 0.99*Foam::GREAT
            )
            {
                recoilClampTarget = Foam::min(recoilClampTarget, configuredClampMax);
            }

            // Raise pressure clamp to accommodate recoil if needed
            if
            (
                recoilClampTarget > Foam::SMALL
             && maxPressure.value() + Foam::SMALL < recoilClampTarget
            )
            {
                if (Foam::Pstream::master())
                {
                    Info<< "    pressureClamp: raising maxPressure from "
                        << maxPressure.value() << " Pa to recoil limit "
                        << recoilClampTarget << " Pa" << Foam::nl;
                }

                maxPressure = Foam::dimensionedScalar
                (
                    "maxPressure",
                    Foam::dimPressure,
                    recoilClampTarget
                );
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

if (pressureClamp)
{
    // Read recoil configuration consistently (same logic as above and UEqn.H)
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
    if
    (
        configuredMaxPressure > Foam::SMALL
     && configuredMaxPressure < 0.99*Foam::GREAT
    )
    {
        recoilClampTarget = Foam::min(recoilClampTarget, configuredMaxPressure);
    }

    // Raise pressure clamp to accommodate recoil if needed
    if (recoilClampTarget > Foam::SMALL && maxPressure.value() + Foam::SMALL < recoilClampTarget)
    {
        maxPressure = Foam::dimensionedScalar
        (
            "maxPressure",
            Foam::dimPressure,
            recoilClampTarget
        );
    }
}

const Foam::dimensionedScalar fallbackMinPressure
(
    "fallbackMinPressure",
    Foam::dimPressure,
    -1e7
);
const Foam::dimensionedScalar fallbackMaxPressure
(
    "fallbackMaxPressure",
    Foam::dimPressure,
    1e9
);

if (minPressure.value() <= -0.99*Foam::GREAT)
{
    minPressure = fallbackMinPressure;
}

if (maxPressure.value() >= 0.99*Foam::GREAT)
{
    maxPressure = fallbackMaxPressure;
}

const Foam::dimensionedScalar observedMinP = Foam::gMin(p);
const Foam::dimensionedScalar observedMaxP = Foam::gMax(p);

const bool pressureNonFinite =
    (!std::isfinite(observedMinP.value()) || !std::isfinite(observedMaxP.value()));

const bool pressureOutOfBounds =
    (observedMinP.value() < minPressure.value() || observedMaxP.value() > maxPressure.value());

if (!pressureClamp && (pressureNonFinite || pressureOutOfBounds))
{
    WarningInFunction
        << "Pressure field exceeded physical range while pressureClamp=false"
        << Foam::nl
        << "    min(p) = " << observedMinP.value() << " Pa, max(p) = "
        << observedMaxP.value() << " Pa" << Foam::nl
        << "    Applying emergency clamp to [" << minPressure.value() << ", "
        << maxPressure.value() << "] Pa" << Foam::endl;

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
