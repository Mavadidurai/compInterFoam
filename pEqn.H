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
            const Foam::scalar minPressure = clampCoeffs.lookupOrDefault<Foam::scalar>
            (
                "minPressure",
                -Foam::GREAT
            );
            const Foam::scalar maxPressure = clampCoeffs.lookupOrDefault<Foam::scalar>
            (
                "maxPressure",
                Foam::GREAT
            );

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

const Foam::dictionary& solverCoeffs = compInterFoamCoeffsDict(mesh);
const Foam::scalar maxPhysicalVelocity = solverCoeffs.lookupOrDefault<Foam::scalar>
(
    "maxVelocity",
    500.0
);

const Foam::scalar maxU = gMax(mag(U)());
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

if (maxPhysicalVelocity > 0 && maxU > maxPhysicalVelocity)
{
    FatalErrorInFunction
        << "Velocity field has DIVERGED!" << nl
        << "Max velocity: " << maxU << " m/s" << nl
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

K = 0.5*magSqr(U);
