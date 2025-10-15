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
const scalar maxPhysicalVelocity = solverCoeffs.lookupOrDefault<Foam::scalar>
(
    "maxVelocity",
    500.0
);
const Foam::scalar maxU = gMax(mag(U)());

if (maxU > maxPhysicalVelocity)
{
    FatalErrorInFunction
        << "Velocity field has DIVERGED!" << nl
        << "Max velocity: " << maxU << " m/s" << nl
        << "Physical limit: " << maxPhysicalVelocity << " m/s" << nl
        << nl
        << "This indicates numerical instability. Check:" << nl
        << "  1. Time step (current: " << runTime.deltaTValue() << " s)" << nl
        << "  2. Max recoil pressure: "
        << (recoilPressurePtr ? gMax(*recoilPressurePtr)/1e6 : 0) << " MPa" << nl
        << "  3. Courant number" << nl
        << "  4. Temperature gradients" << nl
        << abort(FatalError);
}

K = 0.5*magSqr(U);
