/*---------------------------------------------------------------------------*\\
  Standard pressure equation adapted for LIFT extensions
\\*---------------------------------------------------------------------------*/

extern Foam::Switch verbose;

{
    if (verbose)
    {
        Info<< "Solving pressure equation" << endl;
    }

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
            fvc::interpolate(mixture.sigmaK())
          * fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    constrainPressure(p_rgh, U, phiHbyA, rAUf);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
        );

        pEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

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

    K = 0.5*magSqr(U);

    if (verbose)
    {
        Info<< "Pressure solve complete: max p = " << max(p).value()
            << " Pa, min p = " << min(p).value() << " Pa" << endl;
    }
}
