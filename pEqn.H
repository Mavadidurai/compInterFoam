/*---------------------------------------------------------------------------*\\
  Standard pressure equation adapted for LIFT extensions
\\*---------------------------------------------------------------------------*/

extern Foam::Switch verbose;
extern const bool master;
{
    if (verbose && master)
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
    {
        const dictionary& compCoeffs = compInterFoamCoeffsDict(mesh);
        const Switch enableClamp = compCoeffs.lookupOrDefault<Switch>("enableRAUClamp", false);

        if (enableClamp)
        {
            const scalar minRAUf = compCoeffs.lookupOrDefault<scalar>("minRAUf", 1e-14);
            const scalar maxRAUf = compCoeffs.lookupOrDefault<scalar>("maxRAUf", GREAT);

            auto clampValue = [&](scalar& value)
            {
                if (!std::isfinite(value))
                {
                    value = minRAUf;
                }
                else
                {
                    value = Foam::max(minRAUf, Foam::min(maxRAUf, value));
                }
            };

            scalarField& rAUfInternal = rAUf.primitiveFieldRef();
            forAll(rAUfInternal, faceI)
            {
                clampValue(rAUfInternal[faceI]);
            }

            surfaceScalarField::Boundary& rAUfBf = rAUf.boundaryFieldRef();
            forAll(rAUfBf, patchI)
            {
                scalarField& patchField = rAUfBf[patchI];
                forAll(patchField, faceI)
                {
                    clampValue(patchField[faceI]);
                }
            }
        }
    }
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
            fvm::laplacian(rAUf, p_rgh)
         == fvc::div(phiHbyA)
        );

        SolverPerformance<scalar> solverPerf =
            pEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if
        (
            !solverPerf.converged()
         && solverPerf.finalResidual() > 1e-2
        )
        {
            FatalErrorInFunction
                << "Pressure solve failed: converged="
                << solverPerf.converged()
                << ", final residual=" << solverPerf.finalResidual()
                << exit(FatalError);
        }

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();

            U = HbyA + rAU*fvc::reconstruct((phig + pEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
            const scalar maxRealisticVel = 2000.0;
            forAll(U, cellI)
            {
                const scalar magU = mag(U[cellI]);
                if (std::isfinite(magU))
                {
                    if (magU > maxRealisticVel)
                    {
                        U[cellI] *= maxRealisticVel/magU;
                    }
                }
                else
                {
                    U[cellI] = vector::zero;
                }
            }
            phi = fvc::flux(U);
            if (verbose && master)
            {
                Info<< "DIAGNOSTIC AFTER PRESSURE CORRECTION:" << nl
                    << "  max(U) = " << max(mag(U)).value() << " m/s" << nl
                    << "  max(phi) = " << max(mag(phi)).value() << " m^3/s" << nl;
	    }           
        }
    }

    rho = alpha1*rho1 + alpha2*rho2;

    rhoPhi = fvc::interpolate(rho)*phi;
    fvc::makeRelative(rhoPhi, U);

    p = p_rgh + rho*gh;
    p_rgh.correctBoundaryConditions();
    p.correctBoundaryConditions();

    K = 0.5*magSqr(U);

    if (verbose && master)
    {
        Info<< "Pressure solve complete: max p = " << max(p).value()
            << " Pa, min p = " << min(p).value() << " Pa" << endl;
    }
}
