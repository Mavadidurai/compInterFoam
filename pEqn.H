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

    const dictionary& compCoeffs = compInterFoamCoeffsDict(mesh);
    const scalar velocityClamp =
        compCoeffs.lookupOrDefault<scalar>("maxVelocity", 2000.0);
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
        const scalar diagFloorInput =
            compCoeffs.lookupOrDefault<scalar>("minUEqnDiag", 1e-9);
        const scalar diagFloor = Foam::max(diagFloorInput, Foam::SMALL);
        const scalar safeRAUf = 1.0/diagFloor;
        bool adjustedRAUf = false;

        scalarField& rAUfInternal = rAUf.primitiveFieldRef();
        forAll(rAUfInternal, faceI)
        {
            scalar& val = rAUfInternal[faceI];
            if (!std::isfinite(val) || val <= 0)
            {
                val = safeRAUf;
                adjustedRAUf = true;
            }
        }

        surfaceScalarField::Boundary& rAUfBf = rAUf.boundaryFieldRef();
        forAll(rAUfBf, patchI)
        {
            scalarField& patchField = rAUfBf[patchI];
            forAll(patchField, faceI)
            {
                scalar& val = patchField[faceI];
                if (!std::isfinite(val) || val <= 0)
                {
                    val = safeRAUf;
                    adjustedRAUf = true;
                }
            }
        }

        if (adjustedRAUf && verbose && master)
        {
            WarningInFunction << "Clamped non-physical rAUf entries" << endl;
        }
    }
    {
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
    /*
    if (recoilPressurePtr)
    {
        tmp<volVectorField> tRecoilGrad = fvc::grad(*recoilPressurePtr);
        volVectorField& recoilGrad = tRecoilGrad.ref();

        const scalar maxAllowedGrad = 1e4; // Pa/m
        forAll(recoilGrad, cellI)
        {
            const scalar gradMag = mag(recoilGrad[cellI]);
            if (gradMag > maxAllowedGrad)
            {
                recoilGrad[cellI] *= maxAllowedGrad/gradMag;
            }
        }

        HbyA += rAU*recoilGrad;
    }
    */
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
        
    if (recoilPressurePtr)
    {
        const scalar recoilRelaxFactor = 0.01;  // 1% per timestep
        pEqn.source() -= recoilRelaxFactor * (*recoilPressurePtr) * mesh.V();
    }

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

            // Reconstruct the pressure correction contribution only.  The
            // gravitational/capillary term (phig) is already embedded in
            // phiHbyA, so adding it again double-counts the body-force
            // contribution and drives the pressure field into the clamp
            // limits.  Recover the standard OpenFOAM formulation by using
            // just the pressure-correction flux.
            U = HbyA + rAU*fvc::reconstruct(pEqn.flux()/rAUf);
            fvOptions.correct(U);
            if (velocityClamp > 0)
            {
                vectorField& Uinternal = U.primitiveFieldRef();
                forAll(Uinternal, cellI)
                {
                    vector& Ui = Uinternal[cellI];
                    const scalar magU = mag(Ui);
                    if (std::isfinite(magU))
                    {
                        if (magU > velocityClamp)
                        {
                            Ui *= velocityClamp/magU;
                        }
                    }
                    else
                    {
                        Ui = vector::zero;
                    }
                }
            }
            U.correctBoundaryConditions();
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
   const dictionary& pressureCoeffs = compInterFoamCoeffsDict(mesh);
    const Switch pressureClamp = pressureCoeffs.lookupOrDefault<Switch>
    (
        "pressureClamp",
        false
    );

    if (pressureClamp)
    {
        const scalar minPressure = pressureCoeffs.lookupOrDefault<scalar>
        (
            "minPressure",
            -GREAT
        );
        const scalar maxPressure = pressureCoeffs.lookupOrDefault<scalar>
        (
            "maxPressure",
            GREAT
        );

        auto clampField =
        [
            &minPressure,
            &maxPressure
        ](volScalarField& fld)
        {
            bool clamped = false;

            scalarField& internal = fld.primitiveFieldRef();
            forAll(internal, cellI)
            {
                scalar& val = internal[cellI];
                const scalar limited = Foam::max
                (
                    minPressure,
                    Foam::min(maxPressure, val)
                );

                if (limited != val)
                {
                    val = limited;
                    clamped = true;
                }
            }

            volScalarField::Boundary& bFld = fld.boundaryFieldRef();
            forAll(bFld, patchI)
            {
                scalarField& patchField = bFld[patchI];
                forAll(patchField, faceI)
                {
                    scalar& val = patchField[faceI];
                    const scalar limited = Foam::max
                    (
                        minPressure,
                        Foam::min(maxPressure, val)
                    );

                    if (limited != val)
                    {
                        val = limited;
                        clamped = true;
                    }
                }
            }

            return clamped;
        };

        const bool pClamped = clampField(p);
        const bool pRghClamped = clampField(p_rgh);

        if (pClamped || pRghClamped)
        {
            p.correctBoundaryConditions();
            p_rgh.correctBoundaryConditions();

            if (master)
            {
                WarningInFunction
                    << "Pressure clamp applied to "
                    << (pClamped && pRghClamped
                        ? "p and p_rgh"
                        : (pClamped ? "p" : "p_rgh"))
                    << " with limits [" << minPressure
                    << ", " << maxPressure << "]" << endl;
            }
        }
    }
    K = 0.5*magSqr(U);

    if (verbose && master)
    {
        Info<< "Pressure solve complete: max p = " << max(p).value()
            << " Pa, min p = " << min(p).value() << " Pa" << endl;
    }
}
