/*---------------------------------------------------------------------------*\
  Robust pressure equation with effective numerical safeguards
\*---------------------------------------------------------------------------*/
#include <cmath>
{
    Info<< "Solving pressure equation" << endl;
    
    // Create robust interpolation of rAU
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    
// Apply bounds to rAUf to prevent division by zero or infinities
    // Clamp values can be specified in compInterFoamCoeffs or controlDict
    const dictionary& solverControls = mesh.solutionDict().subOrEmptyDict
    ("compInterFoamCoeffs");
    const dictionary& controlDict = mesh.time().controlDict();
    scalar minRAUf = solverControls.lookupOrDefault<scalar>
    (
        "minRAUf",
        controlDict.lookupOrDefault<scalar>("minRAUf", 1e-10)
    );
    scalar maxRAUf = solverControls.lookupOrDefault<scalar>
    (
        "maxRAUf",
        controlDict.lookupOrDefault<scalar>("maxRAUf", 1e6)
    );
    const Switch enableRAUfClamp = solverControls.lookupOrDefault<Switch>
    (
        "enableRAUfClamp",
        controlDict.lookupOrDefault<Switch>("enableRAUfClamp", true)
    );

    if (!std::isfinite(minRAUf))
    {
        minRAUf = 1e-10;
    }

    if (!std::isfinite(maxRAUf))
    {
        maxRAUf = 1e6;
    }

    label minRAUfCount = 0;
    label maxRAUfCount = 0;

    // Internal field
    scalarField& rAUfIn = rAUf.primitiveFieldRef();
    forAll(rAUfIn, faceI)
    {
        if (!std::isfinite(rAUfIn[faceI]) || (enableRAUfClamp && rAUfIn[faceI] < minRAUf))
        {
            if (enableRAUfClamp)
            {
                rAUfIn[faceI] = minRAUf;
                minRAUfCount++;
            }
        }
        else if (enableRAUfClamp && rAUfIn[faceI] > maxRAUf)
        {
            rAUfIn[faceI] = maxRAUf;
            maxRAUfCount++;
        }
    }

    // Boundary fields
    forAll(rAUf.boundaryField(), patchI)
    {
        fvsPatchScalarField& pRauf = rAUf.boundaryFieldRef()[patchI];
        forAll(pRauf, faceI)
        {
            if (!std::isfinite(pRauf[faceI]) || (enableRAUfClamp && pRauf[faceI] < minRAUf))
            {
                if (enableRAUfClamp)
                {
                    pRauf[faceI] = minRAUf;
                    minRAUfCount++;
                }
            }
            else if (enableRAUfClamp && pRauf[faceI] > maxRAUf)
            {
                pRauf[faceI] = maxRAUf;
                maxRAUfCount++;
            }
        }
    }

    reduce(minRAUfCount, sumOp<label>());
    reduce(maxRAUfCount, sumOp<label>());
    if (enableRAUfClamp && (minRAUfCount > 0 || maxRAUfCount > 0))
    {
        Info<< "rAUf clamping stats: " << minRAUfCount
            << " faces to minRAUf (" << minRAUf << "), "
            << maxRAUfCount << " faces to maxRAUf (" << maxRAUf << ")" << endl;
    }
    else if (!enableRAUfClamp &&
        (min(rAUf).value() < minRAUf || max(rAUf).value() > maxRAUf))
    {
        Info<< "rAUf outside expected range [" << minRAUf << ", "
            << maxRAUf << "]: min=" << min(rAUf).value()
            << " max=" << max(rAUf).value() << endl;
    }

 // Construct HbyA and corresponding flux
    volVectorField HbyA(
        IOobject("HbyA", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        constrainHbyA(rAU*UEqn.H(), U, p_rgh)
    );

    surfaceScalarField phiHbyA(
        IOobject("phiHbyA", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        (fvc::interpolate(HbyA) & mesh.Sf())
    );

    constrainPressure(p_rgh, U, phiHbyA, rAUf);

    fvScalarMatrix pEqn(fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA));
    
    // Set reference and apply adaptive under-relaxation
    pEqn.setReference(pRefCell, pRefValue);
    static scalar relaxP
    (
        mesh.solutionDict().lookupOrDefault<scalar>("relaxP", 0.3)
    );
    const scalar targetResP = mesh.solutionDict().lookupOrDefault<scalar>(
        "targetResP", 0.1);
    pEqn.relax(relaxP);
    
    // Create backups in case solver diverges
    volScalarField p_rgh_old(p_rgh);
    surfaceScalarField phi_old(phi);
    volVectorField U_old(U);

    // Solve and inspect initial residual
    SolverPerformance<scalar> solverPerf = pEqn.solve();
    scalar residual = solverPerf.initialResidual();
    Info<< "pEqn initial residual: " << residual << endl;
    relaxP = max(min(relaxP * targetResP/(residual + SMALL), 0.7), 0.1);
    Info<< "  Updated relaxP: " << relaxP << endl;

    // Revert to previous solution if residual too large
    if (residual >= 1e3)
    {
        Warning << "Pressure solve appears unstable - reverting to previous solution" << endl;
        p_rgh = p_rgh_old;
        phi = phi_old;
        U = U_old;
    }

    
    // Check for invalid values after solve
    bool hasInvalidValues = false;
    forAll(p_rgh, cellI)
    {
        if (!std::isfinite(p_rgh[cellI]))
        {
            hasInvalidValues = true;
            p_rgh[cellI] = p_rgh_old[cellI];
        }
    }
    
    if (hasInvalidValues)
    {
        Warning << "Invalid values detected in p_rgh after solve - fixed" << endl;
    }
    
    // Add bounds to pressure to prevent divergence
    const bool pressureClamp = runTime.controlDict().lookupOrDefault<bool>
    (
        "pressureClamp",
        true
    );

    const scalar maxPressure = runTime.controlDict().lookupOrDefault<scalar>
    (
        "maxPressure",
        1e6
    );

    const scalar minPressure = runTime.controlDict().lookupOrDefault<scalar>
    (
        "minPressure",
        1e3
    );

   // Apply pressure limits only when enabled and count clamped cells
    label minPressureCount = 0;
    label maxPressureCount = 0;
    if (pressureClamp)
    {
        dimensionedScalar maxP("maxP", dimPressure, maxPressure);
        dimensionedScalar minP("minP", dimPressure, minPressure);

        scalar maxPVal = maxP.value();
        scalar minPVal = minP.value();

        // Internal field
        scalarField& pRghIn = p_rgh.primitiveFieldRef();
        forAll(pRghIn, cellI)
        {
            if (pRghIn[cellI] < minPVal)
            {
                pRghIn[cellI] = minPVal;
                minPressureCount++;
            }
            else if (pRghIn[cellI] > maxPVal)
            {
                pRghIn[cellI] = maxPVal;
                maxPressureCount++;
            }
        }

        // Boundary fields
        forAll(p_rgh.boundaryField(), patchI)
        {
            fvPatchScalarField& pPatch = p_rgh.boundaryFieldRef()[patchI];
            forAll(pPatch, faceI)
            {
                if (pPatch[faceI] < minPVal)
                {
                    pPatch[faceI] = minPVal;
                    minPressureCount++;
                }
                else if (pPatch[faceI] > maxPVal)
                {
                    pPatch[faceI] = maxPVal;
                    maxPressureCount++;
                }
            }
        }

        reduce(minPressureCount, sumOp<label>());
        reduce(maxPressureCount, sumOp<label>());

        Info<< "Pressure clamping stats: " << minPressureCount
            << " cells to minPressure (" << minPressure << "), "
            << maxPressureCount << " cells to maxPressure (" << maxPressure
            << ")" << endl;
    }

    // Update pressure field
    p = p_rgh + rho*gh;
    
    // Ensure boundary conditions
    p_rgh.correctBoundaryConditions();
    p.correctBoundaryConditions();
    
    // Update density
    rho = alpha1*rho1 + alpha2*rho2;

    // Apply limits to density for safety and make bounds configurable
    const bool densityClamp = runTime.controlDict().lookupOrDefault<bool>
    (
        "densityClamp",
        true
    );

    const scalar maxDensity = runTime.controlDict().lookupOrDefault<scalar>
    (
        "maxDensity",
        10000.0
    );

    const scalar minDensity = runTime.controlDict().lookupOrDefault<scalar>
    (
        "minDensity",
        1.0
    );

    label minRhoCount = 0;
    label maxRhoCount = 0;

    if (densityClamp)
    {
        // Internal field
        scalarField& rhoIn = rho.primitiveFieldRef();
        forAll(rhoIn, cellI)
        {
            if (rhoIn[cellI] < minDensity)
            {
                rhoIn[cellI] = minDensity;
                minRhoCount++;
            }
            else if (rhoIn[cellI] > maxDensity)
            {
                rhoIn[cellI] = maxDensity;
                maxRhoCount++;
            }
        }

        // Boundary fields
        forAll(rho.boundaryField(), patchI)
        {
            fvPatchScalarField& rhoPatch = rho.boundaryFieldRef()[patchI];
            forAll(rhoPatch, faceI)
            {
                if (rhoPatch[faceI] < minDensity)
                {
                    rhoPatch[faceI] = minDensity;
                    minRhoCount++;
                }
                else if (rhoPatch[faceI] > maxDensity)
                {
                    rhoPatch[faceI] = maxDensity;
                    maxRhoCount++;
                }
            }
        }

        reduce(minRhoCount, sumOp<label>());
        reduce(maxRhoCount, sumOp<label>());

        Info<< "Density clamping stats: " << minRhoCount
            << " cells to minDensity (" << minDensity << "), "
            << maxRhoCount << " cells to maxDensity (" << maxDensity
            << ")" << endl;
    }

    // Update kinetic energy
    K = 0.5*magSqr(U);
    
    Info<< "Pressure solve complete: max p = " << max(p).value() 
        << " Pa, min p = " << min(p).value() << " Pa" << endl;
}
