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
    
    // Create pressure equation with robust interpolation
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUf, p_rgh)
    );
    
    // Set reference and apply adaptive under-relaxation
    pEqn.setReference(pRefCell, pRefValue);
    static scalar relaxP
    (
        mesh.solutionDict().lookupOrDefault<scalar>("relaxP", 0.3)
    );
    const scalar targetResP = mesh.solutionDict().lookupOrDefault<scalar>(
        "targetResP", 0.1);
    pEqn.relax(relaxP);
    
    // Create backup of p_rgh in case solver fails
    volScalarField p_rgh_old(p_rgh);
    
    // Solve
    bool solveOK = true;
    try
    {
        // Use check on initial residual instead of good()
        SolverPerformance<scalar> solverPerf = pEqn.solve();
        scalar residual = solverPerf.initialResidual();
        Info<< "pEqn initial residual: " << residual << endl;
        solveOK = (residual < 1e3);
        relaxP = max(min(relaxP * targetResP/(residual + SMALL), 0.7), 0.1);
        Info<< "  Updated relaxP: " << relaxP << endl;

    }
    catch(...)
    {
        solveOK = false;
        Warning << "Exception caught in pressure solver - reverting to previous solution" << endl;
        p_rgh = p_rgh_old;
    }
    
    if (!solveOK)
    {
        Warning << "Pressure solve appears unstable - reverting to previous solution" << endl;
        p_rgh = p_rgh_old;
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

    // Apply pressure limits only when enabled
    if (pressureClamp)
    {
        dimensionedScalar maxP("maxP", dimPressure, maxPressure);
        dimensionedScalar minP("minP", dimPressure, minPressure);
        p_rgh = max(min(p_rgh, maxP), minP);
    }
    
    // Update pressure field
    p = p_rgh + rho*gh;
    
    // Ensure boundary conditions
    p_rgh.correctBoundaryConditions();
    p.correctBoundaryConditions();
    
    // Update density
    rho = alpha1*rho1 + alpha2*rho2;
    
    // Apply limits to density for safety
    dimensionedScalar minRho("minRho", dimDensity, 1.0);
    dimensionedScalar maxRho("maxRho", dimDensity, 10000.0);
    rho = max(min(rho, maxRho), minRho);
    
    // Update kinetic energy
    K = 0.5*magSqr(U);
    
    Info<< "Pressure solve complete: max p = " << max(p).value() 
        << " Pa, min p = " << min(p).value() << " Pa" << endl;
}
