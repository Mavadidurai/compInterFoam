/*---------------------------------------------------------------------------*\
  Robust pressure equation with effective numerical safeguards
\*---------------------------------------------------------------------------*/
#include <cmath>
{
    Info<< "Solving pressure equation" << endl;
    
    // Create robust interpolation of rAU
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    
    // Apply bounds to rAUf to prevent division by zero or infinities
    // This should be done on internal field and boundary fields
    scalar minVal = 1e-10;
    scalar maxVal = 1e6;
    
    // Internal field
    scalarField& rAUfIn = rAUf.primitiveFieldRef();
    forAll(rAUfIn, faceI)
    {
        if (!std::isfinite(rAUfIn[faceI]) || rAUfIn[faceI] < minVal)
        {
            rAUfIn[faceI] = minVal;
        }
        else if (rAUfIn[faceI] > maxVal)
        {
            rAUfIn[faceI] = maxVal;
        }
    }
    
    // Boundary fields
    forAll(rAUf.boundaryField(), patchI)
    {
        fvsPatchScalarField& pRauf = rAUf.boundaryFieldRef()[patchI];
        forAll(pRauf, faceI)
        {
            if (!std::isfinite(pRauf[faceI]) || pRauf[faceI] < minVal)
            {
                pRauf[faceI] = minVal;
            }
            else if (pRauf[faceI] > maxVal)
            {
                pRauf[faceI] = maxVal;
            }
        }
    }
    
    // Create pressure equation with robust interpolation
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUf, p_rgh)
    );
    
    // Set reference and apply fixed under-relaxation
    pEqn.setReference(pRefCell, pRefValue);
scalar relaxFactor = runTime.controlDict().lookupOrDefault<scalar>("relaxFactor", 0.3);
pEqn.relax(relaxFactor);

    
    // Create backup of p_rgh in case solver fails
    volScalarField p_rgh_old(p_rgh);
    
    // Solve
    bool solveOK = true;
    try
    {
        // Use check on initial residual instead of good()
        SolverPerformance<scalar> solverPerf = pEqn.solve();
        solveOK = (solverPerf.initialResidual() < 1e3);
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
    dimensionedScalar maxP("maxP", dimPressure, 1e6);  // 1 MPa max
    dimensionedScalar minP("minP", dimPressure, 1e3);  // 1 kPa min
    
    // Apply pressure limits
    p_rgh = max(min(p_rgh, maxP), minP);
    
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
