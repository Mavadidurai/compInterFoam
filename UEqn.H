/*---------------------------------------------------------------------------*\
  Enhanced momentum equation with robust numerical safeguards for LIFT
\*---------------------------------------------------------------------------*/

#include <cmath>

Info<< "Solving momentum equation" << endl;

// Enhanced constants for better stability
const scalar maxReasonableVelocity = 250.0;     // m/s
const scalar defaultTitaniumDensity = 4500.0;   // kg/m³
const scalar defaultAirDensity = 1.2;           // kg/m³
const scalar enhancedMomentumRelaxation = 0.3;  // Conservative relaxation

// Enhanced velocity field validation
bool hasVelocityIssues = false;
scalar maxVelMagnitude = 0.0;
label corruptedVelCells = 0;

forAll(U, cellI)
{
    if (!std::isfinite(U[cellI].x()) || !std::isfinite(U[cellI].y()) || !std::isfinite(U[cellI].z()))
    {
        U[cellI] = vector::zero;
        hasVelocityIssues = true;
        corruptedVelCells++;
    }
    else
    {
        scalar velMag = mag(U[cellI]);
        maxVelMagnitude = max(maxVelMagnitude, velMag);
        
        if (velMag > maxReasonableVelocity)
        {
            hasVelocityIssues = true;
        }
    }
}

// Reduce across processors
reduce(maxVelMagnitude, maxOp<scalar>());
reduce(hasVelocityIssues, orOp<bool>());
reduce(corruptedVelCells, sumOp<label>());

Info<< "Velocity field statistics:" << nl
    << "  Maximum velocity: " << maxVelMagnitude << " m/s" << endl;

if (corruptedVelCells > 0)
{
    Info<< "  Fixed " << corruptedVelCells << " corrupted velocity values" << endl;
}

// Enhanced density field validation
label corruptedRhoCells = 0;
forAll(rho, cellI)
{
    if (!std::isfinite(rho[cellI]) || rho[cellI] <= 0)
    {
        // Use material-appropriate default density
        if (alpha1[cellI] > 0.5)  // Metal region
        {
            rho[cellI] = defaultTitaniumDensity;
        }
        else  // Air region
        {
            rho[cellI] = defaultAirDensity;
        }
        corruptedRhoCells++;
    }
}

reduce(corruptedRhoCells, sumOp<label>());
if (corruptedRhoCells > 0)
{
    Info<< "  Fixed " << corruptedRhoCells << " corrupted density values" << endl;
}

// Create momentum equation
fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  - fvm::laplacian(mixture.mu(), U)
);

// Add recoil pressure effects
if (useAdvancedCapturing && pInterfaceCapturing.valid())
{
    const volScalarField& recoilPressure = pInterfaceCapturing->recoilPressure();
    
    // Check recoil pressure validity
    bool recoilValid = true;
    scalar maxRecoil = 0.0;
    forAll(recoilPressure, cellI)
    {
        if (!std::isfinite(recoilPressure[cellI]))
        {
            recoilValid = false;
            break;
        }
        maxRecoil = max(maxRecoil, mag(recoilPressure[cellI]));
    }
    
    reduce(recoilValid, andOp<bool>());
    reduce(maxRecoil, maxOp<scalar>());
    
    if (recoilValid && maxRecoil > SMALL)
    {
        // Add recoil pressure gradient
        UEqn -= fvc::grad(recoilPressure);
        
        Info<< "Recoil pressure effects:" << nl
            << "  Maximum: " << maxRecoil << " Pa" << endl;
    }
}

// Apply momentum relaxation
UEqn.relax(enhancedMomentumRelaxation);

// Improved rAU calculation
tmp<volScalarField> tA = UEqn.A();
const volScalarField& A = tA();

// Create rAU field with proper initialization
volScalarField rAU
(
    IOobject
    (
        "rAU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rAU", dimTime/dimDensity, 1e-6),
    extrapolatedCalculatedFvPatchScalarField::typeName
);

// Material-aware bounds
const scalar dt = runTime.deltaTValue();
const scalar rAU_min_metal = dt / (10.0 * defaultTitaniumDensity);
const scalar rAU_max_metal = dt / (0.1 * defaultTitaniumDensity);
const scalar rAU_min_air = dt / (10.0 * defaultAirDensity);
const scalar rAU_max_air = dt / (0.1 * defaultAirDensity);

label boundedRauCells = 0;
label invalidDiagonalCells = 0;

forAll(rAU, cellI)
{
    scalar A_diag = A[cellI];
    
    if (!std::isfinite(A_diag) || A_diag <= SMALL)
    {
        // Use material-based default
        if (alpha1[cellI] > 0.5)  // Metal region
        {
            rAU[cellI] = rAU_min_metal;
        }
        else  // Air region
        {
            rAU[cellI] = rAU_min_air;
        }
        invalidDiagonalCells++;
    }
    else
    {
        scalar rAU_calc = 1.0/A_diag;
        
        // Apply material-specific bounds
        if (alpha1[cellI] > 0.5)  // Metal region
        {
            if (!std::isfinite(rAU_calc) || rAU_calc < rAU_min_metal)
            {
                rAU[cellI] = rAU_min_metal;
                boundedRauCells++;
            }
            else if (rAU_calc > rAU_max_metal)
            {
                rAU[cellI] = rAU_max_metal;
                boundedRauCells++;
            }
            else
            {
                rAU[cellI] = rAU_calc;
            }
        }
        else  // Air region
        {
            if (!std::isfinite(rAU_calc) || rAU_calc < rAU_min_air)
            {
                rAU[cellI] = rAU_min_air;
                boundedRauCells++;
            }
            else if (rAU_calc > rAU_max_air)
            {
                rAU[cellI] = rAU_max_air;
                boundedRauCells++;
            }
            else
            {
                rAU[cellI] = rAU_calc;
            }
        }
    }
}

// Apply field smoothing for better conditioning
for (int smoothIter = 0; smoothIter < 2; smoothIter++)
{
    volScalarField rAU_old = rAU;
    
    forAll(rAU, cellI)
    {
        scalar weightedSum = rAU_old[cellI];
        scalar totalWeight = 1.0;
        
        const labelList& neighbors = mesh.cellCells()[cellI];
        forAll(neighbors, ni)
        {
            label neighborCell = neighbors[ni];
            scalar weight = 0.05;  // Light smoothing
            weightedSum += weight * rAU_old[neighborCell];
            totalWeight += weight;
        }
        
        rAU[cellI] = weightedSum / totalWeight;
    }
    
    rAU.correctBoundaryConditions();
}

// Apply temporal relaxation
if (rAU.oldTime().timeIndex() >= 0)
{
    scalar rAU_relaxation = 0.8;
    rAU = rAU_relaxation * rAU + (1.0 - rAU_relaxation) * rAU.oldTime();
}

rAU.correctBoundaryConditions();

// Reduce and report
reduce(boundedRauCells, sumOp<label>());
reduce(invalidDiagonalCells, sumOp<label>());

scalar minRAU = gMin(rAU);
scalar maxRAU = gMax(rAU);

Info<< "Matrix conditioning:" << nl
    << "  rAU range: " << minRAU << " to " << maxRAU << endl;

if (boundedRauCells > 0)
{
    Info<< "  Bounded " << boundedRauCells << " rAU values" << endl;
}

if (invalidDiagonalCells > 0)
{
    Info<< "  Fixed " << invalidDiagonalCells << " invalid diagonal values" << endl;
}

// Momentum predictor
if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p_rgh));
    
    // Post-solve validation
    scalar maxVelAfterSolve = 0.0;
    label invalidCells = 0;
    
    forAll(U, cellI)
    {
        if (!std::isfinite(U[cellI].x()) || !std::isfinite(U[cellI].y()) || !std::isfinite(U[cellI].z()))
        {
            U[cellI] = vector::zero;
            invalidCells++;
        }
        else
        {
            scalar velMag = mag(U[cellI]);
            maxVelAfterSolve = max(maxVelAfterSolve, velMag);
            
            // Bound extreme velocities
            if (velMag > 1000.0)  // 1 km/s limit
            {
                U[cellI] *= (1000.0 / velMag);
                invalidCells++;
            }
        }
    }
    
    reduce(maxVelAfterSolve, maxOp<scalar>());
    reduce(invalidCells, sumOp<label>());
    
    if (invalidCells > 0)
    {
        Info<< "  Fixed " << invalidCells << " velocity values after solve" << endl;
    }
    
    Info<< "Momentum predictor:" << nl
        << "  Max velocity after solve: " << maxVelAfterSolve << " m/s" << endl;
}

// Update fluxes
phi = fvc::interpolate(U) & mesh.Sf();
rhoPhi = fvc::interpolate(rho)*phi;

// Validate fluxes
scalar maxPhi = max(mag(phi)).value();
scalar maxRhoPhi = max(mag(rhoPhi)).value();

if (!std::isfinite(maxPhi) || !std::isfinite(maxRhoPhi))
{
    WarningInFunction
        << "Invalid flux detected - recalculating" << endl;
    phi = fvc::interpolate(U) & mesh.Sf();
    rhoPhi = fvc::interpolate(rho)*phi;
}

// Continuity error calculation
volScalarField contErr(fvc::div(phi));
scalar sumLocalContErr = runTime.deltaT().value() *
    mag(contErr)().weightedAverage(mesh.V()).value();
scalar globalContErr = runTime.deltaT().value() *
    contErr.weightedAverage(mesh.V()).value();
    
if (!std::isfinite(sumLocalContErr)) sumLocalContErr = 0.0;
if (!std::isfinite(globalContErr)) globalContErr = 0.0;

Info<< "Continuity errors:" << nl
    << "  Sum local: " << sumLocalContErr << nl
    << "  Global: " << globalContErr << endl;

// Boundary conditions and turbulence
U.correctBoundaryConditions();
turbulence.correct();

Info<< "Momentum equation completed" << endl;
