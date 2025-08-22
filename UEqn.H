/*---------------------------------------------------------------------------*\
  Enhanced momentum equation with robust numerical safeguards for LIFT
\*---------------------------------------------------------------------------*/

#include <cmath>

Info<< "Solving momentum equation" << endl;

// Solver controls dictionary
const dictionary& solverControls = mesh.solutionDict().subOrEmptyDict("compInterFoamCoeffs");

// Enhanced constants for better stability
const scalar maxReasonableVelocity = solverControls.lookupOrDefault<scalar>(
    "maxReasonableVelocity", 100.0);
const scalar defaultMetalDensity = solverControls.lookupOrDefault<scalar>(
    "defaultMetalDensity", 4500.0);
const scalar defaultAirDensity = solverControls.lookupOrDefault<scalar>(
    "defaultAirDensity", 1.2);
    const scalar maxVelocity = solverControls.lookupOrDefault<scalar>(
    "maxVelocity", 100.0);
// Dynamic momentum relaxation control
static scalar relaxU = 0.3;
const scalar targetResU = solverControls.lookupOrDefault<scalar>(
    "targetResU", 0.1);

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
       // Use phase-appropriate density based on alpha1
        if (alpha1[cellI] > 0.5)
        {
            // Predominantly phase 1 (metal)
            if (std::isfinite(rho1[cellI]) && rho1[cellI] > 0)
            {
                rho[cellI] = rho1[cellI];
            }
            else
            {
                rho[cellI] = defaultMetalDensity;
            }
        }
        else
        {
            // Predominantly phase 2 (air)
            if (std::isfinite(rho2[cellI]) && rho2[cellI] > 0)
            {
                rho[cellI] = rho2[cellI];
            }
            else
            {
                rho[cellI] = defaultAirDensity;
            }
        }
        corruptedRhoCells++;
    }
}

reduce(corruptedRhoCells, sumOp<label>());
if (corruptedRhoCells > 0)
{
    Info<< "  Fixed " << corruptedRhoCells << " corrupted density values" << endl;
}
// Update thermophysical properties and phase-change source
mixture.correct();
// Refresh only the recoil pressure using latest phase-change source.
// The alpha field has already been corrected in the alpha subcycle,
// so avoid a second full interface-capturing solve here.
if (useAdvancedCapturing && pInterfaceCapturing.valid())
{
    pInterfaceCapturing->calculateRecoilPressure();
}
// Create momentum equation
fvVectorMatrix UEqn(fvm::ddt(rho, U) + fvm::div(rhoPhi, U));

UEqn += turbulence.divDevRhoReff(U);

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
        
        // Add recoil pressure-induced momentum source
        UEqn += (recoilPressure / rho) * fvc::grad(alpha1);
        
        Info<< "Recoil pressure effects:" << nl
            << "  Maximum: " << maxRecoil << " Pa" << endl;
    }
}

// Apply adaptive momentum relaxation
UEqn.relax(relaxU);

// Obtain diagonal coefficients and compute reciprocal with optional clamping
tmp<volScalarField> tA = UEqn.A();


const dictionary& ctrl = runTime.controlDict();

const Switch enableRAUClamp = solverControls.lookupOrDefault<Switch>
(
    "enableRAUClamp",
    ctrl.lookupOrDefault<Switch>("enableRAUClamp", true)
);

const scalar minRau = solverControls.lookupOrDefault<scalar>
(
    "minRAU",
    ctrl.lookupOrDefault<scalar>("minRAU", 1e-10)
);

const scalar maxRau = solverControls.lookupOrDefault<scalar>
(
    "maxRAU",
    ctrl.lookupOrDefault<scalar>("maxRAU", 1e6)
);

// Create rAU field
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
    dimensionedScalar("rAU", dimTime/dimDensity, 0),
    extrapolatedCalculatedFvPatchScalarField::typeName
);

label minRauCells = 0;
label maxRauCells = 0;

forAll(rAU, cellI)
{
    scalar Acoeff = tA()[cellI];

    if (!std::isfinite(Acoeff) || mag(Acoeff) <= SMALL)
    {
        rAU[cellI] = minRau;
        if (enableRAUClamp) minRauCells++;
    }
    else
    {
        scalar value = 1.0 / Acoeff;
        if (enableRAUClamp)
        {
            if (!std::isfinite(value) || value < minRau)
            {
                rAU[cellI] = minRau;
                minRauCells++;
            }
            else if (value > maxRau)
            {
                rAU[cellI] = maxRau;
                maxRauCells++;
            }
            else
            {
                rAU[cellI] = value;
            }
        }
        else
        {
            rAU[cellI] = value;
        }
    }
}

reduce(minRauCells, sumOp<label>());
reduce(maxRauCells, sumOp<label>());
if (enableRAUClamp && (minRauCells > 0 || maxRauCells > 0))
{
    Info<< "rAU clamping stats: " << minRauCells
        << " cells to minRAU (" << minRau << "), "
        << maxRauCells << " cells to maxRAU (" << maxRau << ")" << endl;
}
else if (!enableRAUClamp &&
         (min(rAU).value() < minRau || max(rAU).value() > maxRau))
{
    Info<< "rAU outside expected range [" << minRau << ", " << maxRau
        << "]: min=" << min(rAU).value()
        << " max=" << max(rAU).value() << endl;
}

rAU.correctBoundaryConditions();


// Momentum predictor
if (pimple.momentumPredictor())
{
    SolverPerformance<vector> solverPerf = solve(UEqn == -fvc::grad(p_rgh));
    scalar residual = mag(solverPerf.initialResidual());
    relaxU = max(min(relaxU * targetResU/(residual + SMALL), 0.7), 0.1);
    
    // Post-solve validation
    scalar maxVelAfterSolve = 0.0;
    label invalidCells = 0;
    label limitedCells = 0;

    
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
            // Bound extreme velocities
            if (velMag > maxVelocity)
            {
                U[cellI] *= (maxVelocity/velMag);
                limitedCells++;
            }
        }
    }

    maxVelAfterSolve = gMax(mag(U)().internalField());
    reduce(invalidCells, sumOp<label>());
    reduce(limitedCells, sumOp<label>());

    if (invalidCells > 0)
    {
        Info<< "  Fixed " << invalidCells << " corrupted velocity values after solve" << endl;
    }

    if (limitedCells > 0)
    {
        Info<< "  Limited " << limitedCells << " velocity cells exceeding "
            << maxVelocity << " m/s" << endl;
    }

    Info<< "Momentum predictor:" << nl
        << "  Max velocity after limiting: " << maxVelAfterSolve << " m/s" << endl;
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
