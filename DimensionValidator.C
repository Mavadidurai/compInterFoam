#include "DimensionValidator.H"
#include <cmath>

namespace Foam
{
// Initialize static dimension sets
const dimensionSet DimensionValidator::dimVelocity(0, 1, -1, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimArea    (0, 2,  0, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimPressure(1,-1, -2, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimTemperature(0, 0, 0, 1, 0, 0, 0);
const dimensionSet DimensionValidator::dimEnergy   (1, 2, -2, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimTempRate (0, 0, -1,1, 0, 0, 0);
const dimensionSet DimensionValidator::dimDensity  (1,-3,  0, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimSpecificHeat(0, 2, -2,-1, 0, 0, 0);
const dimensionSet DimensionValidator::dimThermalCond (1, 1, -3,-1, 0, 0, 0);
const dimensionSet DimensionValidator::dimHeatSource  (1,-1, -3, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimLatentHeat  (0, 2, -2, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimAcceleration(0, 1, -2, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimMomentum    (1, 1, -2, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimSpecificMomentum(0, 1, -2, 0, 0, 0, 0);
const dimensionSet DimensionValidator::dimInverseA    (-1, 3,  1, 0, 0, 0, 0);

void DimensionValidator::validateTEqnTerms
(
    const volScalarField& T,
    const volScalarField& rho,
    const surfaceScalarField& phi,
    const volScalarField& alphaEff,
    const volScalarField& Q_laser,
    const volScalarField& phaseChangeSource
)
{
    // Check dimensions of each term
    checkVolFieldDimensions(T, dimTemperature, "Temperature");
    checkVolFieldDimensions(rho, dimDensity, "Density");
    checkSurfaceFieldDimensions(phi, dimVelocity*dimArea, "Flux");
    checkVolFieldDimensions(alphaEff, dimThermalCond/(dimDensity*dimSpecificHeat), "Thermal Diffusivity");
    checkVolFieldDimensions(Q_laser, dimHeatSource, "Laser Source");
    checkVolFieldDimensions(phaseChangeSource, dimTempRate, "Phase Change Source");

    // Log validation
    Info<< "Temperature equation terms validated successfully" << endl;
}

void DimensionValidator::validateSourceTerms
(
    const volScalarField& Q_laser,
    const volScalarField& phaseChangeSource,
    const dimensionedScalar& Cp,
    const dimensionedScalar& rho,
    const dimensionedScalar& latentHeat
)
{
    // Check source term dimensions
    checkVolFieldDimensions(Q_laser, dimHeatSource, "Laser Source");
    checkVolFieldDimensions(phaseChangeSource, dimTempRate, "Phase Change Source");

    // Check material property dimensions
    if (Cp.dimensions() != dimSpecificHeat)
    {
        FatalErrorInFunction
            << "Invalid specific heat dimensions" << abort(FatalError);
    }

    if (rho.dimensions() != dimDensity)
    {
        FatalErrorInFunction
            << "Invalid density dimensions" << abort(FatalError);
    }

    if (latentHeat.dimensions() != dimLatentHeat)
    {
        FatalErrorInFunction
            << "Invalid latent heat dimensions" << abort(FatalError);
    }

    Info<< "Source terms validated successfully" << endl;
}

void DimensionValidator::checkPhysicalConstraints
(
    volScalarField& T,
    volScalarField& p,
    volScalarField& alpha1
)
{
    // Traditional OpenFOAM: Reference values for reporting only
    const scalar roomTemperature = 300.0;        // K
    const scalar vaporTemperature = 3560.0;      // K - Titanium vaporization
    const scalar atmosphericPressure = 101325.0; // Pa
    const scalar minPhysicalAlpha = 0.0;         // Volume fraction bounds
    const scalar maxPhysicalAlpha = 1.0;
    
    // Traditional validation: Statistical reporting
    scalar minT = min(T).value();
    scalar maxT = max(T).value();
    scalar meanT = gAverage(T);
    
    scalar minP = min(p).value();
    scalar maxP = max(p).value();
    scalar meanP = gAverage(p);
    
    scalar minAlpha = min(alpha1).value();
    scalar maxAlpha = max(alpha1).value();
    scalar meanAlpha = gAverage(alpha1);
    
    // Traditional reporting
    Info<< "Field validation statistics:" << nl
        << "Temperature field:" << nl
        << "  Range: " << minT << " to " << maxT << " K" << nl
        << "  Mean: " << meanT << " K" << nl;
    
    Info<< "Pressure field:" << nl
        << "  Range: " << minP << " to " << maxP << " Pa" << nl
        << "  Mean: " << meanP << " Pa" << nl;
    
    Info<< "Phase fraction field:" << nl
        << "  Range: " << minAlpha << " to " << maxAlpha << nl
        << "  Mean: " << meanAlpha << endl;
    
    // Traditional approach: Only fix truly unphysical values
    label corruptedTempCells = 0;
    label corruptedPressureCells = 0;
    label corruptedAlphaCells = 0;
    
    forAll(T, cellI)
    {
        // Only fix corrupted temperature values (NaN/Inf)
        if (!std::isfinite(T[cellI]))
        {
            T[cellI] = roomTemperature;
            corruptedTempCells++;
        }
        
        // Only fix corrupted pressure values (NaN/Inf)
        if (!std::isfinite(p[cellI]))
        {
            p[cellI] = atmosphericPressure;
            corruptedPressureCells++;
        }
        
        // Only fix corrupted alpha values (NaN/Inf or severely out of bounds)
        if (!std::isfinite(alpha1[cellI]) || alpha1[cellI] < -0.1 || alpha1[cellI] > 1.1)
        {
            alpha1[cellI] = max(min(alpha1[cellI], maxPhysicalAlpha), minPhysicalAlpha);
            corruptedAlphaCells++;
        }
    }
    
    // Reduce across processors
    reduce(corruptedTempCells, sumOp<label>());
    reduce(corruptedPressureCells, sumOp<label>());
    reduce(corruptedAlphaCells, sumOp<label>());
    
    // Traditional reporting of fixes
    if (corruptedTempCells > 0)
    {
        WarningInFunction
            << "Fixed " << corruptedTempCells << " corrupted temperature values" << endl;
    }
    
    if (corruptedPressureCells > 0)
    {
        WarningInFunction
            << "Fixed " << corruptedPressureCells << " corrupted pressure values" << endl;
    }
    
    if (corruptedAlphaCells > 0)
    {
        WarningInFunction
            << "Fixed " << corruptedAlphaCells << " corrupted phase fraction values" << endl;
    }
    
    // Traditional warnings for extreme but physically possible values
    if (maxT > vaporTemperature)
    {
        Info<< "Note: Temperatures above vaporization point detected (" 
            << maxT << " K > " << vaporTemperature << " K)" << nl
            << "This is physically possible for LIFT process" << endl;
    }
    
    if (minT < roomTemperature * 0.9) // 10% below room temp
    {
        Info<< "Note: Low temperatures detected (" << minT << " K)" << nl
            << "Check initial conditions if unexpected" << endl;
    }
    
    if (maxP > 10.0 * atmosphericPressure) // 10 atm
    {
        Info<< "Note: High pressures detected (" << maxP/1e6 << " MPa)" << nl
            << "This may be normal for LIFT recoil pressure" << endl;
    }
    
    // Update boundary conditions after any fixes
    T.correctBoundaryConditions();
    p.correctBoundaryConditions();
    alpha1.correctBoundaryConditions();
}

void DimensionValidator::validateMaterialProperties
(
    const dimensionedScalar& Cp,
    const dimensionedScalar& k,
    const dimensionedScalar& rho,
    const dimensionedScalar& latentHeat
)
{
    if (Cp.value() <= 0 || k.value() <= 0 || rho.value() <= 0 || latentHeat.value() < 0)
    {
        FatalErrorInFunction
            << "Invalid material properties detected:" << nl
            << "Cp = " << Cp.value() << nl
            << "k = " << k.value() << nl
            << "rho = " << rho.value() << nl
            << "latentHeat = " << latentHeat.value()
            << abort(FatalError);
    }

    if (Cp.dimensions() != dimSpecificHeat)
    {
        FatalErrorInFunction
            << "Invalid specific heat dimensions" << abort(FatalError);
    }

    if (k.dimensions() != dimThermalCond)
    {
        FatalErrorInFunction
            << "Invalid thermal conductivity dimensions" << abort(FatalError);
    }

    if (rho.dimensions() != dimDensity)
    {
        FatalErrorInFunction
            << "Invalid density dimensions" << abort(FatalError);
    }

    if (latentHeat.dimensions() != dimLatentHeat)
    {
        FatalErrorInFunction
            << "Invalid latent heat dimensions" << abort(FatalError);
    }
}

bool DimensionValidator::validateMomentumEquation
(
    const fvVectorMatrix& UEqn,
    const volVectorField& HbyA,
    const volScalarField& rAU,
    const surfaceScalarField& surfaceTensionForce
)
{
    bool valid = true;
    
    // Check expected dimensions
    if (UEqn.dimensions() != dimMomentum)
    {
        WarningInFunction
            << "UEqn has unexpected dimensions: " << UEqn.dimensions()
            << ", expected " << dimMomentum << endl;
        valid = false;
    }
    
    if (HbyA.dimensions() != dimVelocity)
    {
        WarningInFunction
            << "HbyA has unexpected dimensions: " << HbyA.dimensions()
            << ", expected " << dimVelocity << endl;
        valid = false;
    }
    
    if (rAU.dimensions() != dimInverseA)
    {
        WarningInFunction
            << "rAU has unexpected dimensions: " << rAU.dimensions()
            << ", expected " << dimInverseA << endl;
        valid = false;
    }
    
    if (surfaceTensionForce.dimensions() != dimPressure)
    {
        WarningInFunction
            << "surfaceTensionForce has unexpected dimensions: " 
            << surfaceTensionForce.dimensions()
            << ", expected " << dimPressure << endl;
        valid = false;
    }
    
    return valid;
}

bool DimensionValidator::validatePressureEquation
(
    const volScalarField& p_rgh,
    const surfaceScalarField& phi,
    const surfaceScalarField& phiHbyA, 
    const surfaceScalarField& rAUf
)
{
    bool valid = true;
    
    // Check expected dimensions
    if (p_rgh.dimensions() != dimPressure)
    {
        WarningInFunction
            << "p_rgh has unexpected dimensions: " << p_rgh.dimensions()
            << ", expected " << dimPressure << endl;
        valid = false;
    }
    
    if (phi.dimensions() != dimVelocity*dimArea)
    {
        WarningInFunction
            << "phi has unexpected dimensions: " << phi.dimensions()
            << ", expected " << dimVelocity*dimArea << endl;
        valid = false;
    }
    
    if (phiHbyA.dimensions() != dimVelocity*dimArea)
    {
        WarningInFunction
            << "phiHbyA has unexpected dimensions: " << phiHbyA.dimensions()
            << ", expected " << dimVelocity*dimArea << endl;
        valid = false;
    }
    
    if (rAUf.dimensions() != dimInverseA)
    {
        WarningInFunction
            << "rAUf has unexpected dimensions: " << rAUf.dimensions()
            << ", expected " << dimInverseA << endl;
        valid = false;
    }
    
    return valid;
}

} // End namespace Foam
