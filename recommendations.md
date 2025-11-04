# Specific Code Locations and Recommendations

## Critical Fixes Required

### 1. ADD PRESSURE GRADIENT LIMITING (HIGH PRIORITY)

**File**: `/home/user/compInterFoam/UEqn.H`
**Lines**: 349-430

**Current Code** (UNCLAMPED):
```cpp
const Foam::scalar gradLimit = solverCoeffs.lookupOrDefault<Foam::scalar>
(
    "maxPressureGradient",
    Foam::GREAT  // ← PROBLEM: Defaults to no limit!
);

bool nonFiniteGrad = false;
bool gradClamped = false;

vectorField& gradInternal = gradPrgh.internalFieldRef();

forAll(gradInternal, celli)
{
    vector& g = gradInternal[celli];
    
    if (gradLimit < Foam::GREAT)  // ← Only clamps if explicitly set!
    {
        const scalar magG = Foam::mag(g);
        if (magG > gradLimit && magG > SMALL)
        {
            g *= gradLimit/magG;
            gradClamped = true;
        }
    }
}
```

**Recommended Fix**:
```cpp
// Set a REASONABLE DEFAULT for pressure gradient limiting
const Foam::scalar gradLimitDefault = 
    solverCoeffs.lookupOrDefault<Foam::scalar>
    (
        "maxPressureGradient",
        1e9  // DEFAULT: 1e9 Pa/m (instead of GREAT!)
    );

// Also check for explicit disable via negative value
const bool clampGradient = (gradLimitDefault > Foam::SMALL);

bool nonFiniteGrad = false;
bool gradClamped = false;

vectorField& gradInternal = gradPrgh.internalFieldRef();

forAll(gradInternal, celli)
{
    vector& g = gradInternal[celli];
    
    // Always check for non-finite values
    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        scalar& val = g[cmpt];
        if (!std::isfinite(val))
        {
            val = 0.0;
            nonFiniteGrad = true;
        }
    }
    
    // Clamp magnitude if enabled
    if (clampGradient)
    {
        const scalar magG = Foam::mag(g);
        if (magG > gradLimitDefault && magG > SMALL)
        {
            g *= gradLimitDefault/magG;
            gradClamped = true;
        }
    }
}
```

**Suggested controlDict Entry**:
```
compInterFoamCoeffs
{
    maxPressureGradient  1e9;     // Pa/m (typical: 1e6 to 1e9)
    // To disable: use -1 or comment out
}
```

---

### 2. BOUND GAS-METAL COUPLING COEFFICIENT (HIGH PRIORITY)

**File**: `/home/user/compInterFoam/TEqn.H`
**Lines**: 121-147

**Current Code**:
```cpp
const dimensionedScalar maxCoeff =
    twoTemperatureDict.lookupOrDefault<dimensionedScalar>
    (
        "maxGasMetalExchangeCoeff",
        maxCoeffDefault  // ← Problem: maxCoeffDefault itself is GREAT!
    );
```

**Recommended Fix**:
```cpp
// Define a physically-reasonable upper bound for coupling coefficient
// Literature suggests values > 1e16 lead to stiff coupling
const scalar maxCoeffDefault_Value = 
    ctrl.lookupOrDefault<scalar>(
        "maxGasMetalExchangeCoeff_default", 
        1e15  // W/m³/K (moderate coupling, avoids stiffness)
    );

const dimensionedScalar maxCoeff =
    twoTemperatureDict.lookupOrDefault<dimensionedScalar>
    (
        "maxGasMetalExchangeCoeff",
        dimensionedScalar
        (
            "maxGasMetalExchangeCoeff",
            interfaceExchange.dimensions(),
            maxCoeffDefault_Value  // ← BOUNDED default
        )
    );

// Warn if configured value exceeds safety threshold
if (maxCoeff.value() > 1e16 && verbose && master)
{
    WarningInFunction
        << "Gas-metal exchange coefficient is very large: " 
        << maxCoeff.value() << " W/m³/K (> 1e16)" << nl
        << "This may cause stiff coupling. Consider reducing or "
        << "enabling thermalLimiterAlphaMin." << endl;
}
```

**Suggested controlDict Entry**:
```
twoTemperatureProperties
{
    maxGasMetalExchangeCoeff  1e15;  // W/m³/K (prevents stiff coupling)
    // For stiff coupling: use values 1e16-1e17
    // For loose coupling: use values 1e12-1e14
}
```

---

### 3. ADD VELOCITY LIMITING BEFORE MOMENTUM SOLVE (MEDIUM PRIORITY)

**File**: `/home/user/compInterFoam/UEqn.H`
**Lines**: 56-82

**Current Code** (ONLY CHECKS, doesn't prevent):
```cpp
const Foam::scalar kineticEnergyCeiling =
    solverCoeffs.lookupOrDefault<Foam::scalar>
    (
        "maxKineticEnergyDensity",
        1e12
    );

if (kineticEnergyCeiling > 0 && kineticEnergyCeiling < Foam::GREAT)
{
    tmp<volScalarField> tKineticEnergyDensity(0.5*rho*magSqr(U));
    const volScalarField& kineticEnergyDensity = tKineticEnergyDensity();
    const dimensionedScalar maxKineticEnergyDensity = gMax(kineticEnergyDensity);
    const scalar maxKineticEnergyDensityValue = maxKineticEnergyDensity.value();

    if (!std::isfinite(maxKineticEnergyDensityValue) || 
        maxKineticEnergyDensityValue > kineticEnergyCeiling)
    {
        FatalErrorInFunction  // ← Fails here instead of preventing
            << "Maximum kinetic energy density " << ... << exit(FatalError);
    }
}
```

**Recommended Fix** (Prevent rather than fail):
```cpp
const Foam::scalar kineticEnergyCeiling =
    solverCoeffs.lookupOrDefault<Foam::scalar>
    (
        "maxKineticEnergyDensity",
        1e12
    );

const scalar kineticEnergyWarningThreshold =
    solverCoeffs.lookupOrDefault<scalar>(
        "kineticEnergyWarningThreshold",
        0.8 * kineticEnergyCeiling  // Warn at 80% of ceiling
    );

if (kineticEnergyCeiling > 0 && kineticEnergyCeiling < Foam::GREAT)
{
    tmp<volScalarField> tKineticEnergyDensity(0.5*rho*magSqr(U));
    const volScalarField& kineticEnergyDensity = tKineticEnergyDensity();
    const dimensionedScalar maxKineticEnergyDensity = gMax(kineticEnergyDensity);
    const scalar maxKineticEnergyDensityValue = maxKineticEnergyDensity.value();

    // Check for non-finite values (early warning)
    if (!std::isfinite(maxKineticEnergyDensityValue))
    {
        WarningInFunction
            << "Non-finite kinetic energy density detected (NaN or Inf)!" << nl
            << "This indicates numerical instability upstream." << nl
            << "Check pressure gradients and temperature coupling." << endl;
    }
    
    // Warn at threshold
    if (maxKineticEnergyDensityValue > kineticEnergyWarningThreshold)
    {
        WarningInFunction
            << "Kinetic energy density " << maxKineticEnergyDensityValue 
            << " J/m³ exceeds warning threshold " 
            << kineticEnergyWarningThreshold << " J/m³" << nl
            << "Approaching limit of " << kineticEnergyCeiling << " J/m³" << endl;
    }
    
    // HARD LIMIT: Fail to prevent divergence
    if (maxKineticEnergyDensityValue > kineticEnergyCeiling)
    {
        FatalErrorInFunction
            << "Maximum kinetic energy density "
            << maxKineticEnergyDensityValue
            << " J/m^3 exceeds limit "
            << kineticEnergyCeiling
            << " J/m^3. Velocity field has become unphysical." << nl
            << "Root causes typically:" << nl
            << "  1. Unclamped pressure gradients (set maxPressureGradient)" << nl
            << "  2. Stiff thermal coupling (reduce maxGasMetalExchangeCoeff)" << nl
            << "  3. Temperature-driven density discontinuities" << nl
            << "  4. Time step too large" << exit(FatalError);
    }
}
```

---

### 4. MOVE PRESSURE CLAMPING BEFORE GRADIENT COMPUTATION (MEDIUM PRIORITY)

**File**: `/home/user/compInterFoam/pEqn.H`
**Lines**: 72-334

**Current Sequence** (WRONG ORDER):
```cpp
// Line 72-79: Pressure equation solved
fvScalarMatrix pEqn(fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA));
pEqn.setReference(pRefCell, pRefValue);
pEqn.solve();

// ... other code ...

// Line 290-334: PRESSURE CLAMPED (too late!)
if (!pressureClamp && (pressureNonFinite || pressureOutOfBounds))
{
    // Clamp p field
}
```

**Problem**: Clamping after solve but then using clamped pressure in next iteration's gradient.

**Recommended Sequence**:
```cpp
// After pressure solve, immediately clamp
fvScalarMatrix pEqn(fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA));
pEqn.setReference(pRefCell, pRefValue);
pEqn.solve();

// NEW: Immediate clamping after solve
const Foam::dictionary& clampCoeffs = compInterFoamCoeffsDict(mesh);
const Foam::Switch pressureClamp =
    clampCoeffs.lookupOrDefault<Foam::Switch>("pressureClamp", false);

if (pimple.finalNonOrthogonalIter())
{
    if (pressureClamp)
    {
        Foam::dimensionedScalar minPressure = ...;
        Foam::dimensionedScalar maxPressure = ...;
        
        // Clamp p_rgh immediately
        scalarField& pInternal = p_rgh.primitiveFieldRef();
        forAll(pInternal, celli)
        {
            pInternal[celli] = Foam::min(
                Foam::max(pInternal[celli], minPressure.value()),
                maxPressure.value()
            );
        }
        p_rgh.correctBoundaryConditions();
    }
    
    // NOW compute fluxes with clamped pressure
    phi = phiHbyA + pEqn.flux();
    U = HbyA + rAU*fvc::reconstruct((phig + pEqn.flux())/rAUf);
}

// ... rest of code ...

// Line 290: Remove the emergency clamp (no longer needed)
```

---

### 5. ADD DENSITY-CONSISTENT PRESSURE SOURCE (ADVANCED)

**File**: `/home/user/compInterFoam/pEqn.H`
**Lines**: 72-79

**Current Code** (Missing density source):
```cpp
fvScalarMatrix pEqn
(
    fvm::laplacian(rAUf, p_rgh)
 == fvc::div(phiHbyA)
);
```

**Theory**: When density changes with temperature, pressure equation should include ∂ρ/∂t term.

**Enhanced Version**:
```cpp
// Account for temperature-driven density changes
tmp<volScalarField> tRhoDt = fvc::ddt(rho);
const volScalarField& rhoDt = tRhoDt();

// Create enhanced pressure equation
fvScalarMatrix pEqn
(
    fvm::laplacian(rAUf, p_rgh)
 == fvc::div(phiHbyA)
  + fvc::div(rhoDt)  // ← NEW: Density time derivative source
);
```

**Caveat**: This requires careful tuning and may change convergence behavior.

---

### 6. LOOSEN THERMAL COUPLING TIME STEP RESTRICTION (MEDIUM PRIORITY)

**File**: `/home/user/compInterFoam/TEqn.H`
**Lines**: 946-1010

**Current Code**:
```cpp
if (maxThermalCourant > SMALL && maxThermalCourant < 1.0 - SMALL)
{
    // Compute thermal Courant number
    tmp<volScalarField> tThermalCo
    (
        (cappedInterfaceExchange*runTime.deltaT())
      / Foam::max(CvolLimited, minCvol)
    );
    
    const dimensionedScalar maxThermalCourantDim("maxThermalCourant", dimless, maxThermalCourant);
    
    tmp<volScalarField> tLimiter
    (
        maxThermalCourantDim/(thermalCo + epsCo)
    );
    
    // THIS LIMITER CAN BECOME VERY SMALL
    limiter = Foam::min(limiter, oneDim);
}
```

**Recommended Fix**:
```cpp
// Add a floor to prevent complete elimination of coupling
const scalar couplingSafetyFactor =
    ctrl.lookupOrDefault<scalar>(
        "couplingSafetyFactor",
        0.1  // ← Keep at least 10% of original coupling
    );

if (maxThermalCourant > SMALL && maxThermalCourant < 1.0 - SMALL)
{
    tmp<volScalarField> tThermalCo(
        (cappedInterfaceExchange*runTime.deltaT())
      / Foam::max(CvolLimited, minCvol)
    );
    
    tmp<volScalarField> tLimiter(
        maxThermalCourantDim/(thermalCo + epsCo)
    );
    
    volScalarField& limiter = tLimiter.ref();
    
    // Ensure coupling is not completely eliminated
    limiter = Foam::max(
        limiter,
        dimensionedScalar("couplingSafetyFactor", dimless, couplingSafetyFactor)
    );
    
    limiter = Foam::min(limiter, oneDim);
    limiter = Foam::max(limiter, dimensionedScalar("zero", dimless, 0.0));
    
    cappedInterfaceExchange *= limiter;
}
```

---

### 7. MONITOR AND DIAGNOSE COUPLING ISSUES (DIAGNOSTIC)

**File**: `/home/user/compInterFoam/compInterFoam.C`
**After Line**: 960 (after recoil calculation)

**Add Diagnostic Output**:
```cpp
// After recoil pressure calculation
if (useAdvancedCapturing && pInterfaceCapturing.valid())
{
    pInterfaceCapturing->calculateRecoilPressure();
    
    // NEW: Diagnostic output
    if (verbose && master && runTime.timeIndex() % 10 == 0)
    {
        const volScalarField& recoilP = pInterfaceCapturing->recoilPressure();
        
        tmp<volScalarField> tGradP(fvc::grad(p));
        const volScalarField& gradP = tGradP();
        
        const dimensionedScalar maxGradP = gMax(mag(gradP));
        const dimensionedScalar maxVelFromP = sqrt(2.0*gMax(mag(p))/gMax(rho));
        const dimensionedScalar maxKE = fvc::domainIntegrate(0.5*rho*magSqr(U));
        
        Info<< "    Numerical diagnostics:" << nl
            << "      max(|∇p|) = " << maxGradP.value() << " Pa/m" << nl
            << "      Max velocity from pressure alone: " << maxVelFromP.value() << " m/s" << nl
            << "      Kinetic energy integral: " << maxKE.value() << " J" << nl
            << "      KE density max: " << (maxKE.value()/gSum(mesh.V())).value() << " J/m³" << endl;
    }
}
```

---

## Recommended controlDict Entries

Add these to your system/controlDict for stability:

```
compInterFoamCoeffs
{
    // PRESSURE GRADIENT LIMITING (CRITICAL)
    maxPressureGradient     1e9;        // Pa/m (prevent unbounded gradients)
    maxVelocity             500;        // m/s (hard limit on velocity)
    
    // VELOCITY FIELD LIMITING (IMPORTANT)
    maxUEqnVelocity         500;        // m/s (before momentum solve)
    minUEqnDiag             1e-9;       // Diagonal floor
    
    // KINETIC ENERGY LIMITS (IMPORTANT)
    maxKineticEnergyDensity 1e12;       // J/m³
    
    // PRESSURE LIMITS (IMPORTANT)
    pressureClamp           true;       // Enable immediate clamping
    minPressure            -1e7;        // Pa
    maxPressure             1e9;        // Pa
    
    // VELOCITY DIAGNOSTICS (DIAGNOSTIC)
    velocityAlphaThreshold  0.01;       // Ignore <1% metal cells for diagnostics
}

twoTemperatureProperties
{
    // GAS-METAL COUPLING (CRITICAL)
    maxGasMetalExchangeCoeff   1e15;    // W/m³/K (bounded, no stiff coupling)
    minGasMetalExchangeCoeff   1e10;    // W/m³/K (lower bound)
    
    // INTERFACE DETECTION
    metalFractionFloor         1e-6;    // Minimum metal volume fraction
    metalAmbientBlendWidth     1e-3;    // Smooth transition width
    
    // TEMPERATURE CLAMPING (OPTIONAL)
    enableTClamp               false;   // Set to true if temperatures exceed bounds
    minTe                      100;     // K (electron temp floor)
    maxTe                      50000;   // K (electron temp ceiling)
    minTl                      100;     // K
    maxTl                      8000;    // K (typical for metals near Tvap)
}

PIMPLE
{
    nOuterCorrectors        3;
    nCorrectors             2;
    nNonOrthogonalCorrectors 1;
    maxCo                   0.5;        // Reduce from default 1.0
    maxAlphaCo              0.2;
}
```

---

## Testing Sequence

1. **Add maxPressureGradient first** (default 1e9)
2. **Verify with zero recoil** - should still be stable
3. **Reduce maxGasMetalExchangeCoeff** to 1e15 if still unstable
4. **Enable pressureClamp = true** 
5. **Reduce maxCo** to 0.5
6. **Monitor diagnostics** with verbose = true

