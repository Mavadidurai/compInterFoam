# Extreme Pressure and Kinetic Energy Overflow Analysis
## compInterFoam Numerical Stability Issues

---

## Executive Summary

The extreme pressure (min(p) = -8.8e+13 Pa) and kinetic energy overflow occur due to **unconstrained coupling between the two-temperature model and momentum equations**. With zero recoil pressure, the simulation still experiences numerical instability through:

1. **Stiff thermal coupling** driving uncontrolled gas temperature gradients
2. **Pressure gradient forces** without adequate stabilization
3. **Discontinuous pressure post-correction** causing gradient spikes
4. **Missing material property stabilization** at the metal/gas interface

---

## 1. KINETIC ENERGY CHECK (UEqn.H, Lines 56-82)

### Location and Trigger
```cpp
// File: /home/user/compInterFoam/UEqn.H
// Lines 49-82: Kinetic energy ceiling check

const Foam::scalar kineticEnergyCeiling =
    solverCoeffs.lookupOrDefault<Foam::scalar>
    (
        "maxKineticEnergyDensity",
        1e12  // DEFAULT: 1e12 J/m³
    );

if (kineticEnergyCeiling > 0 && kineticEnergyCeiling < Foam::GREAT)
{
    // 0.5*rho*|U|^2 has SI units of J/m^3
    tmp<volScalarField> tKineticEnergyDensity(0.5*rho*magSqr(U));
    const volScalarField& kineticEnergyDensity = tKineticEnergyDensity();
    const dimensionedScalar maxKineticEnergyDensity = gMax(kineticEnergyDensity);
    const scalar maxKineticEnergyDensityValue = maxKineticEnergyDensity.value();

    if (!std::isfinite(maxKineticEnergyDensityValue) || 
        maxKineticEnergyDensityValue > kineticEnergyCeiling)
    {
        FatalErrorInFunction  // LINE 74: FATAL ERROR OCCURS HERE
            << "Maximum kinetic energy density " 
            << maxKineticEnergyDensityValue << " J/m^3 exceeds limit "
            << kineticEnergyCeiling << " J/m^3..."
            << exit(FatalError);
    }
}
```

### Root Cause
The check occurs **BEFORE momentum solver** and detects velocities that are already unphysical. 
The velocity field becomes uncontrolled due to:
- Large pressure gradients from temperature-induced density variations
- Insufficient velocity limiting/regularization
- Unstabilized momentum source terms

---

## 2. PRESSURE EQUATION & CLAMPING (pEqn.H, Lines 72-334)

### Pressure Equation Assembly (Lines 72-79)
```cpp
// File: /home/user/compInterFoam/pEqn.H
// Lines 72-79

fvScalarMatrix pEqn
(
    fvm::laplacian(rAUf, p_rgh)
 == fvc::div(phiHbyA)
);

pEqn.setReference(pRefCell, pRefValue);
pEqn.solve();
```

**Problem**: The pressure equation is solved WITHOUT any source terms related to temperature or phase change. 
This creates a mismatch with the momentum equation where temperature gradients drive density variations.

### Pressure Clamping (Lines 290-334)

**Location**: Emergency clamp applied AFTER successful pressure solve
```cpp
// Lines 284-298: Check if pressure is out of bounds
const bool pressureNonFinite =
    (!std::isfinite(observedMinP.value()) || 
     !std::isfinite(observedMaxP.value()));

const bool pressureOutOfBounds =
    (observedMinP.value() < minPressure.value() || 
     observedMaxP.value() > maxPressure.value());

if (!pressureClamp && (pressureNonFinite || pressureOutOfBounds))
{
    WarningInFunction
        << "Pressure field exceeded physical range while pressureClamp=false"
        << nl
        << "    min(p) = " << observedMinP.value() << " Pa, max(p) = "
        << observedMaxP.value() << " Pa" << nl
        << "    Applying emergency clamp to [" << minPressure.value() << ", "
        << maxPressure.value() << "] Pa" << endl;
        
    // Lines 300-334: Clamp pressure field to safe bounds
    scalarField& pInternal = p.primitiveFieldRef();
    forAll(pInternal, celli)
    {
        scalar& val = pInternal[celli];
        if (!std::isfinite(val))
        {
            val = 0.0;  // Reset non-finite values to zero
        }
        // Clamp to [minPressure, maxPressure]
        val = Foam::max(minPressure.value(), 
                        Foam::min(val, maxPressure.value()));
    }
}
```

### Critical Issue: Post-Correction Clamping

**Why this is problematic**:

1. **Occurs AFTER pressure solve** → Pressure field is discontinuous
2. **Creates artificial pressure gradients** → Large ∇p drives velocity spikes
3. **Magnitude**: Pressure clamped from -8.8e+13 Pa → [-1e7, 1e9] Pa is a **jump of 10^20+ Pa over cell widths**
4. **Affects momentum indirectly** → Velocity computed from pressure gradient on next iteration

**Pressure gradient impact** (UEqn.H, lines 334-430):
```cpp
// Lines 334-346: Pressure gradient computation
tmp<volVectorField> tGradPrgh = fvc::grad(p_rgh);
volVectorField gradPrgh(...);

// Lines 349-430: Gradient clamping (IF configured)
const Foam::scalar gradLimit = solverCoeffs.lookupOrDefault<Foam::scalar>
(
    "maxPressureGradient",
    Foam::GREAT  // DEFAULT: NO LIMIT!
);

// If maxPressureGradient is not set, unlimited gradient can develop
if (gradLimit < Foam::GREAT)
{
    // Clamp gradient magnitude
    if (magG > gradLimit && magG > SMALL)
    {
        g *= gradLimit/magG;
    }
}
```

**With default settings**: `maxPressureGradient` is NOT configured → pressure gradients are **unclamped**

---

## 3. TWO-TEMPERATURE MODEL COUPLING (TEqn.H, Lines 1097-1108)

### Temperature Equation
```cpp
// File: /home/user/compInterFoam/TEqn.H
// Lines 1097-1108: Gas temperature equation with interface coupling

fvScalarMatrix TEqn
(
    fvm::ddt(CvolLimited, T)
  + fvm::div(CvolPhi, T)
  - fvm::Sp(contErr, T)
  - fvm::laplacian(transportModel.kappaEff(), T)
  ==
    fvOptions(CvolLimited, T)
);

// Add gas-metal interface exchange:
TEqn += fvm::Sp(relaxedInterfaceExchange, T);  // Line 1107
TEqn.source() += relaxedInterfaceExchange*Tl;   // Line 1108
```

### Coupling Coefficient (TEqn.H, Lines 91-174)
```cpp
// Lines 91-94: Gas-metal exchange coefficient field
tmp<volScalarField> tInterfaceExchange = ttm.gasMetalExchangeCoeffField();
volScalarField& interfaceExchange = tInterfaceExchange.ref();
tmp<volScalarField> tCappedInterfaceExchange(interfaceExchange.clone());
volScalarField& cappedInterfaceExchange = tCappedInterfaceExchange.ref();

// Lines 121-147: Capping the exchange coefficient
const dimensionedScalar minCoeff =
    twoTemperatureDict.lookupOrDefault<dimensionedScalar>
    (
        "minGasMetalExchangeCoeff",
        minCoeffDefault
    );

const dimensionedScalar maxCoeff =
    twoTemperatureDict.lookupOrDefault<dimensionedScalar>
    (
        "maxGasMetalExchangeCoeff",
        maxCoeffDefault  // DEFAULT: GREAT (no upper limit!)
    );
```

### Critical Coupling Issue

From **twoTemperatureModel.H lines 55-60**:
```
gasMetalExchangeCoeff
    Interfacial gas-metal coupling coefficient [W/m³/K].
    Values in the 1e12-1e15 range provide moderate
    coupling, while coefficients above 1e16 W/m³/K
    enforce near-equilibrium interface temperatures and
    can lead to stiff coupling.
```

**In your simulation with Te = 11800 K**:
- If coupling coefficient > 1e16 W/m³/K, coupling becomes **stiff**
- Stiff coupling forces unrealistic temperature equilibration
- This creates sharp density gradients ∇ρ = ρ(T) · ∇T
- Density gradients drive large pressure gradients: ∇p ~ ∇ρ · sound_speed²

**Result**: Even with **zero recoil pressure**, the two-temperature model can generate pressure oscillations through:
1. Large T gradients at metal/gas interface
2. Unconstrained density variations
3. No regularization of ∇p in momentum equation

---

## 4. MOMENTUM EQUATION SOURCE TERMS (UEqn.H, Lines 1-528)

### Source Term Assembly (Lines 320-528)

**Only gravity, pressure, and recoil are included as momentum sources** (see lines 330-333):
```cpp
// Lines 330-333: Which forces are included
// The capillary force already appears in phig within pEqn.H via
// mixture.surfaceTensionForce(). Including it here would double-count.
// Gravity is already accounted for through -ghf*fvc::snGrad(rho) in phig.
// So we only add: -grad(p_rgh) + recoilForce
```

**Actual momentum source** (lines 432-528):
```cpp
// Lines 432-434: Pressure gradient force (unclamped)
tmp<volVectorField> tPressureForce(-gradPrgh);
volVectorField& pressureForce = tPressureForce.ref();

tmp<volVectorField> tMomentumSource(tPressureForce);

// Lines 439-484: Add recoil force if available
if (recoilForceReady)
{
    recoilForceFieldPtr = &tRecoilForce();
}

// Lines 513-528: Solve momentum with complete source
solve(UEqn == momentumSource);  // momentumSource = ∇p_rgh + recoilForce
```

### The Problem: Unconstrained Pressure Gradient

**Pressure gradient IS clamped** (lines 349-430), BUT:

1. **Default maxPressureGradient = GREAT** → No actual clamping
2. **Gradient computed from post-correction pressure field** → Discontinuities visible as large gradients
3. **With Δp ~ 1e13 Pa over Δx ~ 1e-6 m** → Gradient ~ 1e19 Pa/m (unlimited!)
4. **Acceleration scale**: a = ∇p/ρ ~ 1e13 m/s² (for ρ ~ 1 kg/m³)
5. **Over 1 timestep (Δt ~ 1e-15 s)**: ΔU ~ a·Δt ~ 1e-2 m/s/step, compounds to divergence

---

## 5. INTERFACE DISCONTINUITIES AT METAL/GAS BOUNDARY

### Material Property Discontinuities (TEqn.H, Lines 233-307)

The coupling between metal and gas exploits a **discontinuous masking strategy**:

```cpp
// Lines 238-270: Active metal mask with blending width
tmp<volScalarField> tActiveMask;
if (blendWidth <= SMALL)
{
    tActiveMask = Foam::pos(metalFraction - cutoffDim);  // SHARP DISCONTINUITY
}
else
{
    // Smooth blending across metalAmbientBlendWidth
    const scalar lowerBoundValue =
        Foam::max(metalFractionFloor - blendWidth, scalar(0));
    // ... linear blend ...
}
```

### Interface Identification (Lines 342-369)
```cpp
// Lines 342-369: Identify interface cells
tmp<volScalarField> tInterfaceMask(activeMask.clone());
volScalarField& interfaceMask = tInterfaceMask.ref();
scalarField& interfaceMaskInternal = interfaceMask.primitiveFieldRef();
interfaceMaskInternal = scalar(0);

forAll(interfaceMaskInternal, cellI)
{
    const scalar mf = metalFraction[cellI];
    const scalar gradMag = mag(gradMetal[cellI]);

    if (mf >= boundedMetalLower && mf <= boundedMetalUpper 
     && gradMag >= interfaceGradientCutoff)
    {
        interfaceMask[cellI] = scalar(1);  // This cell is at interface
    }
}

// Lines 366-369: Localize gas-metal exchange to interface
cappedInterfaceExchange = interfaceMask*cappedInterfaceExchange;
cappedInterfaceExchange.correctBoundaryConditions();
```

### Issue: Localized Heating Creates Pressure Spikes

1. **Interface cells get high coupling** (cappedInterfaceExchange > 0)
2. **Gas temperature rapidly rises** to match lattice temperature Tl
3. **Density drops**: ρ_gas ∝ 1/T
4. **Pressure must adjust**: But momentum equation is not synchronized
5. **Result**: Localized pressure pockets → ∇p spikes at interface

---

## 6. COUPLING SEQUENCE IN TIME LOOP (compInterFoam.C, Lines 919-977)

### Execution Order (Lines 919-977)
```cpp
while (pimple.loop())
{
    mixture.correct();  // Update phase fractions
    
    #include "compressibleAlphaEqnSubCycle.H"
    transportModel.correctPhasePhi();
    mixture.correct();
    
    // Lines 927-944: THERMAL COUPLING (independent of momentum)
    for (label thermalIter = 0; thermalIter < nThermalCouplingIter; ++thermalIter)
    {
        #include "TEqn.H"  // ← Solve gas temperature (coupled to Tl but independent of U,p)
        
        ttm.solve(...)     // ← Solve electron/lattice temperatures (independent of U,p)
    }
    
    // Lines 946-947: Update phase properties
    mixture.correct();  // ← Updates density based on NEW temperatures
    
    // Lines 949-960: CALCULATE RECOIL PRESSURE
    if (useAdvancedCapturing && pInterfaceCapturing.valid())
    {
        pInterfaceCapturing->calculateRecoilPressure();  // Uses NEW temperatures
    }
    
    // Lines 966-967: MOMENTUM SOLVE
    #include "UEqn.H"  // ← Momentum predictor (recoil is now available)
    
    // Lines 970-977: PRESSURE CORRECTION
    while (pimple.correct())
    {
        pInterfaceCapturing->calculateRecoilPressure();  // Recalc recoil
        #include "pEqn.H"  // ← Pressure correction (uses updated rhoPhi)
    }
}
```

### The Coupling Problem

**Temperature and momentum are decoupled** in the sequence:

1. **TEqn.H solves for T_gas** using OLD velocity field (U^n)
2. **ttm.solve()** solves for Te, Tl independently
3. **New density** ρ^(n+1) computed from new temperatures
4. **Momentum uses old density** ρ^n with new recoil pressure
5. **Pressure correction tries to reconcile** but happens AFTER momentum equation
6. **Discontinuity**: ∇(ρ^(n+1)) + U^(n+1/2) forces are out of sync

**With temperature-driven density changes**:
- ρ^(n+1) can differ significantly from ρ^n
- This invalidates the SIMPLE/PIMPLE pressure-velocity coupling assumption
- Results in uncontrolled pressure oscillations

---

## 7. SPECIFIC NUMERICAL ISSUES WITH ZERO RECOIL

### Why Simulation Still Fails Despite recoilPressure = 0

The recoil pressure is defined as:
```
recoilPressure = (interface_strength) × (1 + α(T-Tvap)) × normal
```

When recoil = 0:
- The **temperature coupling is still active**
- Gas-metal interface exchange: ` = cappedInterfaceExchange * (Tl - T)`
- This drives T_gas toward Tl
- T_gas increase → ρ_gas decrease → local pressure drop
- Pressure drop → velocity acceleration
- Velocity acceleration → kinetic energy increase

### Mechanism of Instability

1. **Electron temperature Te = 11800 K** (very high)
2. **Lattice temperature Tl ≈ 11800 K** (equilibrated with electrons)
3. **Gas temperature T << 11800 K** (ambient initially)
4. **Interface exchange coefficient** `h = maxCoeff` (possibly very large)
5. **Heat flux into gas**: q = h·(Tl - T) (LARGE!)
6. **Gas temperature rises** → causes density drop
7. **Density drop creates pressure pocket** (continuity + momentum must balance)
8. **Pressure clamping** after-the-fact → discontinuous ∇p
9. **Velocity spikes from ∇p** → kinetic energy overflow

---

## SUMMARY: ROOT CAUSES

| Issue | Location | Mechanism | Impact |
|-------|----------|-----------|--------|
| **Unconstrained pressure gradient** | UEqn.H:349-430 | maxPressureGradient defaults to GREAT | Large a = ∇p/ρ |
| **Stiff thermal coupling** | TEqn.H:1107-1108 | Exchange coeff unbounded (>1e16) | Uncontrolled density variations |
| **Post-correction pressure clamp** | pEqn.H:290-334 | Discontinuous pressure adjustment | Artificial ∇p spikes |
| **Decoupled T-U-p solve** | compInterFoam.C:919-977 | Temperature solves independently | Momentum/continuity mismatch |
| **Missing density sources** | pEqn.H:72-79 | Pressure equation ignores ∂ρ/∂t | Pressure-density uncoupling |
| **Velocity overshoot** | UEqn.H:56-82 | No velocity limiting before solve | Kinetic energy check fails |

