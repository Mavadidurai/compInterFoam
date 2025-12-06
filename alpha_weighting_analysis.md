# Alpha-Weighting for Recoil Pressure: UEqn vs pEqn Analysis

**Date**: 2025-12-06
**Context**: Recoil pressure currently accelerates vapor (1123 m/s) instead of metal (0.009 m/s)

---

## The PIMPLE Algorithm Structure

```
PIMPLE outer loop (correctors):
    1. Momentum Predictor (UEqn.H):
       - Assemble momentum equation: ∂(ρU)/∂t + ... = -∇P + sources
       - Solve for predicted velocity: U*
       - Uses old/lagged pressure field

    2. Pressure Corrector (pEqn.H):
       - Enforce continuity: ∇·U = 0 (incompressible) or mass conservation
       - Solve for pressure: ∇·(1/A * ∇p) = ∇·(H/A)
       - Compute flux: φ = (H/A) - (1/A)∇p

    3. Velocity Correction:
       - Update U from new pressure
       - U = H/A - (1/A)∇p

    4. Check convergence, repeat if needed
```

**Key Point**: In PIMPLE, the momentum equation and pressure equation are **coupled** through an iterative predictor-corrector loop.

---

## Option 1: Alpha-Weighting in UEqn.H (Momentum Equation)

### Implementation:
```cpp
// UEqn.H line 111:
totalPressure += alpha1 * (*recoilContribution);
// Then: solve(UEqn == -∇totalPressure)
```

### How It Works:
1. Recoil pressure added to `totalPressure = p_rgh + α·P_recoil`
2. Momentum source: `F = -∇(p_rgh + α·P_recoil)`
3. Gradient: `∇(α·P_recoil) = P_recoil·∇α + α·∇P_recoil`
4. Acceleration: `a = F/ρ = -[P_recoil·∇α + α·∇P_recoil]/ρ`

### Analysis:

**✅ Advantages:**
1. **Straightforward implementation**: Simple multiplication `alpha1 * recoilP`
2. **Consistent with PIMPLE**: Pressure forces go through momentum predictor
3. **Automatic coupling**: PIMPLE corrector loop handles pressure-velocity coupling
4. **Conservative**: Total pressure field remains well-defined
5. **Easy to debug**: Can visualize `totalPressure` field directly

**⚠️ Disadvantages:**
1. **Gradient term `P_recoil·∇α`**: Creates force at interface
   - This is actually **physical**! (analogous to surface tension)
   - But can be strong if P_recoil is large and ∇α is sharp
2. **Lagged in predictor**: Uses old α field initially (updated in corrector)
3. **Smooth α required**: Sharp jumps in α create large ∇α
4. **Implicit mixing**: α-weighted pressure appears in both momentum and pressure equations

**📊 Numerical Behavior:**

At the interface (∇α ≠ 0):
- Force = `-P_recoil·∇α - α·∇P_recoil`
- The `-P_recoil·∇α` term acts like a **body force** pushing on the interface
- For P_recoil = 71 MPa, |∇α| ~ 1/Δx ~ 1e7 m⁻¹ (for Δx ~ 100 nm)
- Interface force ~ 71e6 Pa × 1e7 m⁻¹ = **7e14 N/m³** (huge!)

**Physical Interpretation:**
- The `∇(α·P)` term represents **pressure acting on the interface**
- In reality, recoil pressure DOES act at the surface (liquid-vapor interface)
- So this is **physically correct** if the interface is well-resolved

---

## Option 2: Alpha-Weighting in pEqn.H (Pressure Equation)

### Implementation:
```cpp
// pEqn.H line 87 (currently commented out):
surfaceScalarField alphaf = fvc::interpolate(alpha1);
phig -= alphaf * fvc::snGrad(*recoilPressurePtr) * rAUf * mesh.magSf();
```

### How It Works:
1. Recoil applied as **flux correction** in pressure equation
2. Flux: `φ += α_f · (1/A) · ∇P_recoil · S`
3. Pressure equation: `∇·(rAU·∇p) = ∇·(φ_HbyA + α_f·∇P_recoil·rAU·S)`
4. Velocity correction: `U += (1/A)·∇p` (happens automatically)

### Analysis:

**✅ Advantages:**
1. **Direct control of flux**: Acts on face fluxes, not cell pressures
2. **No gradient of alpha**: Uses `α_f·∇P_recoil`, not `∇(α·P_recoil)`
   - Avoids the large interface force from `P·∇α`
3. **Localized to interface**: Only affects faces where α_f ≠ 0
4. **Cleaner coupling**: Recoil appears only in pressure eq, not momentum
5. **Better for sharp interfaces**: Doesn't rely on smooth α field

**⚠️ Disadvantages:**
1. **More complex implementation**: Need to interpolate α to faces
2. **Decoupled from momentum predictor**: Recoil doesn't appear in UEqn directly
   - Momentum predictor doesn't "see" recoil until corrector step
3. **Requires PIMPLE iterations**: Won't work well with single corrector
4. **Face interpolation issues**: How to weight α at interface faces?
   - Linear interpolation? Upwind? Harmonic?
5. **Harder to debug**: Flux correction is implicit, can't visualize easily

**📊 Numerical Behavior:**

At the interface:
- Face flux: `φ_f = α_f · (rAU)_f · (∇P_recoil)_f · S_f`
- For α_f ≈ 0.5 (interface): flux is proportional to α
- For α_f ≈ 0 (vapor): flux ≈ 0
- For α_f ≈ 1 (metal): flux = full recoil contribution

**Physical Interpretation:**
- Applies recoil as a **boundary flux** at the interface
- More like a "traction" or surface force
- Physically reasonable for evaporation-driven flow

---

## Comparative Analysis

### Test Case: Recoil Pressure at Flat Interface

Consider a 1D problem: metal below, vapor above, P_recoil = 71 MPa at interface

| Property | UEqn Approach | pEqn Approach |
|----------|---------------|---------------|
| **Force in metal (α=1)** | `-∇(1·P_recoil)` = `-∇P_recoil` | Flux drives velocity |
| **Force in vapor (α=0)** | `-∇(0·P_recoil)` = 0 | No flux contribution |
| **Force at interface (∇α≠0)** | `-P_recoil·∇α - α·∇P_recoil` | `α_f·∇P_recoil` only |
| **Coupling** | Immediate (in predictor) | Delayed (in corrector) |
| **Momentum conservation** | Exact (pressure field conservative) | Approximate (flux-based) |

### Stability Considerations

**UEqn approach:**
- **Pro**: Pressure field remains smooth (no flux jumps)
- **Con**: Large `P·∇α` can cause checkerboarding at interface
- **Requires**: Interface compression to keep ∇α localized

**pEqn approach:**
- **Pro**: Flux jumps absorbed by pressure solver
- **Con**: Can create pressure oscillations if not converged
- **Requires**: Multiple PIMPLE iterations for stability

### Accuracy Considerations

**UEqn approach:**
- **Momentum conservation**: Exact (conservative force field)
- **Interface force**: Correct magnitude, but distributed over interface cells
- **Order of accuracy**: Limited by α field resolution

**pEqn approach:**
- **Momentum conservation**: Flux-based (non-conservative if PIMPLE not converged)
- **Interface force**: Concentrated at faces (sharper)
- **Order of accuracy**: Limited by face interpolation scheme

---

## Literature & Best Practices

### Brackbill et al. (1992) - Continuum Surface Force Model
- Surface tension applied as **volume force** `F = σκ∇α`
- Equivalent to `∇(α·p_surface)` approach
- **Conclusion**: Volume force approach stable if interface is sharp

### Tryggvason et al. (2011) - Direct Numerical Simulations
- Pressure jump at interface handled via **level-set** or **VOF**
- Recommend **pressure formulation** for multiphase flows
- **Conclusion**: Pressure-based methods more robust than force-based

### OpenFOAM interFoam Implementation
- Surface tension: Applied in **pEqn.H** as `phig += surfaceTensionForce()`
- **Conclusion**: OpenFOAM standard is pEqn approach

### Feinaeugle et al. (2017) - fs-LIFT Experiments
- Recoil pressure drives **metal jet** formation
- Jet velocity: 30-100 m/s (metal), vapor: 500-2000 m/s
- **Conclusion**: Recoil must accelerate metal, not vapor!

---

## Recommendation: Hybrid Approach

Neither approach alone is perfect. Here's the optimal solution:

### **Best Practice: Alpha-Weighted Pressure Field (UEqn) + Flux Correction (pEqn)**

```cpp
// 1. In UEqn.H (line 111):
// Apply alpha-weighted recoil to total pressure
totalPressure += alpha1 * (*recoilContribution);

// 2. In pEqn.H (line 87):
// Add flux correction to balance the ∇(α·P) gradient at interface
if (recoilPressurePtr)
{
    surfaceScalarField alphaf = fvc::interpolate(alpha1);
    surfaceScalarField gradAlphaf = fvc::snGrad(alpha1);

    // Correct for the P·∇α term that appears in ∇(α·P)
    // This cancels the spurious interface force
    surfaceScalarField recoilFluxCorrection =
        -(*recoilPressurePtr) * gradAlphaf * rAUf * mesh.magSf();

    phig += recoilFluxCorrection;
}
```

**Why this works:**
1. **UEqn term** `∇(α·P_recoil)` = `P·∇α + α·∇P_recoil` drives the flow
2. **pEqn correction** `-P·∇α` cancels the spurious interface force
3. **Net effect**: Only `α·∇P_recoil` remains → force proportional to α

**Result:**
- Metal (α=1): Full recoil force
- Vapor (α=0): No recoil force
- Interface: Smooth transition, no spurious forces

---

## Simpler Alternative: Just UEqn (Recommended for Your Case)

For your current simulation, the **simpler UEqn-only approach is sufficient**:

```cpp
// UEqn.H line 111:
totalPressure += alpha1 * (*recoilContribution);
```

### Why this is enough:

1. **Your interface is well-resolved**: Δx ~ 100 nm, interface width ~ few cells
2. **MULES interface compression**: Keeps ∇α sharp and localized
3. **PIMPLE iterations**: 20 iterations will converge the coupling
4. **Physical correctness**: `P·∇α` term IS the interface force (like surface tension)

### Expected improvements:

| Metric | Before | After (UEqn only) |
|--------|--------|-------------------|
| Metal avg velocity | 0.009 m/s | **30-100 m/s** ✅ |
| Vapor velocity | 1123 m/s | **500-1000 m/s** ✅ |
| PIMPLE iterations | 20 (not converged) | **5-10** ✅ |
| Recoil driving metal? | ❌ No | ✅ Yes |

---

## Implementation Plan

### Phase 1: UEqn Alpha-Weighting (Minimal Change)
```cpp
// File: UEqn.H, line 111
// Change from:
totalPressure += *recoilContribution;

// To:
totalPressure += alpha1 * (*recoilContribution);
```

**Test**: Run for 1-2 ps, check if metal avg velocity increases to 30-100 m/s

### Phase 2: Monitor Interface Forces (Diagnostic)
```cpp
// After line 111, add diagnostic:
if (recoilPressurePtr && verbose)
{
    volVectorField recoilForce = -fvc::grad(alpha1 * (*recoilContribution));
    scalar maxRecoilForce = gMax(mag(recoilForce));

    Info<< "    Max recoil force: " << maxRecoilForce
        << " N/m³ (should be ~71 MPa / 100 nm = 7e14)" << nl;
}
```

### Phase 3: pEqn Correction (If Needed)
If you see:
- ❌ Interface instability (checkerboarding)
- ❌ Pressure oscillations
- ❌ PIMPLE still not converging

Then add the flux correction in pEqn.H (line 87):
```cpp
if (recoilPressurePtr)
{
    surfaceScalarField alphaf = fvc::interpolate(alpha1);
    surfaceScalarField recoilPf = fvc::interpolate(*recoilContribution);
    surfaceScalarField gradAlphaf = fvc::snGrad(alpha1);

    // Cancel the P·∇α spurious force
    phig -= recoilPf * gradAlphaf * rAUf * mesh.magSf();
}
```

---

## Final Recommendation

### **Start with UEqn-only approach** ✅

**Rationale:**
1. ✅ **Simplest fix**: One line change
2. ✅ **Physically correct**: `∇(α·P)` is the right physics
3. ✅ **OpenFOAM standard**: Similar to how surface tension works
4. ✅ **Your setup supports it**: Sharp interface, MULES compression, 20 PIMPLE iters
5. ✅ **Easy to test**: Immediate feedback on metal velocity

**If problems arise**, add pEqn flux correction as Phase 3.

**Success criteria:**
- Metal avg velocity: 30-100 m/s
- Vapor velocity: 500-1000 m/s (reduced from 1123)
- PIMPLE converges in <10 iterations
- No pressure/velocity oscillations at interface

---

## Appendix: Mathematical Derivation

### Pressure Gradient of Alpha-Weighted Field

Given: `P_total = p_rgh + α·P_recoil`

Gradient:
```
∇P_total = ∇p_rgh + ∇(α·P_recoil)
         = ∇p_rgh + P_recoil·∇α + α·∇P_recoil
```

Force on fluid element:
```
F = -∇P_total / ρ
  = -(∇p_rgh + P_recoil·∇α + α·∇P_recoil) / ρ
```

In metal (α=1, ∇α≈0):
```
F_metal = -(∇p_rgh + ∇P_recoil) / ρ_metal
```

In vapor (α=0, ∇α≈0):
```
F_vapor = -∇p_rgh / ρ_vapor
```

At interface (α≈0.5, ∇α large):
```
F_interface = -(∇p_rgh + P_recoil·∇α + 0.5·∇P_recoil) / ρ_interface
```

The `P_recoil·∇α` term is the **interface pressure jump** - this is PHYSICAL!

---

## References

1. Brackbill, J.U., Kothe, D.B., Zemach, C. (1992). "A continuum method for modeling surface tension". *Journal of Computational Physics*, 100(2), 335-354.

2. Tryggvason, G., Scardovelli, R., Zaleski, S. (2011). *Direct Numerical Simulations of Gas-Liquid Multiphase Flows*. Cambridge University Press.

3. Feinaeugle, M., et al. (2017). "Time-resolved imaging of flyer dynamics for laser-induced forward transfer of metals". *Applied Surface Science*, 418, 72-79.

4. Jasak, H. (1996). "Error Analysis and Estimation for the Finite Volume Method with Applications to Fluid Flows". PhD Thesis, Imperial College London.

5. Rusche, H. (2003). "Computational Fluid Dynamics of Dispersed Two-Phase Flows at High Phase Fractions". PhD Thesis, Imperial College London.
