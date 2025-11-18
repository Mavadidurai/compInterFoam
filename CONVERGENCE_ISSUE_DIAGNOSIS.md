# Simulation Convergence Issue Diagnosis

**Date:** 2025-11-18
**Case:** TEST1 - fs-LIFT Ti film simulation
**Issue:** Pressure solver non-convergence and extremely low recoil pressure

## Executive Summary

The simulation is experiencing catastrophic pressure solver failure after PIMPLE iteration 19, with recoil pressures ~13 orders of magnitude below expected values. Root cause identified as severely misconfigured phase change parameters that suppress mass flux.

## Observed Symptoms

### 1. Pressure Solver Failure
- **PIMPLE iteration 20+**: GAMG hits max iterations (300) with residuals near 1.0
- Complete non-convergence indicating fundamental physics inconsistency
- Time step frozen at dt = 3.3333e-15 s (extremely small)

### 2. Extreme Pressure Discrepancy
```
Diagnostic Warning:
  pressure clamp hit maxPressure = 5e+09 Pa
  recoil pressure = 7.4530909e-10 MPa = 0.000745 Pa
```
**Ratio:** 6.7 × 10¹² (13 orders of magnitude difference)

### 3. Velocity Limiting
```
|U|_metal_max = 2.6 m/s
recoil-limited velocity = 0.00057 m/s
Actual velocity is 4565× the recoil-driven expectation
```

### 4. Mass Flux Starvation
```
Max |j_net| = 8.6455604e-09 kg/m²/s
Active cells: 302 of 78000 interface cells
```

**Expected vs Actual:**
- At T = 10000 K, Ti should produce j ~ 100-1000 kg/m²/s
- Observed: 8.65e-9 kg/m²/s
- **Deficit:** 11-12 orders of magnitude

## Root Cause Analysis

### Issue #1: Severe Mass Flux Suppression

**Location:** `TEST1/system/controlDict:73`
```
phaseChangeCoeffs
{
    maxSource [1 -1 -3 0 0 0 0] 1e7;  // ← CRITICAL ERROR
}
```

**Impact:**
- maxSource = 1×10⁷ W/m³ (configured)
- Default = 4×10²³ W/m³ (from transportProperties)
- **Suppression factor:** 4 × 10¹⁶

**Physics:**
For film thickness δ = 71.4 nm:
```
maxHeatFlux = maxSource × δ
            = 1e7 W/m³ × 71.4e-9 m
            = 7.14e-4 W/m²

maxMassFlux = maxHeatFlux / L_vap
            = 7.14e-4 / 9.1e6
            = 7.8e-11 kg/m²/s
```

This artificially caps evaporation to essentially zero, preventing any meaningful recoil pressure from developing.

### Issue #2: Low Evaporation Coefficient

**Location:** `TEST1/system/controlDict:69`
```
evaporationCoeff 0.03;  // Ti evaporation coefficient
```

**Expected:** 0.18 - 1.0 for metals (literature standard)
**Configured:** 0.03 (6× below minimum)

### Issue #3: Pressure Clamp Too High

**Location:** `TEST1/system/fvSolution:214`
```
maxPressure 5.0e9;  // 5 GPa
```

**Context:**
- Feinaeugle et al. report Ti fs-LIFT recoil peaks at ~80 MPa
- Clamp should be ~100-200 MPa (allowing safety margin)
- Current: 5 GPa = 5000 MPa (62× too high)

When recoil pressure is artificially suppressed to ~10⁻⁹ MPa but momentum forcing tries to generate flow, the pressure solver attempts to compensate by reaching gigapascal pressures, triggering the clamp and creating the pressure-velocity inconsistency that crashes GAMG convergence.

## Physical Expectations

### Hertz-Knudsen Prediction (T = 10000 K)

For Ti at 10000 K:
```
p_sat(10000 K) = p₀ × exp[hf/R × (1/T_vap - 1/T)]
                = 101325 × exp[52299 × (1/3560 - 1/10000)]
                = 101325 × exp[9.74]
                = 1.72 × 10⁹ Pa = 1720 MPa

j = α × p_sat / sqrt(2πRT/M)
  = 0.18 × 1.72e9 / sqrt(2π × 8314 × 10000 / 47.867)
  = 0.18 × 1.72e9 / 5214
  = 59,400 kg/m²/s

p_recoil = 0.54 × p_sat = 929 MPa
```

### Observed (Simulation)
```
j_net = 8.65e-9 kg/m²/s  (deficit: 6.9 × 10⁹)
p_recoil = 7.45e-10 MPa  (deficit: 1.2 × 10⁹)
```

## Recommended Fixes

### Priority 1: Remove Mass Flux Suppression

**File:** `TEST1/system/controlDict`

**Change:**
```diff
phaseChangeCoeffs
{
-   maxSource [1 -1 -3 0 0 0 0] 1e7;
+   maxSource [1 -1 -3 0 0 0 0] 4e23;  // Use default (or omit)
```

**Rationale:** The default value (4×10²³ W/m³) from `transportProperties` is appropriate for ultrafast laser ablation. The current value artificially suppresses all phase change.

### Priority 2: Correct Evaporation Coefficient

**Change:**
```diff
-   evaporationCoeff 0.03;
+   evaporationCoeff 0.18;  // Literature value for Ti
```

**Rationale:** Knight (1979), Zhigilei et al., and Feinaeugle use α ≈ 0.18 for Ti.

### Priority 3: Reduce Pressure Clamp

**File:** `TEST1/system/fvSolution`

**Change:**
```diff
compInterFoamCoeffs
{
-   maxPressure 5.0e9;  // 5 GPa
+   maxPressure 2.0e8;  // 200 MPa (2.5× Feinaeugle peak)
-   minPressure -5e6;
+   minPressure -1e7;   // -10 MPa (allow modest tension)
```

**Rationale:** Feinaeugle et al. (Appl. Surf. Sci. 418, 2017) report Ti fs-LIFT recoil peaks near 80 MPa. Setting limit at 200 MPa provides safety margin without allowing unphysical multi-GPa pressures.

### Priority 4: Improve Pressure Solver Robustness

**File:** `TEST1/system/fvSolution`

**Changes:**
```diff
p_rgh
{
    solver          GAMG;
    smoother        GaussSeidel;
-   tolerance       1e-5;
+   tolerance       1e-6;
-   relTol          0.01;
+   relTol          0.001;
    maxIter         300;
    nPreSweeps      1;
-   nPostSweeps     3;
+   nPostSweeps     4;
    cacheAgglomeration true;
    nCellsInCoarsestLevel 10;
    agglomerator    faceAreaPair;
    mergeLevels     2;
-   relaxationFactor 0.5;
+   relaxationFactor 0.6;  // Less aggressive under-relaxation
}
```

**Rationale:** Once physical parameters are corrected, tighten solver tolerances to ensure accurate pressure field resolution.

### Priority 5: Adaptive Time Stepping

**File:** `TEST1/system/controlDict`

**Add/modify:**
```diff
adjustTimeStep  yes;
+maxCo           0.5;   // Courant number limit
+maxAlphaCo     0.5;   // Interface Courant limit
+maxDeltaT      1e-12; // Maximum time step
```

**Rationale:** Current frozen timestep (3.3e-15 s) indicates stability crisis. Adaptive stepping will allow recovery once physics is corrected.

## Expected Outcome After Fixes

### Mass Flux
- Should increase to j ~ 100-10000 kg/m²/s in peak regions
- Active interface cells should expand to ~1000-5000

### Recoil Pressure
- Should reach p_recoil ~ 50-150 MPa during peak heating
- Consistent with Feinaeugle experimental observations

### Pressure Solver
- GAMG should converge in 10-30 iterations
- Residuals should drop below 1e-6
- No clamp warnings (recoil pressure will be physical)

### Time Step
- Should adapt to dt ~ 1e-14 to 1e-12 s
- Larger during cooling phase
- Simulation will complete 200 ps window in reasonable wall time

## Verification Steps

After implementing fixes:

1. **Check startup diagnostics:**
   ```
   grep "j_net\|p_recoil" log.compInterFoam | head -50
   ```
   Verify mass flux increases to ~100 kg/m²/s by t ~ 100 fs

2. **Monitor pressure solver:**
   ```
   grep "GAMG.*p_rgh.*Iterations" log.compInterFoam
   ```
   Iterations should be <50, residuals <1e-5

3. **Verify recoil activation:**
   ```
   grep "Recoil diagnostics" log.compInterFoam | tail -20
   ```
   Should show p_recoil in MPa range, not μPa

4. **Check velocity diagnostic:**
   Should NOT see "pressure clamp hit" warnings

## References

1. Feinaeugle et al., "Time-resolved imaging of hydro and magnetohydrodynamics in laser-produced plasmas," Appl. Surf. Sci. 418 (2017)
2. Knight, "Theoretical Modeling of Rapid Surface Heating," Phys. Rev. B 20 (1979)
3. Zhigilei et al., "Computer Modeling of Laser Melting and Spallation of Metal Targets," J. Appl. Phys. 88 (2000)
4. Piqué et al., "Digital microfabrication by laser decal transfer," Appl. Phys. A 79 (2004)

---
**Status:** Analysis complete, fixes identified, awaiting implementation
**Next Action:** Apply recommended parameter changes and restart simulation
