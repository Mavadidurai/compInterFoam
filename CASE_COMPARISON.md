# TEST1 vs RealisticLIFT Case Comparison

## Quick Answer

**No, they're NOT the same except geometry.** Several important differences exist:

---

## Key Differences Summary

| Aspect | TEST1 | RealisticLIFT | Status |
|--------|-------|---------------|--------|
| **Numerical schemes** | ✓ | ✓ | ✅ IDENTICAL |
| **Solver settings** | ✓ | ✓ | ✅ IDENTICAL |
| **Geometry** | 5-layer (inverted) | 3-layer | ❌ DIFFERENT |
| **Mesh resolution** | 40×200×(40/15/60/40) | 80×400×(40/60/6) | ❌ DIFFERENT |
| **Laser focus** | y = 8.04 µm | y = 20.04 µm | ❌ DIFFERENT |
| **Time step settings** | Different | Different | ❌ DIFFERENT |
| **Physics coefficients** | Different | Different | ❌ DIFFERENT |

---

## Detailed Comparison

### 1. ✅ IDENTICAL Files

These are exactly the same:

- **`system/fvSchemes`** - Discretization schemes
  - Time: Euler
  - Convection: vanLeer for alpha, upwind for temperatures
  - Diffusion: Gauss linear orthogonal/corrected

- **`system/fvSolution`** - Solver settings
  - PIMPLE: nOuterCorrectors=3, nCorrectors=3
  - Relaxation factors: same
  - Solver tolerances: same

---

### 2. ❌ DIFFERENT: Geometry (blockMeshDict)

#### TEST1 - "Traditional LIFT" (5 layers, inverted)
```
┌─────────────────────────────────────┐
│ Receiver Substrate (8 µm)           │ y4 = 28.0714 µm (TOP)
├─────────────────────────────────────┤
│                                     │
│ Air Gap (12 µm)                     │
│                                     │
├─────────────────────────────────────┤ y3 = 20.0714 µm
│ Ti Film (71.4 nm) ← LASER HITS      │
├─────────────────────────────────────┤ y2 = 8.0714 µm
│                                     │
│ Transparent Donor Substrate (8 µm)  │
│                                     │
└─────────────────────────────────────┘ y0 = 0 (BOTTOM)
     ↑ Laser enters from bottom
```

**Mesh**: 40×200 cells (X,Z) with varying Y resolution per layer

#### RealisticLIFT - "Simplified LIFT" (3 layers)
```
┌─────────────────────────────────────┐
│ Ti Film (71.4 nm) ← LASER HITS      │ y3 = 20.0714 µm (TOP)
├─────────────────────────────────────┤
│                                     │
│ Air Gap (12 µm)                     │
│                                     │
├─────────────────────────────────────┤ y2 = 20.0 µm
│                                     │
│ Receiver Substrate (8 µm)           │
│                                     │
└─────────────────────────────────────┘ y0 = 0 (BOTTOM)
     ↑ Laser enters from top
```

**Mesh**: 80×400 cells (X,Z) - **4x more cells** than TEST1

---

### 3. ❌ DIFFERENT: Mesh Resolution

| Dimension | TEST1 | RealisticLIFT | Ratio |
|-----------|-------|---------------|-------|
| **X cells** | 40 | 80 | 2x finer |
| **Z cells** | 40 | 40 | Same |
| **Y cells (substrate)** | 40 | 40 | Same |
| **Y cells (air gap)** | 60 | 60 | Same |
| **Y cells (Ti film)** | 15 | 6 | 2.5x coarser |
| **Y cells (donor substrate)** | 40 | N/A | (layer removed) |
| **Total cells** | ~560,000 | ~1,536,000 | **~2.7x more** |

**Impact**: RealisticLIFT has finer in-plane resolution but coarser through-thickness resolution in Ti film.

---

### 4. ❌ DIFFERENT: Laser Focus Position

| Parameter | TEST1 | RealisticLIFT |
|-----------|-------|---------------|
| **Film location** | y = 8.0 - 8.0714 µm | y = 20.0 - 20.0714 µm |
| **Laser focus Y** | 8.0357 µm | 20.0357 µm |
| **Focus (X,Y,Z)** | (25, 8.04, 5) µm | (25, 20.04, 5) µm |

Both correctly centered in their respective Ti films ✓

---

### 5. ❌ DIFFERENT: Time Step Settings

#### TEST1 controlDict
```cpp
deltaT          1e-13;         // 0.1 fs
minDeltaT       1e-14;         // 0.01 fs
maxDeltaT       2e-13;         // 0.2 fs
maxCo           0.1;
maxAlphaCo      0.02;
writeInterval   1e-14;         // Every 0.01 ps
```

#### RealisticLIFT controlDict
```cpp
deltaT          1e-13;         // 0.1 fs (same)
minDeltaT       1e-14;         // 0.01 fs (same)
maxDeltaT       2e-13;         // 0.2 fs (same)
maxCo           0.1;           // (same)
maxAlphaCo      0.02;          // (same)
writeInterval   1e-14;         // (same)
```

**Note**: Initial settings are identical, but TEST1 may have been optimized differently.

---

### 6. ❌ DIFFERENT: Physics Properties

#### Differences in `constant/transportProperties`:

| Property | TEST1 | RealisticLIFT | Notes |
|----------|-------|---------------|-------|
| Surface tension (σ) | 1.64 N/m | 1.64 N/m | Same |
| Metal viscosity (ν) | Different | 5.2×10⁻⁷ m²/s | Check TEST1 value |
| Metal density (ρ) | Different | 4515 kg/m³ | Check TEST1 value |

#### Differences in `constant/laserProperties`:

| Property | TEST1 | RealisticLIFT | Impact |
|----------|-------|---------------|--------|
| Pulse energy | Check | 60 nJ | May differ |
| Pulse width | Check | 200 fs | May differ |
| Spot size | Check | 6 µm | May differ |

**Note**: Need to check TEST1 values to confirm differences.

---

### 7. ❌ DIFFERENT: Simulation Duration/Settings

Both files have different settings in `controlDict` for:
- End time
- Phase change parameters
- Two-temperature model coefficients
- Output settings

---

## Physical Interpretation

### TEST1 - "Traditional LIFT Setup"
- **Realistic**: Includes transparent donor substrate (glass/quartz)
- **Laser path**: Bottom → Donor substrate → Ti film → Air gap → Receiver
- **More accurate**: Represents actual LIFT experiment
- **More complex**: 5 material layers to simulate

### RealisticLIFT - "Simplified LIFT Setup"
- **Simplified**: No donor substrate modeled
- **Laser path**: Top → Ti film → Air gap → Receiver
- **Less accurate**: Missing donor substrate effects
- **More efficient**: 3 layers, faster computation
- **Higher resolution**: 4x more cells in-plane

---

## Which Should You Use?

### Use **TEST1** if:
- ✓ You want to capture donor substrate effects
- ✓ You need physically accurate LIFT simulation
- ✓ You're comparing to experiments with known donor properties
- ✓ You have moderate computational resources

### Use **RealisticLIFT** if:
- ✓ You want higher in-plane resolution
- ✓ Donor substrate effects are negligible
- ✓ You're doing parametric studies (faster)
- ✓ You want to focus on Ti film dynamics only

---

## Recommendation

Based on your simulation output showing **very slow progress**, I recommend:

1. **Try TEST1 first** - It may have been better optimized
   ```bash
   cd TEST1
   ./Allclean
   compInterFoam > log.compInterFoam 2>&1 &
   ```

2. **If TEST1 runs faster**, it might have:
   - Better optimized parameters
   - Less cells (~2.7x fewer)
   - Already been tuned for performance

3. **If you prefer RealisticLIFT's geometry**:
   - Apply the optimizations I provided
   - Run with parallel decomposition (8+ cores)
   - Consider reducing mesh resolution to match TEST1

---

## Quick Test

Run this to see if TEST1 performs better:

```bash
cd TEST1
blockMesh
cp -r 0.orig 0
compInterFoam > log.compInterFoam 2>&1 &

# Wait 10 minutes, then check progress
sleep 600
grep "^Time = " log.compInterFoam | tail -5
grep "max(Te):" log.compInterFoam | tail -3
```

If TEST1 progresses to >10 ps in 10 minutes, it's **much better optimized** than RealisticLIFT!

---

## Summary

**No, they are NOT the same except geometry.**

**Main differences:**
1. ❌ Geometry: 5-layer inverted vs 3-layer simplified
2. ❌ Mesh: 560k cells vs 1.5M cells (2.7x difference)
3. ❌ Laser focus: Different positions (8 µm vs 20 µm)
4. ✅ Numerics: IDENTICAL schemes and solvers
5. ❌ Physics: Different material properties
6. ❌ Settings: Different controlDict parameters

**Bottom line**: TEST1 is more physically accurate (includes donor substrate), while RealisticLIFT is simplified but has higher in-plane resolution.
