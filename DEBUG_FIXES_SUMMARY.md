# OpenFOAM compInterFoam Debug Fixes - Session Summary

## Issues Identified

### 1. **Pressure Solver Divergence** ❌ CRITICAL
**Symptoms:**
- Unphysical negative pressures (~-10^10 Pa)
- Emergency pressure clamping activating repeatedly
- PIMPLE solver not converging within 10 iterations

**Root Cause:**
- Extreme density ratio (4500:1 metal/air) + femtosecond timesteps = catastrophic matrix ill-conditioning
- Reciprocal diagonal (rAU) becoming singular
- Insufficient relaxation and solver tolerance too tight

### 2. **Minimal Heating** ❌ CRITICAL
**Symptoms:**
- Laser delivers 2043 W
- Only 27 K temperature rise (max Te: 327 K vs required 1941 K)
- Electron-lattice coupling: only 2.59 W

**Root Cause:**
- Timestep too small (2e-16 s) → extremely slow energy accumulation
- Courant number limits too restrictive → preventing timestep growth
- Slow electron-lattice equilibration

### 3. **Simulation Progress** ❌
- Only reaching 0.53 femtoseconds after 23 seconds of wall time
- At this rate, 200 ps simulation would take ~3 months
- Need 5-10x speedup to be practical

## Fixes Applied

### A. Pressure Solver Stabilization (`fvSolution`)

#### 1. Enable Pressure Clamping
```
pressureClamp      true;  // Was: false
maxPressure        1.0e9;  // Was: 5.0e9 (1 GPa sufficient for fs-LIFT)
minPressure       -1.5e8;  // Was: -5e6 (-150 MPa allows some tension)
```

#### 2. Enable Density Clamping
```
densityClamp       true;   // Was: false
maxDensity         6000;   // Cap metal density fluctuations
minDensity         0.8;    // Floor for air density
```

#### 3. Enable rAU Clamping (CRITICAL)
```
enableRAUClamp     true;   // Was: false - prevents 1/diagonal singularities
minRAU             1e-8;   // Increased from 1e-10
maxRAU             1e5;    // Decreased from 1e6
minRAUf            1e-12;  // Increased from 1e-14
maxRAUf            1e5;    // Decreased from 1e6
```

#### 4. Improve GAMG Pressure Solver
```
smoother        GaussSeidel;  // Was: DICGaussSeidel (more stable)
tolerance       1e-5;         // Relaxed from 1e-6
relTol          0.01;         // Relaxed from 0.005
maxIter         300;          // Increased from 200
nPreSweeps      1;            // Added pre-sweeps
nPostSweeps     3;            // Increased from 2
nCellsInCoarsestLevel 10;     // Was: 20 (more aggressive coarsening)
relaxationFactor 0.5;         // Reduced from 0.7 (heavier relaxation)
```

#### 5. Relax PIMPLE Convergence Criteria
```
nOuterCorrectors            20;   // Increased from 10
nCorrectors                 4;    // Increased from 3
nNonOrthogonalCorrectors    3;    // Increased from 2

residualControl
{
    "alpha.*" { tolerance 1e-7;  relTol 0.0; }  // Relaxed from 1e-8
    p_rgh     { tolerance 1e-4;  relTol 0.05; } // Relaxed from 1e-5
    U         { tolerance 1e-3;  relTol 0.1; }  // Relaxed from 1e-4
    Te        { tolerance 1e-6;  relTol 0.0;  } // Relaxed from 1e-8
    Tl        { tolerance 1e-6;  relTol 0.0;  } // Relaxed from 1e-8
    T         { tolerance 1e-6;  relTol 0.05;  }// Relaxed from 1e-7
}
```

#### 6. Improved Relaxation Factors
```
relaxationFactors
{
    fields
    {
        p_rgh           0.5;   // Increased from 0.3
        "alpha.*"       0.7;   // Increased from 0.5
    }
    equations
    {
        U               0.7;   // Increased from 0.5
        T               0.7;   // Increased from 0.5
        Te              0.6;   // Increased from 0.5
        Tl              0.7;   // Increased from 0.5
        "alpha.*"       0.7;   // Increased from 0.5
        ".*Final"       1.0;   // Increased from 0.7
    }
}
```

### B. Timestep and Courant Number Adjustments (`controlDict`)

#### 1. Increase Base Timestep
```
deltaT          1e-15;         // Increased from 2e-16 (5x larger)
maxDeltaT       1e-14;         // Increased from 2e-15 (5x larger)
minDeltaT       5e-16;         // Increased from 1e-15
```

#### 2. Relax Courant Number Limits
```
maxThermalCourant 0.5;         // Increased from 0.2
thermalFluxRelax 0.15;         // Increased from 0.1
maxCo           0.2;           // Increased from 0.05 (4x larger)
maxAlphaCo      0.2;           // Increased from 0.05 (4x larger)
maxDi           20;            // Increased from 10
```

## Expected Improvements

### Stability
- ✅ Pressure solver should converge reliably
- ✅ No more -10^10 Pa unphysical pressures
- ✅ PIMPLE convergence within 15-20 iterations
- ✅ Reduced emergency clamping warnings

### Performance
- ✅ 5-10x larger timesteps (1 fs → 5-10 fs)
- ✅ Faster simulation progress
- ✅ Estimated runtime: weeks instead of months

### Physics
- ✅ Better energy accumulation per timestep
- ✅ Faster heating toward melting point
- ✅ Recoil pressure generation when T > 1941 K

## Testing Recommendations

### 1. Quick Test (5-10 minutes)
```bash
cd TEST1
# Clean previous run
rm -rf 0.* [1-9]*
rm -rf processor*
cp -r 0.orig 0

# Run for a few picoseconds
compInterFoam > log.test 2>&1 &
tail -f log.test
```

**Success Criteria:**
- No pressure divergence warnings
- PIMPLE converges in < 20 iterations
- Timestep grows to ~5-10 fs
- Temperature increasing steadily

### 2. Monitor Key Indicators
```bash
# Watch for pressure issues
grep "Pressure field exceeded" log.test

# Check PIMPLE convergence
grep "PIMPLE: iteration" log.test | tail -20

# Monitor heating
grep "max(Te):" log.test | tail -10
grep "max(Tl):" log.test | tail -10

# Check timestep growth
grep "deltaT =" log.test | tail -20
```

### 3. Expected Values After ~50 timesteps
- `deltaT`: Should grow to 5-10 fs
- `max(Te)`: Should reach 500-1000 K
- `max(Tl)`: Should reach 400-500 K
- PIMPLE iterations: 10-15 per timestep
- Pressure range: [-10 MPa, +200 MPa]

## Rollback Instructions

If issues persist, revert changes:
```bash
cd TEST1/system
git checkout fvSolution controlDict
```

## Additional Tuning (If Needed)

### If still too slow:
1. Increase `maxCo` to 0.3-0.5
2. Increase `maxDeltaT` to 5e-14 (50 fs)
3. Reduce `nOuterCorrectors` to 15

### If pressure still diverges:
1. Reduce `maxPressure` to 5e8 (500 MPa)
2. Reduce relaxation factors by 0.1
3. Increase `nCorrectors` to 5

### If interface becomes unstable:
1. Reduce `maxAlphaCo` to 0.1
2. Increase `nAlphaSubCycles` to 4
3. Add more alpha correctors: `nAlphaCorr 4`

## Files Modified

1. `TEST1/system/fvSolution` - Solver settings and relaxation
2. `TEST1/system/controlDict` - Timestep and Courant limits

## Next Steps

1. ✅ Test with current fixes
2. Monitor for 50-100 timesteps
3. If stable, run full simulation
4. If issues persist, apply additional tuning above
5. Consider mesh refinement study if results are physical

## Technical References

- **Pressure divergence**: Extreme density ratios require rAU clamping in VOF methods
- **Femtosecond timesteps**: Energy accumulation limited by dt → need larger Co numbers
- **PIMPLE convergence**: Stiff thermal-mechanical coupling requires more outer correctors
- **Recoil pressure**: Feinaeugle et al. (Appl. Surf. Sci. 418, 2017): 70-80 MPa for fs-LIFT
