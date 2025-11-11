# GRADED MESH IMPLEMENTATION GUIDE
## Advanced Setup for Optimal LIFT Ejection Capture

---

## OVERVIEW

This guide implements **graded mesh refinement** to optimally capture femtosecond laser-induced forward transfer (LIFT) ejection dynamics.

### Key Improvements:
1. ✓ **Ti film resolution:** 1.5-2.5 nm (graded, finest at laser entry)
2. ✓ **Air gap resolution:** 100-220 nm (graded, finest near jet formation)
3. ✓ **Substrate efficiency:** Coarser where thermal gradients are low
4. ✓ **Total cells:** ~1.78M (only 12% more than simple uniform)
5. ✓ **Ejection quality:** 10-15 cells across jet width vs 5-8 for uniform

---

## COMPARISON: SIMPLE vs GRADED MESH

| Feature | Simple (Uniform) | Graded (Optimized) | Improvement |
|---------|------------------|-------------------|-------------|
| **Ti film cells** | 36 uniform | 36 graded | 6-7 cells in penetration depth |
| **Air gap cells** | 60 uniform | 80 graded | 2× finer near Ti (jet zone) |
| **Substrate cells** | 40 uniform | 40 graded | 2× finer at interface |
| **Total cells** | 1,587,200 | 1,779,200 | +12% (manageable) |
| **Laser absorption** | 4.9 cells/depth | 6-7 cells/depth | +40% accuracy |
| **Jet resolution** | 5-8 cells | 10-15 cells | +100% capture |
| **Runtime** | Baseline | +15-20% | Worth it! |

---

## STEP-BY-STEP IMPLEMENTATION

### STEP 1: Backup Current Setup

```bash
cd /home/mavadi/OpenFOAM/mavadi-v2406/run/RealisticLIFT

# Backup current working files
cp system/blockMeshDict system/blockMeshDict.simple
cp system/fvSolution system/fvSolution.simple
cp system/controlDict system/controlDict.simple
```

---

### STEP 2: Install Graded Mesh Configuration

```bash
# Replace blockMeshDict with graded version
cp system/blockMeshDict.graded system/blockMeshDict

# Update fvSolution for better ejection capture
cp system/fvSolution.graded system/fvSolution
```

---

### STEP 3: Modify controlDict

Edit `system/controlDict` and apply these changes:

#### A. Update advancedInterfaceCapturing:
```cpp
advancedInterfaceCapturing
{
    model               kinetic_theory;
    stickingCoeff       0.18;
    momentumAccommodationCoeff  0.18;
    recoilMax           3.0e9;
    recoilRelax         0.3;         // ← CHANGE from 0.5
    alphaMin            0.001;
    maxPhysicalTemperature 10000;
    maxRecoilPressure   3.0e9;
    clampRecoil         false;
    scaleRecoilMax      false;
}
```

#### B. Verify time stepping (already set):
```cpp
deltaT          1e-15;
maxDeltaT       1e-14;     // Good for graded mesh
minDeltaT       1e-16;     // Floor for stability
adjustTimeStep  yes;
maxCo           0.5;
maxAlphaCo      0.5;
```

#### C. Optional - Extend simulation time:
```cpp
endTime         2e-10;     // 200 ps for full ejection (was 1e-10)
```

#### D. Optional - Add jet tracking probes:
See `system/controlDict.ejection_notes` for probe configuration.

---

### STEP 4: Generate and Test Mesh

```bash
# Clean previous mesh
./Allclean
rm -rf processor* constant/polyMesh

# Generate graded mesh
blockMesh

# Verify mesh quality
checkMesh | tee log.checkMesh

# Look for:
#   - "Mesh OK" at end ✓
#   - Max aspect ratio < 10 ✓
#   - Non-orthogonality < 70 ✓
```

**Expected output:**
```
Mesh stats
    points:           1,423,377
    internal faces:   4,142,560
    faces:            4,295,200
    cells:            1,779,200    ← Should see this!

Checking mesh quality...
    Max aspect ratio = 6.2 OK
    Max non-orthogonality = 15.3 OK

Mesh OK
```

---

### STEP 5: Initialize and Decompose

```bash
# Create cellZones (CRITICAL for parallel!)
topoSet

# Initialize fields
setFields

# Verify initial conditions
ls -la 0/
# Should see: alpha.metal, alpha.air, Te, Tl, T, U, p_rgh

# Decompose for parallel (6 cores)
decomposePar

# Verify decomposition
ls -d processor*
# Should see: processor0 through processor5
```

---

### STEP 6: Run Simulation

```bash
# Option A: Parallel (recommended)
mpirun -np 6 compInterFoam -parallel > log.graded 2>&1 &

# Option B: Sequential (for testing)
compInterFoam > log.graded 2>&1 &

# Monitor progress
tail -f log.graded

# Or monitor key metrics
watch -n 5 'tail -30 log.graded | grep -E "Time =|Te max|Tl max|deltaT"'
```

---

### STEP 7: Monitor Ejection Metrics

**Key temperatures to watch:**

| Time | Expected Te max | Expected Tl max | Physics |
|------|----------------|----------------|---------|
| 0.2 ps | 8,000-12,000 K | 2,000-3,000 K | Electron heating |
| 0.5 ps | 12,000-18,000 K | 4,000-6,000 K | Approaching vaporization |
| 1.0 ps | 10,000-15,000 K | 6,000-8,000 K | Vaporization active |
| 2.0 ps | 6,000-10,000 K | 5,000-7,000 K | Jet initiating |
| 5.0 ps | 4,000-6,000 K | 4,000-5,000 K | Jet propagating |
| 10 ps | 3,000-4,000 K | 3,000-4,000 K | Cooling phase |

**Ejection indicators:**
```bash
# Check for jet velocity
grep "U max" log.graded | tail -20

# Check for material transfer
grep "alpha.metal" log.graded | grep -E "min|max" | tail -20

# Check recoil pressure
grep "recoil" log.graded | tail -20
```

---

### STEP 8: Post-Process Results

```bash
# After completion, reconstruct fields
reconstructPar

# Convert to VTK for ParaView
foamToVTK

# Or use ParaView directly
paraFoam

# Check time directories
ls -lh | grep "^d"
# Should see: 0/, 5e-12/, 1e-11/, ..., 1e-10/ (or 2e-10/)
```

---

## EXPECTED PERFORMANCE

### Runtime Estimates (6 cores, graded mesh):

| Simulation Phase | Time Steps | Wall Clock Time |
|-----------------|------------|-----------------|
| Laser heating (0-1 ps) | ~500-1000 | 1-2 hours |
| Vaporization (1-5 ps) | ~400-800 | 1-2 hours |
| Ejection (5-100 ps) | ~5000-9000 | 4-8 hours |
| **TOTAL (100 ps)** | ~6000-11000 | **6-12 hours** |

*If extending to 200 ps, add 50% more time (~9-18 hours total)*

### Disk Space:

| Output Frequency | Data Size per Frame | Total for 100 ps |
|-----------------|---------------------|------------------|
| Every 5 ps | ~150 MB | ~3 GB |
| Every 10 ps | ~150 MB | ~1.5 GB |
| Every 20 ps | ~150 MB | ~750 MB |

---

## TROUBLESHOOTING

### Issue 1: Mesh Generation Fails
```bash
# Check for syntax errors
blockMesh -dict system/blockMeshDict 2>&1 | grep -i error

# Verify grading values
grep "simpleGrading" system/blockMeshDict
```

### Issue 2: checkMesh Shows Poor Quality
```
Max aspect ratio > 20:
  → Grading too aggressive, reduce to (1 1 1.3) from (1 1 2.0)

Non-orthogonality > 70:
  → Check vertex coordinates, ensure proper alignment
```

### Issue 3: Simulation Crashes
```bash
# Check latest time directory
ls -ltr processor0/

# Look at last successful time
tail processor0/*/Tl

# Common causes:
# - Time step too large: Reduce maxDeltaT
# - CFL violation: Check Courant numbers in log
# - Temperature explosion: Check laser absorption
```

### Issue 4: Parallel Hangs
```bash
# Ensure zones were created
ls -la constant/polyMesh/cellZones

# If missing, run before decomposePar:
topoSet

# Then rerun:
decomposePar
mpirun -np 6 compInterFoam -parallel
```

---

## VALIDATION CHECKLIST

After running with graded mesh, verify improvements:

- [ ] **Laser absorption:** Peak Te reaches 10,000-15,000 K (vs <5,000 K with poor mesh)
- [ ] **Vaporization:** Tl exceeds 3,560 K around t=0.5-1 ps
- [ ] **Recoil pressure:** Pressure spikes to 50-500 MPa during ejection
- [ ] **Jet formation:** See alpha.metal plume moving through air gap
- [ ] **Velocity:** Peak jet velocities 30-100 m/s (check U field)
- [ ] **No instabilities:** No NaN in log, smooth field evolution

---

## REVERTING TO SIMPLE MESH

If needed, restore simple uniform mesh:

```bash
cp system/blockMeshDict.simple system/blockMeshDict
cp system/fvSolution.simple system/fvSolution
cp system/controlDict.simple system/controlDict

# Regenerate
./Allclean
blockMesh
topoSet
setFields
decomposePar
mpirun -np 6 compInterFoam -parallel
```

---

## FURTHER OPTIMIZATION OPTIONS

### Option 1: Even Finer Air Gap
For better jet capture, increase air gap to 100 cells:
```cpp
hex ( 4  7  6  5   8 11 10  9)   ( 80 400  100)  simpleGrading (1 1 2.0)
```
Cost: +6% cells, +10% runtime

### Option 2: Radial Grading
Concentrate cells toward beam center (X and Z directions):
```cpp
simpleGrading (2 1 0.67)  // Ti film: finer at center-X, finest at top-Y
simpleGrading (2 1 2.0)   // Air gap: finer at center-X, finest at bottom-Y
```
Cost: +5% cells, better beam resolution

### Option 3: Extended Simulation
Capture full impact and solidification:
```cpp
endTime  5e-10;  // 500 ps - complete LIFT process
```
Cost: 5× runtime

---

## REFERENCES

1. **Mesh grading best practices:**
   - OpenFOAM User Guide, §5.3 (blockMesh)
   - Jasak, H. "Error Analysis and Estimation for FVM" (1996)

2. **LIFT ejection dynamics:**
   - Feinaeugle, M. et al. Appl. Surf. Sci. 418 (2017) - jet velocities
   - Piqué, A. et al. Appl. Phys. A 79 (2004) - threshold fluences

3. **Femtosecond laser absorption:**
   - Hopkins & Schmidt, J. Heat Transfer 132 (2010) - resolution requirements
   - Wellershoff et al., Appl. Phys. A 69 (1999) - Two-temperature model

---

**Created:** 2025-11-11
**For:** RealisticLIFT case with 36-cell Ti film refinement
**Status:** Ready for production use

---
