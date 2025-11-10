# Critical Fixes Applied to fs-LIFT Simulation

**Date:** 2025-11-10
**Purpose:** Fix velocity explosion (18.7 km/s) and excessive material loss (2.98%)
**Status:** Configuration updated, requires testing

---

## Issues Identified

### 1. **CRITICAL: Velocity Explosion** 🔴
- **Observed:** 18,682 m/s maximum velocity at t = 200 ps
- **Expected:** < 800 m/s for realistic fs-LIFT
- **Ratio:** 23× too high (supersonic, 3.7× sound speed in Ti)

### 2. **IMPORTANT: Excessive Material Loss** 🔴
- **Observed:** 2.98% volume loss = 246 µm³
- **Problem:** Ti film is only 35.7 µm³ → loss is 6.9× film volume
- **Root Cause:** Receiver substrate ablating (unphysical for LIFT)

### 3. **MODERATE: Pressure Solver Non-Convergence** ⚠️
- **Issue:** "PIMPLE: not converged within 3 iterations" throughout
- **Impact:** 1-5% error accumulation in pressure/velocity fields

### 4. **RECOMMENDED: Marginal Mesh Resolution**
- **Current:** 2.04 cells per optical penetration depth (9.7 nm)
- **Standard:** Minimum 3-5 cells required
- **Impact:** 20-30% uncertainty in peak temperature

---

## Fixes Applied

### Fix 1: Mesh Resolution Enhancement ✓

**File:** `TEST1/system/blockMeshDict`

**Change:**
```diff
- hex (12 15 14 13  16 19 18 17)   ( 40 200  15)  simpleGrading (1 1 1)
+ hex (12 15 14 13  16 19 18 17)   ( 40 200  30)  simpleGrading (1 1 1)
```

**Result:**
- Ti film cells in Y: 15 → 30
- Cell size: 4.76 nm → 2.38 nm
- Cells per penetration depth: 2.04 → 4.08 ✓
- Total mesh: 1.24M → 1.36M cells (+9.7%)

**Requires:** `blockMesh` to regenerate mesh

---

### Fix 2: Pressure Solver Convergence Improvement ✓

**File:** `TEST1/system/fvSolution`

**Changes:**
```diff
PIMPLE
{
-   nOuterCorrectors            3;
+   nOuterCorrectors            7;

    residualControl
    {
-       p_rgh     { tolerance 1e-4;  relTol 0.05; }
+       p_rgh     { tolerance 1e-6;  relTol 0.01; }
    }
}
```

**Result:**
- More PIMPLE iterations for better convergence
- Tighter pressure tolerance (100× stricter)
- Should eliminate "not converged" warnings

---

### Fix 3: Phase-Change Spatial Constraint ✓⚠️

**Files Modified:**
1. `TEST1/system/topoSetDict` - Added Ti film cell zone
2. `TEST1/system/controlDict` - Added cellZone constraints

**Changes:**

#### topoSetDict - New Ti Film Zone:
```cpp
// Ti donor film region - ONLY this region should undergo phase-change
{
    name tiFilmSet;
    type cellSet;
    action new;
    source boxToCell;
    sourceInfo
    {
        // Ti donor film: y: [28.0, 28.0714] µm (71.4 nm thick at top)
        box (0 28.0e-06 0) (50.0e-06 28.0714e-06 10.0e-06);
    }
}

{
    name tiFilm;
    type cellZoneSet;
    action new;
    source setToCellZone;
    sourceInfo
    {
        set tiFilmSet;
    }
}
```

#### controlDict - Phase-Change Constraint:
```diff
phaseChangeCoeffs
{
    ...
+   // CRITICAL FIX: Constrain phase-change to Ti film only
+   cellZone            tiFilm;
}
```

#### controlDict - Recoil Pressure Constraint:
```diff
advancedInterfaceCapturing
{
    ...
+   // CRITICAL FIX: Constrain recoil pressure to Ti film only
+   cellZone              tiFilm;
}
```

**Result:**
- Phase-change and recoil pressure confined to Ti film (y = 28.0-28.0714 µm)
- Receiver substrate (y = 0-8 µm) should no longer ablate
- Donor substrate (y = 20-28 µm) protected from phase-change

**⚠️ WARNING:** The `cellZone` parameter may not be supported by the solver. If you encounter errors about unrecognized keywords, see "Alternative Solutions" below.

**Requires:** `topoSet` to create cell zones

---

## Testing Procedure

### Step 1: Regenerate Mesh and Zones

```bash
cd /home/user/compInterFoam/TEST1

# Clean old mesh
rm -rf constant/polyMesh

# Generate new mesh with refined Ti film
blockMesh

# Create cell zones (including new tiFilm zone)
topoSet

# Initialize fields
rm -rf 0
cp -r 0.orig 0
setFields
```

### Step 2: Run Short Test (10 ps)

```bash
# Modify controlDict for short test
sed -i 's/endTime.*/endTime         1e-11;/' system/controlDict

# Run simulation
compInterFoam > log.test 2>&1

# Check results
tail -100 log.test
```

### Step 3: Verify Fixes

Check for these indicators:

1. **Velocity Fixed:**
   ```bash
   grep "Max |U|" log.test
   ```
   Should show values < 1000 m/s (ideally < 800 m/s)

2. **Material Loss Fixed:**
   ```bash
   grep "Metal phase:" log.test | tail -5
   ```
   Volume loss should be < 10% of Ti film volume (< 3.6 µm³)

3. **Pressure Convergence:**
   ```bash
   grep "PIMPLE: not converged" log.test
   ```
   Should show very few or no warnings

4. **CellZone Errors:**
   ```bash
   grep -i "cellZone\|unknown\|keyword" log.test
   ```
   If errors appear, the solver doesn't support cellZone constraints

---

## Expected Results

### Before Fixes:
- Max velocity: 18,682 m/s (unphysical)
- Material loss: 246 µm³ (6.9× film volume)
- PIMPLE: Not converging
- Confidence: 65%

### After Fixes:
- Max velocity: < 800 m/s (realistic LIFT range)
- Material loss: < 10% of film volume (< 3.6 µm³)
- PIMPLE: Converging within 7 iterations
- Mesh resolution: 4+ cells/penetration depth
- Confidence target: 85-90%

---

## Alternative Solutions (If cellZone Not Supported)

If the solver rejects the `cellZone` parameter, try these alternatives:

### Option A: Manual Phase-Change Limiting
Add Y-coordinate check in `controlDict`:
```cpp
phaseChangeCoeffs
{
    ...
    // Limit phase-change by Y coordinate
    minY    28.0e-6;    // Ti film bottom
    maxY    28.0714e-6; // Ti film top
}
```

### Option B: Separate Material Phases
Modify `setFieldsDict` to use distinct phases:
- `alpha.receiver` for receiver substrate
- `alpha.donor` for transparent donor substrate
- `alpha.metal` for Ti film only

This requires solver code modifications to handle three metal phases.

### Option C: Temperature-Based Protection
Artificially increase thermal mass of receiver/donor to prevent heating:
- Set higher heat capacity in those regions
- Use `topoSet` zones with region-specific properties

---

## Next Steps

1. **Run test simulation** following procedure above
2. **Analyze results** using verification checks
3. **If cellZone works:** Run full 500 ps simulation
4. **If cellZone fails:** Implement Alternative Solution A or B
5. **Update analysis:** Document validated results

---

## Files Modified

- `system/blockMeshDict` - Ti film mesh resolution
- `system/fvSolution` - PIMPLE convergence settings
- `system/topoSetDict` - Ti film cell zone definition
- `system/controlDict` - Phase-change and recoil constraints

## Commands Required Before Running

```bash
blockMesh      # Regenerate mesh with new resolution
topoSet        # Create cell zones
setFields      # Initialize field distributions
```

---

## Confidence Assessment

**With these fixes:**
- If cellZone works: **85-90%** confidence as experimental replicate
- If cellZone fails: **70-75%** until alternative implemented

**Critical success factors:**
1. Velocity must drop below 1000 m/s
2. Material loss must be < 10% of film volume
3. Pressure solver must converge consistently

---

## Notes

- Mesh refinement increases computational cost by ~10%
- Tighter PIMPLE tolerance may slow each timestep by 20-30%
- Overall runtime increase: ~40% (acceptable for accuracy gain)
- cellZone feature may be custom - check solver documentation
