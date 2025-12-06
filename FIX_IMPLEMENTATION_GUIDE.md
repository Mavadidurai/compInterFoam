# Alpha-Weighted Recoil Pressure Fix - Implementation Guide

**Date**: 2025-12-06
**Issue**: Recoil pressure accelerating vapor (1123 m/s) instead of metal (0.009 m/s)
**Solution**: Alpha-weight recoil pressure in momentum equation

---

## Quick Summary

### The Problem:
```
Current behavior:
  - Metal avg velocity: 0.009 m/s  ❌ (should be 30-100 m/s)
  - Vapor velocity:     1123 m/s   ⚠️ (should be 500-1000 m/s)
  - Recoil pressure:    71.37 MPa  ✅ (correct value)
  - PIMPLE:            Not converging (20 iterations) ❌
```

**Root cause**: Recoil pressure gradient `-∇P` accelerates **low-density vapor** much more than **high-density metal** because `acceleration ∝ 1/ρ`.

### The Fix:
Apply alpha-weighting to recoil pressure: `P_total = p_rgh + α·P_recoil`

This makes the force **proportional to α**:
- Metal (α=1): Full recoil force ✅
- Vapor (α=0): No recoil force ✅
- Interface: Pressure traction (physical) ✅

---

## Files Modified

1. **UEqn.H** - Line 111 (one line change + documentation)

---

## Step-by-Step Implementation

### Step 1: Verify Current State

```bash
cd /home/user/compInterFoam
grep -n "totalPressure += \*recoilContribution" UEqn.H
```

**Expected output**: Should show line 111 (if already fixed, will show alpha1)

### Step 2: Apply the Fix

The fix has already been applied to `UEqn.H`. Verify it:

```bash
grep -A 2 "totalPressure +=" UEqn.H | grep alpha1
```

**Expected**: Should see `totalPressure += alpha1 * (*recoilContribution);`

### Step 3: Recompile

```bash
cd /home/user/compInterFoam
wmake
```

**Watch for**:
- ✅ Compilation succeeds
- ❌ Any errors (check syntax if errors occur)

**Expected output** (final lines):
```
'-Wl,-rpath,/opt/openfoam/...'
-o /home/user/OpenFOAM/mavadi-v2406/platforms/.../compInterFoam
wmake: done
```

### Step 4: Test the Fix

```bash
cd TEST3

# Option A: Resume stopped simulation
fg

# Option B: Restart from latest time
compInterFoam

# Option C: Restart from specific time
compInterFoam -time 13.28e-12
```

### Step 5: Monitor Key Metrics

Watch the log output for these changes:

#### ✅ Success Indicators:

**Velocity Diagnostic:**
```
Before fix:
  Velocity diagnostic: |U|_metal_max = 1123 m/s (6.32× recoil-limited 177.8 m/s)
  avg(|U|) (metal): 0.0093 m/s

After fix (expected):
  Velocity diagnostic: |U|_metal_max = 30-100 m/s (0.2-0.6× recoil-limited)
  avg(|U|) (metal): 5-20 m/s  ← Should increase significantly!
```

**PIMPLE Convergence:**
```
Before fix:
  PIMPLE: not converged within 20 iterations

After fix (expected):
  PIMPLE: converged in 5-10 iterations
```

**Phase Velocities:**
```
Before fix:
  max(|U|): 1123.69 m/s
  avg(|U|) (metal): 0.0093 m/s

After fix (expected):
  max(|U|): 500-1000 m/s  (vapor thermal expansion)
  avg(|U|) (metal): 30-100 m/s  (metal jet)
```

#### ❌ Potential Issues:

**Interface Oscillations:**
```
Symptom: Rapid pressure fluctuations
Cause: P·∇α term creating large interface forces
Solution: Add pEqn flux correction (Phase 3 below)
```

**Checkerboard Pattern:**
```
Symptom: Alternating high/low pressure cells
Cause: Pressure-velocity decoupling
Solution: Increase PIMPLE iterations, reduce relaxation factors
```

---

## Validation Tests

### Test 1: Check Metal Velocity Increase

```bash
# While simulation is running, monitor velocity:
tail -f log.compInterFoam | grep "avg(|U|) (metal)"
```

**Success**: Value increases from 0.009 to 5-50 m/s range

### Test 2: Check PIMPLE Convergence

```bash
tail -f log.compInterFoam | grep "PIMPLE:"
```

**Success**: See "PIMPLE: iteration X" where X < 15 (not hitting 20 limit)

### Test 3: Visualize in ParaView

```bash
paraFoam &
```

**Load latest timestep, check**:
1. **Velocity field** (color by U magnitude):
   - High velocity (30-100 m/s) in metal regions (α > 0.5)
   - Lower velocity in vapor (α < 0.1)
   - Upward jet from film surface

2. **Alpha field** (color by alpha.metal):
   - Metal jet/droplets visible (α > 0.5)
   - Traveling away from surface

3. **Pressure field** (color by p_rgh):
   - Smooth distribution
   - No oscillations or checkerboard

---

## Understanding the Physics

### Before Fix (Wrong):

```
Recoil pressure: P_recoil = 71 MPa
Gradient: ∇P = constant everywhere
Force: F = -∇P
Acceleration: a = F/ρ

Metal (ρ=4500 kg/m³):
  a = ∇P / 4500 = small → U = 0.009 m/s ❌

Vapor (ρ=1 kg/m³):
  a = ∇P / 1 = huge → U = 1123 m/s ❌
```

### After Fix (Correct):

```
Effective pressure: P_eff = α · P_recoil
Gradient: ∇(α·P) = P·∇α + α·∇P
Force: F = -∇(α·P)
Acceleration: a = F/ρ

Metal (α=1, ∇α≈0):
  F = -∇P
  a = ∇P / 4500 = proper → U = 30-100 m/s ✅

Vapor (α=0, ∇α≈0):
  F = 0
  a = 0 (only thermal) → U = thermal only ✅

Interface (α=0.5, ∇α large):
  F = -P·∇α (interface traction)
  a = interface force / ρ_mixed ✅
```

---

## Configuration Checks

### Vapor Velocity Limiter

**File**: `TEST3/system/fvSolution`
**Setting**: `maxVaporVelocity` (not currently set)
**Default**: 2000 m/s (from UEqn.H line 450)

**Your current value**: 2000 m/s (default)
**Status**: ✅ OK - 1123 m/s is below limit

**To adjust** (if needed):
```cpp
// Add to compInterFoamCoeffs section in fvSolution:
maxVaporVelocity   1500;  // Limit vapor to 1500 m/s
```

### Global Velocity Limiter

**File**: `TEST3/system/fvSolution`
**Setting**: `maxVelocity` (line 168)
**Current value**: 0 (DISABLED)

**Status**: ✅ CORRECT - limiter should be disabled for physics

**Do NOT change** this - it's correctly disabled to allow proper physics.

---

## Troubleshooting

### Issue 1: Compilation Fails

**Error**: `error: 'alpha1' was not declared in this scope`

**Cause**: alpha1 not visible in UEqn.H scope

**Fix**: Check that alpha1 is defined before line 111. Look for:
```cpp
const volScalarField& alpha1 = mixture.alpha1();
```

If not found, add near top of UEqn.H after line 10.

---

### Issue 2: Simulation Crashes After Fix

**Symptom**: `Floating point exception` or `FOAM FATAL ERROR`

**Likely cause**: Large interface forces from P·∇α term

**Fix Option A - Reduce Recoil Strength** (temporary):
```cpp
// In TEST3/system/controlDict, advancedInterfaceCapturing section:
maxRecoilPressure  5e7;  // Reduce from default (71 MPa → 50 MPa)
```

**Fix Option B - Add pEqn Flux Correction** (advanced):
See "Phase 3: Advanced Fix" section below.

---

### Issue 3: PIMPLE Still Not Converging

**Symptom**: Still seeing "PIMPLE: not converged within 20 iterations"

**Fix**: Increase PIMPLE correctors in `fvSolution`:
```cpp
PIMPLE
{
    nOuterCorrectors  30;  // Increase from 20
    nCorrectors       10;  // Increase from 8
}
```

And/or increase relaxation in `compInterFoamCoeffs`:
```cpp
URelaxationFactor  0.5;  // Increase from 0.3
pRelaxationFactor  0.4;  // Increase from 0.2
```

---

## Phase 3: Advanced Fix (If Needed)

If you see **interface oscillations** or **checkerboard patterns**, add flux correction in pEqn.H:

### Step 1: Open pEqn.H
```bash
nano pEqn.H
```

### Step 2: Find line 87 (currently commented out section)

Look for:
```cpp
// Recoil traction removed from phig - recoil pressure handled via p_rgh only
```

### Step 3: Add this code after line 87:

```cpp
// Flux correction to cancel P·∇α spurious interface force
if (recoilPressurePtr)
{
    surfaceScalarField alphaf = fvc::interpolate(alpha1);
    surfaceScalarField recoilPf = fvc::interpolate(*recoilContribution);
    surfaceScalarField gradAlphaf = fvc::snGrad(alpha1);

    // This cancels the P·∇α term from ∇(α·P) in UEqn
    // Leaves only α·∇P (metal acceleration)
    phig -= recoilPf * gradAlphaf * rAUf * mesh.magSf();

    if (verbose && Foam::Pstream::master())
    {
        Info<< "    Applied recoil flux correction to cancel P·∇α interface force"
            << nl;
    }
}
```

### Step 4: Recompile and test
```bash
wmake && cd TEST3 && compInterFoam
```

**This is ONLY needed if you see problems!** Try UEqn-only fix first.

---

## Expected Results Summary

| Metric | Before | After | Target |
|--------|--------|-------|--------|
| **Metal avg velocity** | 0.009 m/s | 30-100 m/s | ✅ Jet formation |
| **Max velocity** | 1123 m/s | 500-1000 m/s | ✅ Thermal expansion |
| **Velocity ratio** | 6.3× recoil | 0.2-0.6× recoil | ✅ Physical |
| **PIMPLE iters** | 20 (fail) | 5-10 | ✅ Converged |
| **Metal loss rate** | 0.0098 µm³/ps | Similar/higher | ✅ Ablation |
| **Alpha.metal** | 0.40212 | ~0.4020-0.4021 | ✅ Slight drop OK |

---

## Success Criteria

After applying the fix, you should see:

✅ **Metal velocity increases** from 0.009 → 30-100 m/s
✅ **Vapor velocity decreases** from 1123 → 500-1000 m/s
✅ **PIMPLE converges** in <10 iterations
✅ **Velocity ratio** drops from 6.3× to 0.2-0.6× recoil-limited
✅ **Jet formation** visible in ParaView (metal parcels ejecting)
✅ **No crashes** or numerical instabilities

---

## References

1. Feinaeugle, M., et al. (2017). "Time-resolved imaging of flyer dynamics for laser-induced forward transfer of metals". *Applied Surface Science*, 418, 72-79.

2. Brackbill, J.U., Kothe, D.B., Zemach, C. (1992). "A continuum method for modeling surface tension". *Journal of Computational Physics*, 100(2), 335-354.

3. Tryggvason, G., Scardovelli, R., Zaleski, S. (2011). *Direct Numerical Simulations of Gas-Liquid Multiphase Flows*. Cambridge University Press.

---

## Contact & Support

**Analysis documents**:
- `alpha_weighting_analysis.md` - Detailed UEqn vs pEqn comparison
- `alpha_metal_analysis.md` - Why alpha.metal appears constant
- `simulation_analysis.md` - Initial simulation health check
- `phase_transition_analysis.md` - Vaporization phase analysis

**Patch file**: `alpha_weighted_recoil_fix.patch`

**Implementation**: Already applied to `UEqn.H` line 111
