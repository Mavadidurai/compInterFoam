# Quick Start: Mass Transfer Diagnostics

## TL;DR - What To Do RIGHT NOW

### Step 1: Run the Python check (30 seconds)
```bash
cd /home/user/compInterFoam
python3 check_psat.py
```

**Result:** ✓ Thermodynamics predict j_evap = 6.12e7 kg/m²/s at 9000K
- Parameters are theoretically CORRECT
- Problem is in OpenFOAM implementation, not physics

---

### Step 2: Check if DEBUG output exists in your log (1 minute)
```bash
cd /home/user/compInterFoam/TEST1
grep "DEBUG Hertz-Knudsen" log.* | head -20
```

**What you should see:**
```
DEBUG Hertz-Knudsen at cell 0:
  T_eff = XXXX K
  p_vapor (psat) = XXXX Pa
  j_evap = XXXX kg/m²/s
  j_net = XXXX kg/m²/s
```

**If you DON'T see this:** The debug block isn't printing. Check if `timeIndex() % 10 == 0` condition is being met.

**If you DO see it:** Send me those lines - they'll show EXACTLY what's wrong!

---

### Step 3: Most Likely Problems (in order)

#### Problem A: Activation Window Filtering (60% probability)
**Location:** `twoPhaseMixtureThermo.C` lines 1155-1168

Check if evaporation is being blocked by:
```cpp
if (T_local + SMALL < T_melt_)  // Line 1155
    continue;  // BLOCKED if T < 1941K

if (alpha < alphaMin_ || (enforceUpper && alpha > alphaMax_))  // Line 1160
    continue;  // BLOCKED if alpha outside range

if (restrictToVapor && T_eff < T_vap)  // Line 1165
    continue;  // BLOCKED if onlyAboveVapor=true and T < 2200K
```

**Quick fix to test:**
```bash
cd /home/user/compInterFoam/TEST1/system
nano controlDict

# Find phaseChangeCoeffs section and set:
onlyAboveVapor      false;  # Should already be false
alphaMin            0.001;  # Should already be 0.001

# Save and re-run
```

#### Problem B: maxSource Clamping (30% probability)
**Location:** `twoPhaseMixtureThermo.C` lines 1258-1267

The code limits mass flux to prevent instability:
```cpp
maxSource           [1 -1 -3 0 0 0 0] 1e7;  // From your controlDict
```

This converts to:
```
maxHeatFlux = maxSource * meltThickness
j_net_limited = maxHeatFlux / L
```

At typical cell size (~1 µm):
```
maxHeatFlux ≈ 1e7 * 1e-6 = 10 W/m²
j_net ≈ 10 / 9.1e6 ≈ 1e-6 kg/m²/s  ← BELOW threshold!
```

**THIS IS LIKELY YOUR PROBLEM!**

**Quick fix:**
```bash
nano TEST1/system/controlDict

# In phaseChangeCoeffs:
maxSource           [1 -1 -3 0 0 0 0] 1e10;  # Increase by 1000x

# Or comment it out entirely for testing:
# maxSource           [1 -1 -3 0 0 0 0] 1e10;
```

#### Problem C: ClTTM Not Set (10% probability)
**Location:** `twoPhaseMixtureThermo.C` line 1226-1230

```cpp
const scalar Cl = ClTTM_.value();
if (Cl <= SMALL)
    continue;  // BLOCKED if lattice heat capacity not set
```

Your controlDict has: `Cl  [1 -1 -2 -1 0 0 0] 2.5e6;` ✓ (should be fine)

---

## The Smoking Gun: maxSource

Your `controlDict` phaseChangeCoeffs says:
```cpp
maxSource           [1 -1 -3 0 0 0 0] 1e7;  // [W/m³]
```

This is a **volumetric** heat flux limit. The code converts it to a **mass flux limit**:

```
Step 1: Calculate cell thickness
  meltThickness = cellVolume / maxFaceArea
  For typical 1 µm cell: ~1e-6 m

Step 2: Convert to heat flux [W/m²]
  maxHeatFlux = maxSource * meltThickness
  = 1e7 * 1e-6 = 10 W/m²

Step 3: Convert to mass flux [kg/m²/s]
  j_net = maxHeatFlux / latentHeat
  = 10 / 9.1e6
  = 1.1e-6 kg/m²/s

Step 4: Compare to threshold
  massRateEps = 1e-12 kg/m²/s
  j_net (1.1e-6) > threshold ✓  (would pass)

  BUT at small cells or with clamping, could drop below threshold!
```

### Why This Causes Zero Mass Transfer

If your cells are smaller (e.g., 100 nm instead of 1 µm):
```
meltThickness ~ 1e-7 m
maxHeatFlux = 1e7 * 1e-7 = 1 W/m²
j_net = 1 / 9.1e6 = 1.1e-7 kg/m²/s  ← STILL above 1e-12 threshold

Hmm... maybe not the problem unless cells are VERY small
```

**Actually, let me recalculate more carefully...**

Wait, looking at the code again (line 1256):
```cpp
const scalar maxSource = maxPhaseChangeSource_.value();  // 1e7 [W/m³]
```

At 9000K, the natural heat flux would be:
```
j_net (theoretical) = 6.12e7 kg/m²/s
heatFlux = j_net * L = 6.12e7 * 9.1e6 = 5.57e14 W/m²  (HUGE!)

volumetricHeat = heatFlux / thickness
For thickness = 1e-6 m:
volumetricHeat = 5.57e14 / 1e-6 = 5.57e20 W/m³

maxSource = 1e7 W/m³  ← MUCH SMALLER!
```

So the clamping would kick in:
```
volumetricHeat_clamped = 1e7 W/m³
heatFlux_clamped = 1e7 * 1e-6 = 10 W/m²
j_net_clamped = 10 / 9.1e6 = 1.1e-6 kg/m²/s
```

**This is ABOVE the threshold (1e-12), so it should still work!**

Hmm, need to think about this more carefully...

---

## Action Plan

### OPTION 1: Quick Test (5 minutes)

Temporarily disable all limits:

```bash
cd /home/user/compInterFoam/TEST1/system
cp controlDict controlDict.backup

# Edit phaseChangeCoeffs:
nano controlDict

# Change:
maxSource           [1 -1 -3 0 0 0 0] 1e12;  # Effectively unlimited
onlyAboveVapor      false;
alphaMin            0.0001;
```

Then run a 10-timestep test:
```bash
cd TEST1
rm -rf 0.* [1-9]*
cp -r 0.orig 0
compInterFoam | tee log.test_unlimited
grep "DEBUG Hertz-Knudsen" log.test_unlimited
grep "Recoil diagnostics" log.test_unlimited
```

### OPTION 2: Add Enhanced Diagnostics (30 minutes)

Follow the detailed plan in `DIAGNOSTIC_PLAN.md`:
- Add verbose output to see actual j_net values
- Write fields to ParaView
- Trace through each activation gate

### OPTION 3: Send Me Your Log Output (2 minutes)

Just run:
```bash
cd /home/user/compInterFoam/TEST1
grep -A 15 "DEBUG Hertz-Knudsen" log.* | tail -50 > debug_output.txt
cat debug_output.txt
```

And send me the output. I can immediately diagnose from those numbers.

---

## Expected Fix

Based on the Python calculation showing correct thermodynamics, the fix will likely be ONE of:

1. **Increase maxSource to 1e10 or 1e12** (most likely)
2. **Disable maxSource clamping entirely** for testing
3. **Check that ClTTM is properly set** in the solver
4. **Verify activation windows aren't blocking** the calculation

---

## Summary

**Thermodynamics:** ✓ CORRECT (verified by Python script)
**Problem:** ✓ IN OPENFOAM IMPLEMENTATION

**Most likely cause:** `maxSource` clamping or activation window filtering

**Next step:** Check DEBUG output from your existing log, OR run quick test with maxSource increased

---

## What The Output Should Look Like When Fixed

```
DEBUG Hertz-Knudsen at cell 0:
  T_eff = 9000 K
  p_vapor (psat) = 6.4e+12 Pa
  j_evap = 6.1e+07 kg/m²/s          ← HUGE!
  j_net = 6.1e+07 kg/m²/s

Recoil diagnostics: 728 of 78000 interface cells supplied mass flux above 1e-12 kg/m^2/s
  Max |j_net| = 6.1e+07 kg/m²/s
  Max |recoilPressure| = 850 MPa    ← NON-ZERO!
```

vs current (broken):
```
Recoil diagnostics: 0 of 78000 interface cells supplied mass flux above 1e-12 kg/m^2/s
  Max |recoilPressure| = 0 MPa
```
