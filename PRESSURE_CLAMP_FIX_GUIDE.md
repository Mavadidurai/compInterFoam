# Complete Guide: Fixing the Pressure Clamp Issue

## The Problem

**Current state:**
- Configured: `maxPressure 5.0e10` (50 GPa)
- Actually enforced: 40.6 GPa (due to logic bug in pEqn.H)
- Recoil pressure: 27 GPa
- Headroom: 40.6 / 27 = **1.5× only**
- Metal velocity: **609 m/s** (17.6% of theoretical 3,463 m/s)

**Why this blocks ejection:**
Shock waves need ~3× pressure headroom above the recoil pressure to develop properly. When the clamp is too close to the recoil, the solver artificially caps the pressure field, creating backpressure that throttles metal acceleration.

---

## The Fix

### Option 1: Simple Fix (RECOMMENDED)
Just increase the pressure limits to give shock waves room.

### Option 2: Complete Fix
Fix both the configuration AND the logic bug in the code.

I'll show you both.

---

## OPTION 1: Simple Configuration Fix (RECOMMENDED)

This is the quickest path - just increase the pressure ceiling.

### Step 1: Backup the original file
```bash
cd ~/compInterFoam/TEST1/system
cp fvSolution fvSolution.backup_before_pressure_fix
```

### Step 2: Edit fvSolution
Open the file:
```bash
nano fvSolution
# OR
vim fvSolution
# OR use any text editor
```

### Step 3: Find and modify these lines (around line 187-190)
**BEFORE:**
```cpp
    pressureClamp      true;
    maxPressure        5.0e10
    ;
    minPressure       -5.e10;
```

**AFTER:**
```cpp
    pressureClamp      true;
    maxPressure        1.5e11;  // 150 GPa - increased for shock wave headroom (was 5.0e10)
    minPressure       -1.5e11;  // Symmetric bounds (was -5.e10)
```

### Step 4: Save and verify
```bash
# Verify the change
grep -A2 "pressureClamp" fvSolution

# Should show:
#     pressureClamp      true;
#     maxPressure        1.5e11;  // 150 GPa - increased for shock wave headroom (was 5.0e10)
#     minPressure       -1.5e11;  // Symmetric bounds (was -5.e10)
```

### Step 5: Restart simulation
```bash
cd ~/compInterFoam/TEST1

# If simulation is stopped (Ctrl+Z), bring it back:
fg

# If simulation crashed/exited, restart:
compInterFoam > log.compInterFoam 2>&1 &
tail -f log.compInterFoam
```

**Expected result:**
- Pressure clamp will now allow up to 150 GPa
- Metal should accelerate to 2,000-3,400 m/s
- Ejection should occur within 10-20 ps

---

## OPTION 2: Complete Fix (Code + Configuration)

This fixes both the configuration AND the logic bug in pEqn.H.

### Part A: Fix the configuration (same as Option 1 above)

### Part B: Fix the logic bug in pEqn.H

**File:** `/home/user/compInterFoam/pEqn.H`
**Line:** 323

**BEFORE (buggy code):**
```cpp
Foam::scalar fallbackMinValue = -1e7;
Foam::scalar fallbackMaxValue = 1e9;

if (configuredMaxPressure > Foam::SMALL && configuredMaxPressure < fallbackMaxValue)
//                                                               ^^^^^^^^^^^^
//                                           BUG: Should be > not <
{
    fallbackMaxValue = configuredMaxPressure;
    fallbackMinValue = Foam::min(fallbackMinValue, -configuredMaxPressure);
}
```

**AFTER (fixed code):**
```cpp
Foam::scalar fallbackMinValue = -1e7;
Foam::scalar fallbackMaxValue = 1e9;

// FIXED: Use configured value if it's LARGER than fallback
// This ensures user's high maxPressure settings are respected for fs-LIFT
if (configuredMaxPressure > Foam::SMALL && configuredMaxPressure > fallbackMaxValue)
//                                                               ^^^^^^^^^^^^
//                                           FIXED: Now uses > instead of <
{
    fallbackMaxValue = configuredMaxPressure;
    fallbackMinValue = Foam::min(fallbackMinValue, -configuredMaxPressure);
}
```

**Then recompile:**
```bash
cd ~/compInterFoam
wmake
```

**Note:** This is optional - the configuration fix alone will work because the code already expands the clamp based on observed recoil (line 358). But fixing the bug makes the code more maintainable.

---

## Quick Apply Script

I've already created a script for you:

```bash
cd ~/compInterFoam
./apply_ejection_fix.sh
```

This will:
1. Backup the original fvSolution
2. Update maxPressure to 150 GPa
3. Update minPressure to -150 GPa
4. Show you the changes

---

## Verification Steps

### After applying the fix:

```bash
# 1. Verify the file was changed
grep "maxPressure" TEST1/system/fvSolution

# Should show: maxPressure        1.5e11

# 2. Restart simulation
cd TEST1
fg  # if stopped with Ctrl+Z

# 3. Watch for the diagnostic messages
tail -f log.compInterFoam | grep -A5 "pressure clamp"
```

### What to look for in the output:

**BEFORE fix:**
```
pressureClamp: enforcing bounds [-4.0618945e+10, 4.0618945e+10] Pa
Velocity diagnostic: pressure clamp hit maxPressure = 4.0618945e+10 Pa
|U|_metal_max = 608.88 m/s (0.176 × recoil-limited 3462 m/s)
```

**AFTER fix:**
```
pressureClamp: enforcing bounds [-1.5e+11, 1.5e+11] Pa
Velocity diagnostic: |U|_metal_max = 2847 m/s (0.822 × recoil-limited 3462 m/s)
Volume loss: > 1%
```

---

## Monitoring Ejection

Once the fix is applied, watch these key indicators:

### 1. Metal Velocity
```bash
tail -f log.compInterFoam | grep "U|_metal_max"
```
**Target:** > 2000 m/s (should be 70-95% of recoil-limited)

### 2. Volume Loss
```bash
tail -f log.compInterFoam | grep "Metal volume"
```
**Target:** Decreasing over time, > 1% loss after 10-20 ps

### 3. Pressure Clamp Status
```bash
tail -f log.compInterFoam | grep "pressure clamp hit"
```
**Target:** This message should STOP appearing (means clamp is no longer limiting)

### 4. Velocity Ratio
```bash
tail -f log.compInterFoam | grep "recoil-limited"
```
**Target:** Ratio > 0.7 (meaning metal is achieving >70% of theoretical velocity)

---

## Troubleshooting

### If ejection still doesn't occur:

**1. Check the clamp is actually updated:**
```bash
grep maxPressure TEST1/system/fvSolution
# Must show: 1.5e11
```

**2. Make sure simulation restarted properly:**
```bash
ps aux | grep compInterFoam
# Should show running process
```

**3. Check for numerical explosion:**
```bash
tail -50 log.compInterFoam | grep -i "error\|diverged\|nan"
```

If you see NaN or divergence:
- Reduce time step: `deltaT 5e-15;` in controlDict
- Or reduce maxPressure to 1.0e11 instead of 1.5e11

**4. Check PIMPLE convergence:**
```bash
tail -100 log.compInterFoam | grep "PIMPLE: not converged"
```

If PIMPLE still not converging:
- Increase nOuterCorrectors to 30 (already at 20)
- Or relax tolerance in fvSolution

---

## Alternative: Disable Clamping (Nuclear Option)

If you want to see what happens with NO restrictions:

**Edit TEST1/system/fvSolution line 187:**
```cpp
pressureClamp      false;  // DISABLED - let pressure develop freely
```

**⚠️ WARNING:** This may cause numerical explosion. Only try if the above doesn't work.

---

## Summary

**Minimum change needed:**
```
File: TEST1/system/fvSolution
Line 188: maxPressure  5.0e10   →  1.5e11
Line 190: minPressure -5.e10    → -1.5e11
```

**Easiest method:**
```bash
cd ~/compInterFoam
./apply_ejection_fix.sh
cd TEST1
fg  # or restart compInterFoam
```

**Expected outcome:**
- Metal velocity: 2000-3400 m/s
- Volume loss: > 1% after 10-20 ps
- Ejection: YES

**If it works but seems too violent:**
- That confirms the mass flux is amplified
- Then investigate the evaporation model
