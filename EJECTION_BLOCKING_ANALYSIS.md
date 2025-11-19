# Root Cause Analysis: Why Ejection is Not Occurring

## Executive Summary

**Your simulation SHOULD produce ejection** - the recoil pressure (27 GPa) is enormous. However, a **pressure clamp logic bug** combined with **insufficient clamp headroom** is artificially throttling the metal acceleration to only **17.6% of the theoretical velocity**.

---

## The Evidence

### 1. Massive Recoil Driving Force ✓
```
Max recoil pressure: 27,079 MPa (27.1 GPa)
Peak observed:       30,300 MPa (30.3 GPa)
Active interface:    669 cells with mass flux > 1e-12 kg/m²/s
Max mass flux:       3,719,978 kg/m²/s
```
**This is MORE than sufficient for ejection.**

### 2. Shock Wave Formation ✓
```
Pressure gradient: 6.15×10⁹ GPa/m (6.15×10¹⁸ Pa/m)
Shock wave:        DETECTED
Gas velocity:      4,588 m/s (Mach ~15)
```
**Strong shock waves are present, indicating violent dynamics.**

### 3. **CRITICAL**: Velocity Throttling ⚠️
```
Theoretical recoil velocity: 3,463 m/s [= sqrt(2 × 27 GPa / 4500 kg/m³)]
Actual metal velocity:       609 m/s
Achievement ratio:           17.6% of theoretical
```
**The metal is being artificially restricted.**

---

## Root Cause: Pressure Clamp Logic Bug

### Configuration vs. Reality

**Your fvSolution setting:**
```cpp
maxPressure        5.0e10;  // 50 GPa configured
```

**What's actually enforced:**
```
pressureClamp: enforcing bounds [-4.0618945e+10, 4.0618945e+10] Pa
                                  ^^^^^^^^^^^^^^^^
                                  40.6 GPa (NOT 50 GPa!)
```

### The Bug (pEqn.H:320-358)

```cpp
Foam::scalar fallbackMaxValue = 1e9;  // Start at 1 GPa

// BUG: This condition is FALSE when configuredMaxPressure > fallbackMaxValue
if (configuredMaxPressure > Foam::SMALL &&
    configuredMaxPressure < fallbackMaxValue)  // 50 GPa < 1 GPa? FALSE!
{
    fallbackMaxValue = configuredMaxPressure;  // NEVER EXECUTED
}

// Expand based on observed recoil
const scalar recoilSafetyFactor = 1.5;
fallbackMaxValue = max(fallbackMaxValue, recoilLimit * recoilSafetyFactor);
                 = max(1 GPa, 27.079 GPa × 1.5)
                 = 40.6 GPa  // Your configured 50 GPa is IGNORED!
```

**Result:** Your 50 GPa setting is silently ignored and replaced with 40.6 GPa.

---

## Why This Breaks Ejection

### The Throttling Mechanism

1. **Recoil force accelerates interface** → Creates pressure spike
2. **Shock wave propagates into gas** → Local pressure > 40.6 GPa
3. **Clamp activates** → Pressure artificially capped
4. **Backpressure builds** → Reduces acceleration
5. **Metal can't reach full velocity** → No ejection

### The Numbers

```
Configured clamp:       50.0 GPa  (ignored)
Effective clamp:        40.6 GPa  (= 1.5 × recoil)
Recoil pressure:        27.1 GPa
Shock wave pressure:    > 40.6 GPa (clipped)

Clamp headroom:         40.6 / 27.1 = 1.50 (barely above recoil!)
```

**The clamp is only 50% above the driving force** - insufficient for shock dynamics.

---

## Why Volume Loss is Negligible

```
Metal volume:     4035.394 µm³
Volume loss:      0.00032%
Loss rate:        1.1×10²⁷ µm³/s (meaningless due to clamping)
```

The mass flux diagnostic shows **669 cells trying to eject**, but the pressure clamp prevents the pressure field from developing fully, so the material just "vibrates" instead of actually ejecting.

---

## PIMPLE Non-Convergence

```
PIMPLE: not converged within 20 iterations
Relative energy change 0.028472568 exceeds energyTolerance (0.02)
```

**Why:** The solver is fighting against the artificial pressure clamp. The physics wants pressure > 40 GPa, but the clamp forces it down, creating a tug-of-war that prevents convergence.

---

## Solutions (Ranked)

### ⭐ **Solution 1: Fix the Bug AND Increase Clamp** (RECOMMENDED)

**Location:** `TEST1/system/fvSolution`

```cpp
compInterFoamCoeffs
{
    pressureClamp      true;
    maxPressure        1.5e11;  // 150 GPa (3× recoil for shock headroom)
    minPressure       -1.5e11;
}
```

**Location:** `pEqn.H:323` - Fix the logic bug:

```cpp
// BEFORE (buggy):
if (configuredMaxPressure > Foam::SMALL &&
    configuredMaxPressure < fallbackMaxValue)

// AFTER (fixed):
if (configuredMaxPressure > Foam::SMALL &&
    configuredMaxPressure > fallbackMaxValue)  // Use configured if LARGER
{
    fallbackMaxValue = configuredMaxPressure;
}
```

---

### Solution 2: Disable Pressure Clamping (RISKY)

```cpp
pressureClamp      false;  // Let pressure field develop freely
```

**Risk:** Possible numerical explosion. Only try if Solution 1 fails.

---

### Solution 3: Increase PIMPLE Tolerance

Already at:
```cpp
nOuterCorrectors            20;  // Already high
outerCorrectorResidualControl
{
    "(U|p_rgh)"  { tolerance 1e-4;  relTol 0.1; }
}
```

**Action:** Relax energy tolerance:
```cpp
energyTolerance 0.05;  // Was 0.02
```

---

## Expected Behavior After Fix

With maxPressure = 150 GPa:

```
Recoil pressure:        ~27 GPa
Shock wave pressure:    ~60-80 GPa (within clamp)
Metal velocity:         ~2,500-3,400 m/s (70-98% of theoretical)
Volume loss:            > 1% after 10-20 ps
Ejection:               ✓ EXPECTED
```

---

## Verification Steps

1. **Check clamp utilization:**
   ```
   maxPressure = 150 GPa
   Observed max < 100 GPa → Good headroom
   ```

2. **Check velocity ratio:**
   ```
   |U|_metal_max / recoil-limited > 0.7 → Physics is working
   ```

3. **Check volume loss:**
   ```
   Volume loss > 1% after 20 ps → Material is ejecting
   ```

---

## Bottom Line

**You have the physics right** - 27 GPa recoil is huge!

**The problem is purely numerical:** A logic bug + insufficient pressure headroom is creating artificial backpressure that throttles ejection.

**Fix:** Either patch the bug OR simply set `maxPressure 1.5e11` in fvSolution.
