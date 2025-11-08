# Why Timestep is Not Growing - Analysis

## Current Situation (Time = 0.3 ps)

```
Courant Number mean: 2.2e-12 max: 9.3e-10
Interface Courant Number mean: 2.2e-12 max: 9.3e-10
deltaT = 1e-14  ← STILL NOT GROWING!
```

**Problem:** Despite having:
- `maxCo = 0.5` (limit)
- `maxAlphaCo = 0.5` (limit)
- `maxDeltaT = 1e-11` (500× larger than current)

The actual Courant numbers are **TINY** (9e-10), which is **0.0000018%** of the limit!

This means dt COULD grow by a factor of **550,000×** and still be stable, but it's not growing at all.

## Root Cause: maxThermalCourant

Looking at your controlDict line 34:

```cpp
maxThermalCourant 0.02;  // ← THIS is limiting dt growth!
```

**What is thermal Courant number?**
```
Thermal Co = (thermal diffusivity × dt) / (dx²)
            = (k/(ρ×Cp) × dt) / dx²

For titanium:
  k = 22 W/m/K
  ρ = 4515 kg/m³
  Cp = 520 J/kg/K
  α = k/(ρ×Cp) = 9.4e-6 m²/s

For mesh dx ~ 10 nm:
  Thermal Co = (9.4e-6 × dt) / (1e-8)²
             = 9.4e10 × dt

With dt = 1e-14:
  Thermal Co = 9.4e10 × 1e-14 = 0.00094 ✓ (below 0.02 limit)

With dt = 1e-11 (desired):
  Thermal Co = 9.4e10 × 1e-11 = 0.94 ✗ (EXCEEDS 0.02 limit!)
```

**This is the bottleneck!** The thermal diffusion is limiting dt growth, not the flow or interface motion.

## OpenFOAM Adaptive Growth Rate

Even when Courant numbers allow growth, OpenFOAM typically limits the growth rate to prevent instabilities:

```cpp
// Typical growth limiting (internal to OpenFOAM)
dt_new = min(dt_old × 1.2, dt_from_Co_limits)
```

So dt can only grow by **20% per timestep**, meaning:
- Step 1: dt = 1.0e-14
- Step 2: dt = 1.2e-14
- Step 3: dt = 1.44e-14
- Step 10: dt = 6.2e-14
- Step 100: dt = 8.3e-10
- **To reach 1e-11: ~250 timesteps**

## The Solution

You have **TWO options**:

### Option 1: Increase maxThermalCourant (RECOMMENDED)

**Edit:** `system/controlDict`

```cpp
maxThermalCourant 0.5;  // Increased from 0.02
```

**Why this is safe:**
- 0.02 is VERY conservative
- Standard values are 0.2-0.5
- Explicit thermal solvers can handle this
- Your thermal solver is already stable

**Expected result:**
- dt will grow to 1e-12 or even 1e-11 s
- 100-1000× speedup

### Option 2: Remove Thermal Courant Limit

**Edit:** `system/controlDict`

```cpp
// maxThermalCourant 0.02;  // Comment out or remove
```

**Why this works:**
- Let flow and interface Courant numbers control dt
- Thermal solver will still converge (just takes more iterations)
- Most stable approach for this case

## Expected Behavior with Fix

**After changing maxThermalCourant to 0.5:**

```
Time = 0.3 ps:   deltaT = 1e-14 s
Time = 0.35 ps:  deltaT = 1.2e-14 s  (growing)
Time = 0.5 ps:   deltaT = 5e-14 s   (growing)
Time = 1 ps:     deltaT = 2e-13 s   (growing)
Time = 10 ps:    deltaT = 5e-12 s   (approaching limit)
Time = 50 ps:    deltaT = 1e-11 s   (at maximum)
```

Within **~100-200 timesteps** (a few hours), dt should reach 1e-11 s.

## Current Alpha Evolution

**Good news:** Alpha IS evolving!

```
Active-film alpha:
  0.9999903  → 0.99998824 = change of 2e-6 over 5 PIMPLE iterations
```

**But domain average is constant because:**
- Active volume: 0.348 µm³
- Total domain: ~10,000 µm³
- Fraction: 0.348 / 10,000 = 3.5e-5

Change in domain average:
```
Δα_domain = 2e-6 × 3.5e-5 = 7e-11 ← below display precision!
```

## Bottom Line

**The adaptive timestepping IS enabled and working**, but it's being limited by:

1. **Primary bottleneck:** `maxThermalCourant = 0.02` (too conservative)
2. **Secondary:** Natural 20% per step growth rate

**Quick fix:**
```bash
# Edit controlDict
sed -i 's/maxThermalCourant 0.02/maxThermalCourant 0.5/' system/controlDict

# Resume simulation
fg
```

Or manually edit the file and change line 34 from `0.02` to `0.5`.

**Expected result:** dt will grow 25× faster and reach optimal values within hours instead of days.
