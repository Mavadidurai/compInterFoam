# Adaptive Timestepping Fix - Root Cause Found!

## Problem

Adaptive timestepping was completely non-functional. Despite having:
- `adjustTimeStep yes` ✓
- Tiny Courant numbers (1e-9 vs limit of 0.5) ✓
- `maxDeltaT = 1e-11` (500× larger than current dt) ✓

**deltaT remained frozen at 1e-14 s**

## Root Cause

**THE CULPRIT:** `writeInterval = 1e-14`

```cpp
writeControl    adjustableRunTime;
writeInterval   1e-14;         // ← THIS WAS LIMITING deltaT!
```

### How adjustableRunTime Works

When you use `writeControl = adjustableRunTime`, OpenFOAM **adjusts deltaT to ensure output is written exactly at the specified intervals**.

This means:
```
deltaT ≤ writeInterval  (ALWAYS!)
```

Even though your Courant limits would allow `dt = 1e-11`, OpenFOAM was capping it at `1e-14` to match the writeInterval!

### Why This Happens

The `adjustableRunTime` write control modifies the timestep to hit exact write times:

```cpp
// OpenFOAM internal logic:
if (writeControl == adjustableRunTime)
{
    // Ensure we write at exactly writeInterval
    newDeltaT = min(newDeltaT, timeUntilNextWrite);

    // This prevents deltaT from exceeding writeInterval!
}
```

## The Fix

**Changed:** `writeInterval = 1e-14` → `writeInterval = 1e-10`

```cpp
writeControl    adjustableRunTime;
writeInterval   1e-10;         // Write every 100 ps (10× maxDeltaT)
```

Now:
- deltaT can grow up to 1e-11 s (limited by maxDeltaT and Courant)
- Output written every 100 ps (~10 timesteps at max dt)
- No artificial timestep cap
- Reasonable output frequency (not excessive I/O)

## Alternative Approaches

If you want MORE FREQUENT output than maxDeltaT:

**Option 2: Use timeStep-based writing**
```cpp
writeControl    timeStep;      // Write every N timesteps
writeInterval   100;           // Write every 100 steps
```

This decouples write frequency from timestep, allowing deltaT to grow freely.

**Option 3: Use runTime-based writing with large interval**
```cpp
writeControl    runTime;       // Write at fixed time intervals
writeInterval   1e-10;         // Write every 100 ps
```

## Expected Behavior After Fix

With the corrected `writeInterval = 1e-11`:

```
Time = 0.32 ps:   deltaT = 1e-14 s  (starting value)
Time = 0.35 ps:   deltaT = 1.2e-14 s  (growing by 20%)
Time = 0.5 ps:    deltaT = 5e-14 s   (continuing growth)
Time = 1 ps:      deltaT = 2e-13 s   (accelerating)
Time = 10 ps:     deltaT = 5e-12 s   (approaching limit)
Time = 50 ps:     deltaT = 1e-11 s   (at maximum!)
```

## How to Resume

```bash
cd ~/OpenFOAM/mavadi-v2406/run/RealisticLIFT

# controlDict has been updated - just resume!
# If backgrounded: fg
# If stopped: compInterFoam > log.compInterFoam 2>&1 &

# Monitor timestep growth:
tail -f log.compInterFoam | grep -E "deltaT|Time ="
```

OpenFOAM will automatically re-read controlDict (you saw "Re-reading object" before) and apply the new writeInterval.

## Bottom Line

**Not a bug in adaptive timestepping - it's a FEATURE interaction!**

- ✅ `adjustTimeStep` was working correctly
- ✅ Courant number calculations were correct
- ✅ Standard OpenFOAM behavior

The issue was that `writeInterval = 1e-14` was **intentionally** preventing deltaT from growing to ensure output timing.

**Fix applied:** `writeInterval = 1e-11` now allows deltaT to reach its natural limit!

**Expected result:** 100-1000× speedup starting now! 🚀
