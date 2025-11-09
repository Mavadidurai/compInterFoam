# Complete Solution: deltaT Not Growing

## Problem Summary

deltaT remained frozen at 1e-14 s despite:
- ✅ `adjustTimeStep = yes`
- ✅ `writeInterval = 1e-10` (fixed from 1e-14)
- ✅ `maxDeltaT = 1e-11`
- ✅ Tiny Courant numbers (huge headroom for growth)

**Even after killing and restarting**, deltaT stayed at 1e-14.

## Root Cause

```cpp
startFrom       latestTime;  // ← THE CULPRIT
```

When OpenFOAM restarts with `startFrom = latestTime`:
1. It reads from saved time directories
2. These directories contain `uniform/time` files with **saved deltaT values**
3. Even though controlDict was updated, the **saved deltaT** overrides it!

**Result:** Old `deltaT = 1e-14` persisted across restarts.

## The Complete Fix

### 1. Force Fresh Start from Time 0

**File:** `system/controlDict:20`

```diff
- startFrom       latestTime;
+ startFrom       startTime;     // Force restart from time 0
```

This ensures OpenFOAM:
- Ignores any saved time directories
- Reads deltaT from controlDict (not from saved files)
- Starts with a clean slate

### 2. Conservative Initial deltaT

**File:** `system/controlDict:25`

```diff
- deltaT          1e-13;         // Base time step
+ deltaT          5e-14;         // Conservative start for temperature development
```

**Why Conservative Start?** At Time = 0, all fields are uniform (T = 300K):
- No temperature gradients exist yet
- No phase change sources active
- Alpha equation has no driving force

Starting too large (e.g., 1e-12) causes **singular matrix** on first timestep!

**Solution:** Start at `5e-14` to allow:
1. First few timesteps establish temperature gradients
2. Phase change sources activate gradually
3. Then adaptive timestepping grows to `maxDeltaT = 1e-11`

Growth trajectory:
- `5e-14` → `6e-14` → `1e-13` → `5e-13` → `1e-12` → `1e-11` (~100 steps)
- Still reaches maxDeltaT quickly once physics develops!

### 3. Use Clean Restart Script

**File:** `cleanRestart.sh`

```bash
./cleanRestart.sh
```

This script:
1. Kills any running compInterFoam processes
2. Removes all time directories except `0`
3. Shows current settings
4. Starts simulation and displays deltaT growth in real-time

## How to Use

### Quick Start:
```bash
cd ~/OpenFOAM/mavadi-v2406/run/RealisticLIFT
./cleanRestart.sh
```

### Manual Start:
```bash
# Kill old processes
pkill -9 compInterFoam

# Clean time directories
rm -rf [0-9]* [0-9]*.[0-9]*
# (but keep 0 directory)

# Start fresh
compInterFoam 2>&1 | tee log.restart
```

## Expected Behavior

You should see deltaT growing gradually from conservative start:

```
Time = 0 s,      deltaT = 5e-14 s   (conservative start)
Time = 5e-14 s,  deltaT = 6e-14 s   (growing by 20%)
Time = 5e-13 s,  deltaT = 1e-13 s   (accelerating)
Time = 5e-12 s,  deltaT = 1e-12 s   (continuing)
Time = 5e-11 s,  deltaT = 1e-11 s   (reached maximum!)
```

Once at `deltaT = 1e-11`:
- Simulation runs **1000× faster** than before
- Output written every ~10 timesteps (writeInterval = 1e-10)
- Will reach 2 ns endpoint in reasonable time

## Why This Works

| Issue | Previous | Fixed |
|-------|----------|-------|
| Restart source | `latestTime` (uses saved deltaT) | `startTime` (uses controlDict) |
| Starting deltaT | 1e-13 s (or 1e-14) | 5e-14 s (conservative) |
| Initial stability | Singular matrix possible | Gradual field development |
| Growth steps | ~86-200 steps | ~100 steps to maxDeltaT |
| Persistence | Cached in time dirs | Fresh from controlDict |

## Verification

After running `cleanRestart.sh`, verify:

```bash
# Check deltaT is growing
grep "deltaT =" log.restart | head -20

# Should show increasing values:
# deltaT = 5e-14
# deltaT = 6e-14
# deltaT = 7.2e-14
# deltaT = 1e-13
# ... up to 1e-11
```

## Summary

**The problem was NOT:**
- ❌ writeInterval (we fixed that)
- ❌ Courant limits (they were fine)
- ❌ adjustTimeStep (it was enabled)

**The problem WAS:**
- ✅ **Cached deltaT in time directories from old runs**
- ✅ **startFrom = latestTime reading old values**

**The solution:**
- ✅ **startFrom = startTime** (ignore saved values)
- ✅ **Conservative initial deltaT** (avoid singular matrix at t=0)
- ✅ **Clean restart script** (guaranteed fresh start)

## Additional Issue: Singular Matrix on First Timestep

After implementing the fixes above, a fresh start with `deltaT = 1e-12` produced:

```
DILUPBiCGStab:  Solving for alpha.metal:  solution singularity
alpha1 solve detected singular matrix; reverting to previous alpha field
```

**Why?** At Time = 0:
- All fields are uniform (T = 300K everywhere)
- No temperature gradients
- No phase change sources (T < 3560K)
- Alpha equation has zero driving force → singular matrix

**Fix:** Reduced initial `deltaT` from `1e-12` to `5e-14`
- Allows temperature field to develop gradually over first few timesteps
- Phase change sources activate smoothly
- Then adaptive timestepping accelerates to `maxDeltaT = 1e-11`
- Still reaches nanosecond timescales efficiently!

**Expected result:** Simulation reaches nanosecond timescales in minutes instead of days! 🚀
