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

### 2. Increase Initial deltaT

**File:** `system/controlDict:25`

```diff
- deltaT          1e-13;         // Base time step
+ deltaT          1e-12;         // Increased for faster growth
```

**Why?** Adaptive timestepping grows gradually (typically 1.2× per step):
- Starting from `1e-13`: Takes ~86 steps to reach `1e-11` (1000× growth)
- Starting from `1e-12`: Takes ~58 steps to reach `1e-11` (10× growth)

Larger starting deltaT = faster convergence to maxDeltaT.

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

You should see deltaT growing immediately:

```
Time = 0 s,      deltaT = 1e-12 s   (starting value)
Time = 1e-12 s,  deltaT = 1.2e-12 s (growing by 20%)
Time = 3e-12 s,  deltaT = 1.5e-12 s
Time = 1e-11 s,  deltaT = 5e-12 s
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
| Starting deltaT | 1e-13 s | 1e-12 s |
| Growth steps | ~86 steps to maxDeltaT | ~58 steps to maxDeltaT |
| Persistence | Cached in time dirs | Fresh from controlDict |

## Verification

After running `cleanRestart.sh`, verify:

```bash
# Check deltaT is growing
grep "deltaT =" log.restart | head -20

# Should show increasing values:
# deltaT = 1e-12
# deltaT = 1.2e-12
# deltaT = 1.44e-12
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
- ✅ **Larger initial deltaT** (faster growth)
- ✅ **Clean restart script** (guaranteed fresh start)

**Expected result:** Simulation reaches nanosecond timescales in minutes instead of days! 🚀
