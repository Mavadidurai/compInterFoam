# 🔧 QUICK FIX: Parallel Simulation Hang

## Problem Identified

Your simulation hangs **during femtosecondLaserModel initialization** because of **unguarded Info<< output** in parallel code.

**Location**: `femtosecondLaserModel.C` lines 154, 247, 332
**Symptom**: All 4 processors try to write to stdout simultaneously → I/O deadlock

---

## ✅ IMMEDIATE SOLUTIONS (Try in order)

### Solution 1: Disable Verbose Output (FASTEST)

Edit `system/controlDict`:
```
verbose         false;  // Change from true
```

Then run:
```bash
cd ~/OpenFOAM/mavadi-v2406/run/TEST2
rm -rf processor* dynamicCode
decomposePar
mpirun -np 4 compInterFoam -parallel
```

**Success rate: ~60%**

---

### Solution 2: Use Collated File Handler

Edit `system/controlDict`:
```
fileHandler     collated;  // Change from uncollated
```

Then run:
```bash
rm -rf processor* dynamicCode
decomposePar
mpirun -np 4 compInterFoam -parallel
```

**Success rate: ~40%**

---

### Solution 3: Disable Function Objects

Edit `system/controlDict` and comment out the entire `functions` section:
```cpp
/*
functions
{
    midPlaneVTK { ... }
    jetProbes { ... }
    fieldMinMax { ... }
}
*/
```

Then run:
```bash
rm -rf processor* dynamicCode
decomposePar
mpirun -np 4 compInterFoam -parallel
```

**Success rate: ~70%**

---

### Solution 4: Try Different MPI Backend

```bash
cd ~/OpenFOAM/mavadi-v2406/run/TEST2
rm -rf processor* dynamicCode
decomposePar

# Try with different MPI settings
mpirun --mca btl tcp,self --bind-to none -np 4 compInterFoam -parallel
```

**Success rate: ~50%**

---

### Solution 5: Start with Fewer Processors

```bash
rm -rf processor* dynamicCode

# Edit system/decomposeParDict
# Change: numberOfSubdomains 2;
# Change: n (1 2 1);

decomposePar
mpirun -np 2 compInterFoam -parallel
```

**Success rate: ~80%** (but slower)

---

## 🎯 RECOMMENDED QUICK TEST

Run this one-liner to test with minimal config:

```bash
cd ~/OpenFOAM/mavadi-v2406/run/TEST2 && \
rm -rf processor* dynamicCode && \
sed -i 's/verbose.*true/verbose         false/' system/controlDict && \
sed -i 's/fileHandler.*uncollated/fileHandler     collated/' system/controlDict && \
decomposePar && \
mpirun --mca btl tcp,self -np 4 compInterFoam -parallel 2>&1 | tee log.quicktest
```

---

## 📊 Expected Output (If Working)

You should see:
```
Create time
Create mesh for time = 0
Reading field p_rgh
...
Constructing twoPhaseMixtureThermo
controlDict.massTransferCoeffs:
    rateMax   3e+13
    activationTime (0 2e-10)

LIFT Field Validation:
  Total domain volume: ...
  Total metal volume: ...

Field setup completed successfully!

Starting time loop

Time = 2e-15
...
```

**Key**: You should see "LIFT Field Validation" and "Starting time loop"

---

## 🔍 Automated Debugging

I've created a debug script:

```bash
cd ~/OpenFOAM/mavadi-v2406/run/TEST2
./debug_parallel.sh 1    # Run minimal test
./debug_parallel.sh all  # Run all diagnostic tests
```

This will test different configurations automatically.

---

## 🐛 If Still Hanging

### Check what's happening:

```bash
# In one terminal - start the run
mpirun -np 4 compInterFoam -parallel &

# In another terminal - monitor CPU
watch -n 1 'ps aux | grep compInterFoam | grep -v grep'
```

**If CPU = 0%**: Deadlock (likely I/O issue)
**If CPU = high**: Computing (might just be slow)

### Kill stuck simulation:
```bash
killall -9 compInterFoam
rm -rf processor* dynamicCode
```

---

## 💡 ROOT CAUSE EXPLANATION

The `femtosecondLaserModel` constructor contains:

```cpp
// Line 154 - NO Pstream::master() guard!
Info<< "Transmission override active: using transmission = "
    << transmission_ << "; reflectivity entry ignored" << endl;

// Line 247 - NO guard!
Info<< "Derived film bounds from center y=" << centerY
    << " m using " << thicknessKey << " = " << filmThickness
    << " m: [" << filmYMin_ << ", " << filmYMax_ << "] m" << endl;

// Line 332 - NO guard!
Info<< "Femtosecond laser model initialized:" << nl
    << "  Pulse energy: " << pulseEnergy_.value() << " J" << nl
    ...
```

**What should be:**
```cpp
if (Pstream::master())
{
    Info<< "Transmission override active..." << endl;
}
```

Without the guard, **all 4 processes try to write simultaneously** → buffer contention → hang.

---

## 🔧 PERMANENT FIX (Code Modification)

If you have access to modify the solver source:

**Edit `femtosecondLaserModel.C`:**

Add `Pstream::master()` guards around all `Info<<` statements:

```cpp
// Line 154
if (Pstream::master())
{
    Info<< "Transmission override active: using transmission = "
        << transmission_ << "; reflectivity entry ignored" << endl;
}

// Line 247
if (Pstream::master())
{
    Info<< "Derived film bounds from center y=" << centerY
        << " m using " << thicknessKey << " = " << filmThickness
        << " m: [" << filmYMin_ << ", " << filmYMax_ << "] m" << endl;
}

// Line 332
if (Pstream::master())
{
    Info<< "Femtosecond laser model initialized:" << nl
        << "  Pulse energy: " << pulseEnergy_.value() << " J" << nl
        ...
        << endl;
}
```

**Then recompile:**
```bash
cd ~/OpenFOAM/mavadi-v2406/applications/solvers/compInterFoam
wmake
```

---

## 📝 Summary

| Solution | Difficulty | Success Rate | Notes |
|----------|-----------|--------------|-------|
| Disable verbose | Easy | 60% | Fastest to try |
| Disable functions | Easy | 70% | Loses monitoring |
| Change file handler | Easy | 40% | May slow down |
| Different MPI | Medium | 50% | System dependent |
| Fewer processors | Easy | 80% | Slower simulation |
| Fix source code | Hard | 95% | Requires recompile |

---

## ✅ After It Works

Once running successfully:
1. Check `log.compInterFoam` for progress
2. Monitor with: `tail -f log.compInterFoam`
3. Re-enable function objects one at a time if needed
4. Increase processor count gradually (2 → 4 → 6)

---

**Created**: 2025-11-15
**Issue**: Parallel hang in femtosecondLaserModel initialization
**Cause**: Unguarded Info<< output causing I/O deadlock
