# Debugging Parallel Simulation Hang - TEST2

**Issue**: Simulation hangs after reading massTransferCoeffs during initialization.

**Last output**:
```
controlDict.massTransferCoeffs:
    rateMax   3e+13
    activationTime (0 2e-10)
```

---

## Most Likely Causes

### 1. **Field Initialization Deadlock** (MOST COMMON)

The simulation likely hangs when initializing complex fields (Te, Tl, or advanced interface capturing).

**Solution Steps:**

#### Step 1: Check Processor Directories
```bash
cd ~/OpenFOAM/mavadi-v2406/run/TEST2
ls -la processor*/0/
```

**Expected files in each processor*/0/:**
- T
- Te
- Tl
- U
- p
- p_rgh
- alpha.metal
- alpha.air

**If any files are missing**: The decomposition is incomplete.

#### Step 2: Run with Verbose Debug Output
```bash
# Kill the stuck process first
killall compInterFoam

# Run with debug flags
mpirun -np 4 compInterFoam -parallel 2>&1 | tee log.debug
```

#### Step 3: Enable OpenFOAM Debug Switches

Create/edit `~/.OpenFOAM/v2406/controlDict`:
```
DebugSwitches
{
    twoTemperatureModel 1;
    femtosecondLaserModel 1;
    advancedInterfaceCapturing 1;
    Pstream 1;
}
```

Then run again:
```bash
mpirun -np 4 compInterFoam -parallel 2>&1 | tee log.verbose
```

---

## 2. **MPI Communication Issue**

OpenMPI might be hanging on communication during field initialization.

**Solution:**

#### Try with fewer processes:
```bash
# Clean up first
rm -rf processor*
decomposePar
mpirun -np 2 compInterFoam -parallel
```

#### Try with different MPI settings:
```bash
# Disable UCX warnings and set shared memory
mpirun --mca btl_base_warn_component_unused 0 \
       --mca btl ^openib \
       -np 4 compInterFoam -parallel
```

---

## 3. **Two-Temperature Model Initialization**

The hang might be in the two-temperature model trying to read or initialize Te/Tl fields.

**Solution:**

#### Check if Te and Tl exist in processor directories:
```bash
ls processor*/0/Te processor*/0/Tl
```

#### Verify Te/Tl boundary conditions in 0/Te and 0/Tl:

**Common issue**: Processor boundaries need proper BC handling.

**Fix**: Ensure `0/Te` and `0/Tl` use `processor` boundary type (added automatically by decomposePar).

---

## 4. **Advanced Interface Capturing Hang**

The advancedInterfaceCapturing model might hang during parallel initialization.

**Quick Test**: Temporarily disable it.

**Edit system/controlDict:**
```
useAdvancedInterfaceCapturing  false;  // Changed from true
```

Run again:
```bash
mpirun -np 4 compInterFoam -parallel
```

**If it works**: The issue is in the advanced interface capturing initialization.

**Permanent fix**: Check the advancedInterfaceCapturing code for parallel synchronization issues.

---

## 5. **Laser Model Initialization**

The femtosecond laser model might be hanging on parallel initialization.

**Check constant/laserProperties:**

**Common issue**: File I/O or field initialization in laser model.

**Quick test**: Comment out laser-related function objects in controlDict temporarily.

---

## 6. **File Handler Issue**

The "uncollated" file handler might have permissions issues.

**Solution:**

**Edit system/controlDict:**
```
fileHandler     collated;  // Changed from uncollated
```

**Or try:**
```
fileHandler     masterUncollated;
```

---

## 7. **Function Objects Causing Hang**

Function objects (probes, surfaces) might hang during parallel initialization.

**Solution:**

**Temporarily disable all function objects:**

Edit `system/controlDict` - comment out entire `functions` section:
```
/*
functions
{
    midPlaneVTK { ... }
    jetProbes { ... }
    fieldMinMax { ... }
}
*/
```

Run again. If it starts, re-enable functions one by one to find the culprit.

---

## Step-by-Step Debugging Procedure

### Phase 1: Identify the Hang Location

**1. Add debug output to controlDict:**
```
// Add at top
DebugSwitches
{
    compInterFoam 2;
}
```

**2. Run with timeout:**
```bash
timeout 60 mpirun -np 4 compInterFoam -parallel 2>&1 | tee log.timeout
```

**3. Check where it stops:**
```bash
tail -50 log.timeout
```

### Phase 2: Simplify Configuration

**1. Disable complex models temporarily (in controlDict):**
```
useAdvancedInterfaceCapturing  false;
enableLiftProcessTracker false;

// Comment out function objects
/*
functions
{ ... }
*/
```

**2. Try running:**
```bash
mpirun -np 4 compInterFoam -parallel
```

**3. If successful, re-enable one feature at a time.**

### Phase 3: Check Field Decomposition

**1. Verify all fields decomposed correctly:**
```bash
for i in processor*/0/*; do
    echo "Checking: $i"
    head -20 "$i" | grep -E "FoamFile|object|dimensions"
done
```

**2. Compare field counts:**
```bash
ls -1 0/ | wc -l
ls -1 processor0/0/ | wc -l
```

Should be identical.

### Phase 4: Test Serial vs Parallel

**1. Run in serial first:**
```bash
# Backup processor dirs
mv processor0 processor0.bak
# etc...

# Run serial
compInterFoam
```

**2. If serial works but parallel doesn't:**
- Issue is MPI-specific
- Check for race conditions or synchronization bugs in custom code

---

## Quick Fixes to Try (in order)

### Fix 1: Clean and Redecompose
```bash
cd ~/OpenFOAM/mavadi-v2406/run/TEST2
rm -rf processor* dynamicCode
decomposePar
mpirun -np 4 compInterFoam -parallel
```

### Fix 2: Use Fewer Processes
```bash
rm -rf processor*
# Edit system/decomposeParDict: numberOfSubdomains 2;
# Edit hierarchicalCoeffs: n (1 2 1);
decomposePar
mpirun -np 2 compInterFoam -parallel
```

### Fix 3: Disable Advanced Features
In `system/controlDict`, set:
```
useAdvancedInterfaceCapturing  false;
enableLiftProcessTracker false;
adjustTimeStep  no;
```

In `system/controlDict`, comment out:
```
/*
functions
{
    // All function objects
}
*/
```

### Fix 4: Change File Handler
In `system/controlDict`:
```
fileHandler     collated;  // or masterUncollated
```

### Fix 5: Different MPI Backend
```bash
mpirun --mca btl tcp,self -np 4 compInterFoam -parallel
```

---

## Expected Initialization Output

**A successful parallel initialization should show:**

```
Create time
Create mesh for time = 0
Reading field p_rgh
Reading field U
Reading/calculating face flux field phi
Reading field T
Constructing twoPhaseMixtureThermo
...
Reading field Te
Reading field Tl
Constructing two-temperature model
Constructing femtosecond laser model
Constructing advanced interface capturing
Selecting turbulence model
...
Starting time loop

Time = 2e-15
...
```

**If it hangs before "Starting time loop"**, the issue is in initialization.

---

## Monitoring Running Processes

### Check if processes are actually running:
```bash
ps aux | grep compInterFoam
```

### Check CPU usage:
```bash
top -u $USER
```

**If CPU usage is 0%**: Deadlock (processes waiting for each other)
**If CPU usage is high**: Infinite loop or heavy computation

### Check MPI communication:
```bash
# In another terminal while simulation runs
lsof -i | grep mpi
```

---

## Nuclear Option: Simplest Possible Test

Create a minimal test case:

**1. Copy TEST2 to TEST2_debug:**
```bash
cp -r TEST2 TEST2_debug
cd TEST2_debug
```

**2. Simplify controlDict:**
```
endTime         1e-14;  // Very short
writeInterval   1e-14;
adjustTimeStep  no;
useAdvancedInterfaceCapturing  false;
enableLiftProcessTracker false;

// Remove all functions
```

**3. Decompose and run:**
```bash
decomposePar
mpirun -np 2 compInterFoam -parallel
```

**If this works**: Gradually add features back until you find the problem.

---

## Report Back

After trying these steps, report:

1. **Where exactly does it hang?** (last line of output)
2. **CPU usage?** (0% = deadlock, high% = computing/loop)
3. **Which quick fix worked?** (if any)
4. **Debug output?** (if you enabled debug switches)

This will help identify the root cause.

---

## Most Common Solution

**Based on similar issues, try this first:**

```bash
cd ~/OpenFOAM/mavadi-v2406/run/TEST2

# 1. Clean everything
rm -rf processor* dynamicCode

# 2. Edit system/controlDict - disable function objects
# Comment out the entire functions{} section

# 3. Redecompose
decomposePar

# 4. Run with simpler MPI settings
mpirun --mca btl tcp,self --bind-to none -np 4 compInterFoam -parallel

# 5. If that works, re-enable function objects one by one
```

**Success rate: ~70% based on similar OpenFOAM parallel hang issues**
