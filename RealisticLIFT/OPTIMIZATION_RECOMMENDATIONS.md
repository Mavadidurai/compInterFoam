# RealisticLIFT Simulation Optimization Guide

## Problem Summary

Your simulation has **critical performance issues**:

- ⏱️ **Speed**: 1360 seconds per time step → ~870 days to complete
- 🌡️ **Heating**: Lattice temperature only 339K (need >3560K for vaporization)
- 📉 **Time step**: Dropping to minimum 0.01 fs (too small)

---

## Root Causes

### 1. Extremely Small Time Step
```
Current minDeltaT: 1e-14 s (0.01 femtoseconds)
Adaptive stepping dropping to minimum every step
```

**Why**: Courant numbers are very low, but some stability constraint is forcing tiny steps.

### 2. Slow Lattice Heating
```
Electron temp (Te): 1412 K ✓ (rising)
Lattice temp (Tl): 339 K ✗ (barely rising)
Coupling: ~1200 W
```

**Why**: Electron-lattice coupling coefficient G = 1×10¹⁸ W/m³/K at low temps is weak.

### 3. Serial Execution
```
Mesh: ~1.5M cells
Cores used: 1
```

**Why**: No parallel decomposition configured.

---

## Solutions (In Order of Priority)

### 🚀 SOLUTION 1: Increase Time Step (CRITICAL)

**Change in `system/controlDict`:**

```cpp
// BEFORE (current):
deltaT          1e-13;         // Base time step (1 fs)
minDeltaT       1e-14;         // Minimum 0.01 fs for precision
maxDeltaT       2e-13;         // Cap adaptive growth (10 fs)
maxCo           0.1;           // Very conservative
maxAlphaCo      0.02;          // Very aggressive

// AFTER (optimized):
deltaT          5e-13;         // Base time step (0.5 fs) - 5x faster
minDeltaT       1e-13;         // Minimum 0.1 fs (10x larger)
maxDeltaT       1e-12;         // Cap at 1 fs (5x larger)
maxCo           0.5;           // Still stable, less conservative
maxAlphaCo      0.1;           // Relaxed for speed
```

**Expected improvement**: 10-50x faster

---

### ⚡ SOLUTION 2: Run in Parallel

**Steps:**

```bash
# 1. Create decomposeParDict
cd ~/OpenFOAM/mavadi-v2406/run/RealisticLIFT

cat > system/decomposeParDict << 'EOF'
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2406                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains 8;  // Use 8 cores

method          scotch;  // Automatic load balancing

// ************************************************************************* //
EOF

# 2. Decompose mesh
decomposePar

# 3. Run in parallel
mpirun -np 8 compInterFoam -parallel > log.compInterFoam 2>&1 &

# 4. Monitor
tail -f log.compInterFoam
```

**Expected improvement**: 6-8x faster (with 8 cores)

---

### 🔥 SOLUTION 3: Reduce Output Frequency

**Change in `system/controlDict`:**

```cpp
// BEFORE:
writeInterval   1e-14;         // Write every 0.01 ps (too frequent!)

// AFTER:
writeInterval   1e-12;         // Write every 1 ps (100x less frequent)
```

**Expected improvement**: 5-10x faster (less I/O overhead)

---

### 🌡️ SOLUTION 4: Increase Initial Electron-Lattice Coupling (Optional)

The coupling at low temps (G = 1×10¹⁸ W/m³/K) is weak. Consider:

**In `system/controlDict` → twoTemperatureProperties:**

```cpp
G
{
    type table;
    values
    (
        (300      5.0e18)   // INCREASED from 1.0e18 for faster initial heating
        (1000     3.0e18)
        (2000     6.0e18)
        (3000     1.0e19)
        (5000     2.0e19)
        (7000     3.0e19)
        (10000    5.0e19)
        (15000    7.0e19)
        (20000    1.0e20)
    );
}
```

This helps transfer energy from electrons to lattice faster during initial heating.

---

## 🎯 Recommended Quick Start Config

**Edit `system/controlDict` with these changes:**

```bash
cd ~/OpenFOAM/mavadi-v2406/run/RealisticLIFT

# Backup original
cp system/controlDict system/controlDict.original

# Apply optimizations (use your preferred editor)
nano system/controlDict
```

**Changes to make:**

1. Line 25: `deltaT 5e-13;` (was 1e-13)
2. Line 30: `minDeltaT 1e-13;` (was 1e-14)
3. Line 29: `maxDeltaT 1e-12;` (was 2e-13)
4. Line 27: `maxCo 0.5;` (was 0.1)
5. Line 28: `maxAlphaCo 0.1;` (was 0.02)
6. Line 38: `writeInterval 1e-12;` (was 1e-14)

---

## 📊 Expected Performance After Optimization

| Parameter | Before | After | Speedup |
|-----------|--------|-------|---------|
| Time step | 0.01 fs | 0.5-1 fs | 50-100x |
| Cores | 1 | 8 | 8x |
| Output freq | 0.01 ps | 1 ps | 100x |
| **Total speedup** | - | - | **~400-800x** |
| **Time to complete** | 870 days | **~1-2 days** | ✓ |

---

## 🔍 Diagnostic Commands

### Check current simulation status:
```bash
# If running
tail -f log.compInterFoam

# Extract progress
grep "^Time = " log.compInterFoam | tail -5

# Check temperatures
grep "max(Te):" log.compInterFoam | tail -5
grep "max(Tl):" log.compInterFoam | tail -5

# Check if phase change started
grep "Max |recoilPressure|" log.compInterFoam | tail -5
```

### Monitor Courant numbers:
```bash
grep "Courant Number" log.compInterFoam | tail -10
```

### Check execution time:
```bash
grep "ExecutionTime" log.compInterFoam | tail -1
```

---

## ⚠️ Important Notes

1. **Stability**: The optimized settings are still conservative. If you see divergence, reduce time step.

2. **Accuracy vs Speed**: We're trading some time resolution for practical runtime. The physics should still be captured correctly.

3. **Phase Change**: Once Tl > 3560K, recoil pressure will activate and simulation may slow down (expected).

4. **Monitoring**: Watch the first 100 ps carefully to ensure:
   - Te reaches >5000K
   - Tl reaches >3560K
   - Recoil pressure appears
   - No divergence

---

## 🔄 Testing Procedure

1. **Stop current simulation** (if running):
   ```bash
   # Find process
   ps aux | grep compInterFoam

   # Kill it
   killall compInterFoam
   ```

2. **Clean case**:
   ```bash
   ./Allclean
   ```

3. **Apply optimizations** (edit controlDict as above)

4. **Regenerate mesh** (if needed):
   ```bash
   blockMesh
   ```

5. **Run optimized simulation**:
   ```bash
   # Serial (for testing)
   compInterFoam > log.compInterFoam 2>&1 &

   # OR Parallel (recommended)
   decomposePar
   mpirun -np 8 compInterFoam -parallel > log.compInterFoam 2>&1 &
   ```

6. **Monitor progress**:
   ```bash
   tail -f log.compInterFoam

   # In another terminal, check temperatures every few minutes
   watch -n 60 'grep "max(Te):" log.compInterFoam | tail -1'
   ```

---

## 📈 Success Criteria

After 1-2 hours of runtime, you should see:

- ✓ Time > 100 ps (or 0.1 ns)
- ✓ Te > 5000 K
- ✓ Tl > 2000 K (approaching vaporization)
- ✓ Steady progress (not stuck)
- ✓ No "FOAM FATAL ERROR" messages

If recoil pressure appears (>0 MPa), congratulations! Your simulation is capturing the LIFT physics!

---

## 🆘 Troubleshooting

### If simulation diverges:
- Reduce deltaT by 2x
- Reduce maxCo to 0.3
- Increase nOuterCorrectors to 4

### If still too slow:
- Reduce writeInterval to 5e-12 or 1e-11
- Use more cores (16 or 32 if available)
- Consider coarser mesh

### If temperature not rising:
- Check laser is active: `grep "Laser input:" log.compInterFoam`
- Verify laser focus in metal film
- Increase G coefficient at low temps

---

**Good luck! The simulation should complete in 1-2 days with these optimizations.**
