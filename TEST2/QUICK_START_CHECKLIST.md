# TEST2 Quick Start Checklist

## Pre-Simulation Setup

### 1. Generate Mesh (REQUIRED)
```bash
cd /home/user/compInterFoam/TEST2
blockMesh
```
**Expected output:** ~4.99M cells created
**Time:** 5-10 minutes

### 2. Check Mesh Quality (REQUIRED)
```bash
checkMesh
```
**Look for:**
- Max non-orthogonality < 70 degrees
- Max skewness < 4
- No failed checks

### 3. Initialize Fields (REQUIRED)
```bash
setFields
```
**This sets:**
- Ti donor film (Y: 20.0-20.0714 μm)
- Air gap (Y: 8.0-20.0 μm)
- Substrate (Y: 0-8.0 μm)

### 4. Optional: Create Allrun Script
```bash
cat > Allrun << 'EOF'
#!/bin/bash
cd "${0%/*}" || exit
. $WM_PROJECT_DIR/bin/tools/RunFunctions

runApplication blockMesh
runApplication checkMesh
runApplication setFields
runApplication compInterFoam
EOF

chmod +x Allrun
```

## Run Simulation

### Option A: Serial Run (Simple)
```bash
compInterFoam
```
**Warning:** Will take 24-72 hours with 5M cells

### Option B: Parallel Run (Recommended)
```bash
# 1. Check decomposition settings
cat system/decomposeParDict

# 2. Decompose domain
decomposePar

# 3. Run in parallel (adjust -np based on your cores)
mpirun -np 8 compInterFoam -parallel

# 4. Reconstruct results
reconstructPar
```

## Monitoring During Run

### Check Progress
```bash
# View log file
tail -f log.compInterFoam

# Monitor residuals
foamMonitor -l postProcessing/residuals/0/residuals.dat
```

### Key Metrics to Watch
- **Time step:** Should vary between 1e-16 to 2e-15 s
- **Courant numbers:** Should stay < limits (0.2-0.3)
- **Temperature:** Max should reach 2000-6000 K
- **Velocity:** Max jet velocity 100-800 m/s
- **Pressure:** Recoil pressure ~70-80 MPa

## Quick Checks

### Verify Configuration
- [x] Config files present and valid
- [ ] Mesh generated (blockMesh)
- [ ] Mesh quality OK (checkMesh)
- [ ] Fields initialized (setFields)
- [ ] Sufficient disk space (~50 GB)

### Critical Parameters
- **Simulation time:** 0 to 1e-10 s (100 ps)
- **Base time step:** 2e-15 s (adaptive)
- **Write interval:** 1e-14 s (~1000 time directories)
- **Cell count:** 4,992,000 cells
- **Laser pulse:** 60 nJ, 200 fs, 343 nm

## Troubleshooting

### Mesh Generation Fails
```bash
# Check blockMeshDict syntax
blockMesh -dict system/blockMeshDict
```

### Simulation Crashes
1. Check log file: `log.compInterFoam`
2. Verify initial conditions: `paraFoam`
3. Reduce time step in controlDict
4. Check disk space: `df -h`

### Slow Performance
1. Use parallel execution (Option B above)
2. Reduce write frequency (increase writeInterval)
3. Enable purgeWrite to save disk space

## Post-Processing

### Visualization
```bash
# ParaView
paraFoam

# Or reconstruct first if parallel
reconstructPar
paraFoam
```

### Probe Data
```bash
# Jet velocity data
cat postProcessing/jetProbes/*/U

# Field extrema
cat postProcessing/fieldMinMax/*/fieldMinMax.dat
```

## Expected Results

**Physical Phenomena:**
1. Laser absorption in Ti film (0-200 fs)
2. Rapid heating to 2000-6000 K
3. Melting and vaporization
4. Recoil pressure buildup (70-80 MPa)
5. Jet formation and ejection (500-1000 fs)
6. Transfer to substrate

**Key Outputs:**
- Mid-plane VTK files
- Jet probe data (8 locations)
- Field min/max values
- Full 3D fields at 10 fs intervals

## Quick Reference

| Parameter | Value | Location |
|-----------|-------|----------|
| Solver | compInterFoam | system/controlDict:18 |
| End time | 1e-10 s | system/controlDict:23 |
| Time step | 2e-15 s | system/controlDict:25 |
| Write interval | 1e-14 s | system/controlDict:28 |
| Cells | 4,992,000 | system/blockMeshDict |
| Ti film | 71.4 nm | system/blockMeshDict:35 |
| Laser energy | 60 nJ | constant/laserProperties:53 |
| Pulse width | 200 fs | constant/laserProperties:54 |
| Wavelength | 343 nm | constant/laserProperties:55 |

---

**For detailed information, see:** PRE_SIMULATION_CHECK.md
