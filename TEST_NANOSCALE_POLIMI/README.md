# Ultra-Fine Nanoscale LIFT Simulation for PoliMi HPC

## Overview

This case contains a complete ultra-fine nanoscale simulation of Laser-Induced Forward Transfer (LIFT) optimized for high-performance computing on the PoliMi cluster.

**Key Features:**
- **Mesh Resolution**: 1-5 nm cells in critical regions (nanoscale!)
- **Total Cells**: ~102 Million cells
- **Physics**: Complete LIFT cycle from heating to solidification
- **Time Range**: 0 → 2 nanoseconds
- **Parallelization**: 128-256 cores recommended

---

## Simulation Details

### Geometry
```
Domain (bottom to top):
├── Substrate (receiver):  8.0 μm thick
├── Air gap:              12.0 μm thick
└── Ti donor film:        71.4 nm thick
```

**Dimensions**: 50 μm × 20.0714 μm × 10 μm (x × y × z)

### Mesh Specifications

| Region       | Dimensions (μm)      | Cells (x×y×z)    | Resolution         | Total Cells |
|--------------|---------------------|------------------|--------------------|-------------|
| Substrate    | 50 × 8.0 × 10       | 250×1600×100     | 200nm × 5nm × 100nm| 40 M        |
| Air gap      | 50 × 12.0 × 10      | 250×2400×100     | 200nm × 5nm × 100nm| 60 M        |
| Ti film      | 50 × 0.0714 × 10    | 250×72×100       | 200nm × 1nm × 100nm| 1.8 M       |
| **TOTAL**    |                     |                  |                    | **102 M**   |

**Critical Achievement**: **1 nm resolution** in Ti film (Y-direction)

### Physical Processes

1. **Laser Heating (0-200 ps)**
   - Femtosecond laser pulse: 200 fs FWHM, 60 nJ, λ=343 nm
   - Two-temperature model (electron-phonon coupling)
   - Gaussian spatial profile: 6 μm spot size
   - Peak absorption in Ti film

2. **Phase Transitions**
   - Melting: T > 1941 K
   - Superheated evaporation: T > 2200 K
   - Clausius-Clapeyron kinetics
   - Recoil pressure up to 3 GPa

3. **Material Transfer (200-500 ps)**
   - High-rate evaporation
   - Recoil-driven ejection
   - Interface tracking with advanced capturing

4. **Cooling & Solidification (500-2000 ps)**
   - Thermal diffusion
   - Re-solidification
   - Final structure formation

### Numerical Settings

**Time Stepping:**
- Initial Δt: 0.1 fs (1×10⁻¹⁶ s)
- Adaptive: yes
- Maximum Δt: 5 fs (5×10⁻¹⁵ s)
- Courant limits: Co=0.15, αCo=0.25, Thermal=0.15

**Solvers:**
- Pressure: GAMG with DIC preconditioner
- Velocity: smoothSolver with GaussSeidel
- Temperature: GAMG with DIC
- High-precision: 10 decimal places

---

## Quick Start Guide

### Prerequisites

1. **OpenFOAM v2406** or compatible version
2. **PoliMi HPC cluster** access with:
   - 128-256 CPU cores
   - 512-1024 GB RAM
   - High-speed parallel filesystem (Lustre/GPFS)
3. **Compiled compInterFoam solver** with custom models

### Setup Instructions

#### 1. Transfer to PoliMi HPC

```bash
# On your local machine
scp -r TEST_NANOSCALE_POLIMI your_username@polimi-hpc.it:~/compInterFoam/

# SSH to PoliMi HPC
ssh your_username@polimi-hpc.it
cd ~/compInterFoam/TEST_NANOSCALE_POLIMI
```

#### 2. Modify Paths

Edit the following files for your HPC environment:

**Allrun:**
```bash
# Line ~10: Update OpenFOAM path
source /path/to/your/openfoam/etc/bashrc
```

**run_nanoscale_lift.slurm:**
```bash
# Lines 6-9: Update SLURM settings
#SBATCH --account=your_account_name
#SBATCH --partition=your_partition

# Line 41: Update OpenFOAM path
FOAM_PATH="/path/to/openfoam/etc/bashrc"
```

#### 3. Generate Mesh & Decompose

```bash
# Make scripts executable
chmod +x Allrun Allclean

# Generate mesh (can run on login node, takes ~5-10 min)
blockMesh | tee log.blockMesh

# Check mesh quality
checkMesh -allTopology -allGeometry | tee log.checkMesh

# Prepare initial conditions
cp -r 0.orig 0

# Decompose for parallel run (takes ~10-20 min)
decomposePar -copyZero | tee log.decomposePar
```

#### 4. Submit to SLURM

```bash
# Submit job
sbatch run_nanoscale_lift.slurm

# Check job status
squeue -u $USER

# Monitor output
tail -f slurm_JOBID.out

# Monitor log
tail -f log.compInterFoam
```

---

## Computational Requirements

### Resource Estimates

| Configuration | Cores | Nodes | Memory  | Cells/Core | Est. Time/Step |
|--------------|-------|-------|---------|------------|----------------|
| Small        | 64    | 2     | 256 GB  | 1.6M       | ~1.0 s         |
| **Recommended** | **128** | **4** | **512 GB** | **800k** | **~0.3 s** |
| Large        | 256   | 8     | 1024 GB | 400k       | ~0.15 s        |
| X-Large      | 512   | 16    | 2048 GB | 200k       | ~0.08 s        |

### Performance Projections

**For 128 cores (recommended):**
- Timesteps needed: ~1,000,000 (for 2 ns with avg Δt ≈ 2 fs)
- Time per step: ~0.3 seconds
- Total compute time: ~300,000 seconds = **83 hours**
- Recommended walltime: **72-96 hours**

**Storage Requirements:**
- Mesh: ~50 GB (decomposed)
- Full time history (binary): ~500 GB - 2 TB
- Reduced output: ~100-200 GB

**Recommended Settings:**
```
writeInterval   2e-13;  // Every 200 fs → ~10,000 output times
purgeWrite      200;    // Keep last 200 times → ~40 time directories
```

---

## Running the Simulation

### Method 1: Using Allrun + SLURM (Recommended)

```bash
# Edit SLURM script first
nano run_nanoscale_lift.slurm

# Submit
sbatch run_nanoscale_lift.slurm
```

### Method 2: Step-by-step

```bash
# 1. Mesh generation
./Allrun  # Generates mesh and decomposes

# 2. Submit to SLURM
sbatch run_nanoscale_lift.slurm

# 3. Monitor progress
tail -f log.compInterFoam
squeue -u $USER
```

### Method 3: Interactive (Testing Only)

```bash
# Request interactive node
srun --nodes=1 --ntasks-per-node=32 --time=1:00:00 --pty bash

# Load environment
source /path/to/openfoam/etc/bashrc

# Run short test
mpirun -np 32 compInterFoam -parallel | tee log.test
```

---

## Monitoring & Diagnostics

### Check Simulation Progress

```bash
# Monitor main log
tail -f log.compInterFoam

# Extract current time
grep "^Time = " log.compInterFoam | tail -5

# Count timesteps
grep "^Time = " log.compInterFoam | wc -l

# Check for errors/warnings
grep -i "error\|warning" log.compInterFoam

# Monitor SLURM output
tail -f slurm_*.out
```

### Key Metrics to Watch

```bash
# Extract Courant numbers
grep "Courant Number" log.compInterFoam | tail -20

# Check temperature extremes
grep "max(T)" log.compInterFoam | tail -20

# Check mass conservation
grep "alpha.metal" log.compInterFoam | grep "volume" | tail -20
```

### Post-Processing

```bash
# Reconstruct parallel data (WARNING: Large files!)
reconstructPar

# Reconstruct specific times
reconstructPar -time '1e-10,2e-10,5e-10,1e-9,2e-9'

# Reconstruct only latest time
reconstructPar -latestTime

# Extract mid-plane data
reconstructPar -fields '(alpha.metal T Tl Te U p_rgh)'

# Load in ParaView
paraFoam
```

---

## Expected Results

### Key Phenomena to Observe

1. **Phase 1: Laser Absorption (0-200 ps)**
   - Rapid temperature rise in Ti film
   - Electron temperature peaks at ~10,000-20,000 K
   - Lattice temperature reaches ~3,000-6,000 K
   - Melting initiates at ~1941 K

2. **Phase 2: Material Response (200-500 ps)**
   - Superheating above 2200 K
   - Recoil pressure builds (up to 3 GPa)
   - Material ejection begins
   - Jet/droplet formation

3. **Phase 3: Transfer (500-1000 ps)**
   - Material crosses air gap
   - Cooling during flight
   - Interface deformation

4. **Phase 4: Solidification (1000-2000 ps)**
   - Deposition on substrate
   - Thermal equilibration
   - Final structure formation

### Output Files

- **Time directories**: Field data (T, Tl, Te, p, U, alpha.metal)
- **postProcessing/**: Probe data, surface data
- **log.compInterFoam**: Solver log with diagnostics
- **VTK/**: Mid-plane visualization (if enabled)

---

## Troubleshooting

### Common Issues

**1. Simulation crashes early**
```
Possible causes:
- Initial conditions not set properly
- Timestep too large
- Check: log.compInterFoam for errors

Solution:
- Reduce maxDeltaT in system/controlDict
- Check initial T, Tl, Te fields in 0/
```

**2. Very slow timestep**
```
Possible causes:
- Courant number limit too strict
- High thermal gradients

Solution:
- Monitor Courant numbers in log
- Check if maxThermalCourant can be increased to 0.2-0.3
```

**3. Temperature instabilities**
```
Possible causes:
- Electron substepping insufficient
- Coupling iterations too few

Solution:
- Increase maxElectronSubCycles to 100 in controlDict
- Increase nInnerCouplingSweeps to 5
```

**4. Memory issues**
```
Possible causes:
- Too many time directories
- Insufficient node memory

Solution:
- Reduce purgeWrite to 50
- Increase writeInterval
- Request more memory in SLURM script
```

**5. I/O bottleneck**
```
Possible causes:
- Too frequent writes
- Slow filesystem

Solution:
- Increase writeInterval to 5e-13
- Use binary format (already set)
- Check filesystem performance
```

---

## Customization

### Adjust Core Count

Edit `system/decomposeParDict`:

```cpp
// For 256 cores:
numberOfSubdomains 256;
hierarchicalCoeffs
{
    n    (16 8 2);  // 16×8×2 = 256
}
```

Then re-decompose:
```bash
reconstructPar
rm -rf processor*
decomposePar -copyZero
```

### Extend Simulation Time

Edit `system/controlDict`:
```cpp
endTime    5e-9;  // Extend to 5 ns
```

### Modify Output Frequency

Edit `system/controlDict`:
```cpp
writeInterval   5e-13;  // Write every 500 fs
purgeWrite      100;    // Keep last 100 times
```

---

## File Structure

```
TEST_NANOSCALE_POLIMI/
├── 0.orig/                      # Initial conditions templates
│   ├── T, Tl, Te               # Temperature fields
│   ├── p, p_rgh                # Pressure fields
│   ├── U                        # Velocity field
│   └── alpha.metal, alpha.air   # Phase fractions
├── constant/                    # Material properties
│   ├── laserProperties          # Laser parameters
│   ├── thermophysicalProperties.*  # Material properties
│   ├── transportProperties      # Transport coefficients
│   └── g                        # Gravity
├── system/                      # Numerical settings
│   ├── controlDict              # Time control & physics
│   ├── blockMeshDict            # Mesh generation
│   ├── decomposeParDict         # Parallel decomposition
│   ├── fvSchemes                # Discretization schemes
│   └── fvSolution               # Linear solvers
├── Allrun                       # Setup script
├── Allclean                     # Clean script
├── run_nanoscale_lift.slurm     # SLURM submission script
└── README.md                    # This file
```

---

## Technical Support

### Getting Help

1. **Check logs**: Start with `log.compInterFoam` and `log.blockMesh`
2. **Mesh quality**: Run `checkMesh -allTopology -allGeometry`
3. **Decomposition**: Verify with `ls -d processor* | wc -l`
4. **SLURM**: Check `squeue`, `sacct`, and SLURM output files

### Key Parameters Summary

| Parameter | File | Value | Purpose |
|-----------|------|-------|---------|
| Cells | blockMeshDict | 102M | Mesh resolution |
| Cores | decomposeParDict | 128 | Parallelization |
| endTime | controlDict | 2e-9 | Simulation duration |
| deltaT | controlDict | 1e-16 | Initial timestep |
| maxDeltaT | controlDict | 5e-15 | Maximum timestep |
| writeInterval | controlDict | 2e-13 | Output frequency |
| pulseEnergy | laserProperties | 6e-8 | Laser energy (60 nJ) |
| spotSize | laserProperties | 6e-6 | Beam diameter (6 μm) |

---

## Performance Tips

1. **Optimize I/O**: Use fast filesystem, minimize write frequency
2. **Load balancing**: Verify with `decomposePar -cellDist`
3. **Memory**: Monitor with `sacct -j JOBID --format=MaxRSS`
4. **Scaling**: Test strong/weak scaling before production runs
5. **Checkpointing**: Consider saving restart points every few hours

---

## Citation

If you use this simulation setup, please cite:
- Your relevant publications
- OpenFOAM: https://www.openfoam.com
- compInterFoam solver development

---

## Version History

- **v1.0** (2024-11-14): Initial nanoscale setup for PoliMi HPC
  - 102M cell mesh
  - Complete LIFT cycle (0-2 ns)
  - Optimized for 128-256 cores

---

## Contact

For questions about this setup:
- Check OpenFOAM documentation: https://www.openfoam.com/documentation
- OpenFOAM forums: https://www.cfd-online.com/Forums/openfoam/

---

**Good luck with your nanoscale LIFT simulations!**
