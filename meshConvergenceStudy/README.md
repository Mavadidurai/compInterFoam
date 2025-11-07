# Mesh Convergence Study for CompInterFoam LIFT Simulation

## Overview

This directory contains a comprehensive mesh convergence study framework for the **compInterFoam** LIFT (Laser-Induced Forward Transfer) simulation. The study implements the **Grid Convergence Index (GCI)** methodology following ASME V&V 20-2009 standards for verification and validation in CFD.

### Purpose

The mesh convergence study determines the mesh independence of the simulation results by:
- Generating multiple mesh resolutions (coarse, medium, fine, very fine)
- Running simulations on each mesh
- Comparing key solution metrics across meshes
- Computing Richardson extrapolation and GCI values
- Assessing convergence quality

## Directory Structure

```
meshConvergenceStudy/
├── README.md                    # This file
├── generateMeshes.py            # Generate blockMeshDict for each mesh level
├── setupCases.py                # Setup complete case directories
├── analyzeConvergence.py        # Post-processing and GCI analysis
├── meshes/                      # Generated blockMeshDict files
│   ├── blockMeshDict.coarse
│   ├── blockMeshDict.medium
│   ├── blockMeshDict.fine
│   └── blockMeshDict.very_fine
├── meshStudy/                   # Case directories (created by setupCases.py)
│   ├── coarse/
│   ├── medium/
│   ├── fine/
│   ├── very_fine/
│   └── runAll.sh                # Master script to run all cases
└── convergenceResults/          # Analysis outputs (created by analyzeConvergence.py)
    ├── convergence_plots.png    # Solution vs mesh size plots
    ├── gci_values.png           # GCI bar chart
    └── convergence_report.txt   # Detailed text report
```

## Mesh Refinement Levels

The study uses 4 refinement levels based on the original **LiftTest1** mesh configuration:

| Level       | Factor | Base Cells (Nx×Ny×Nz) | Total Cells | Ti Film Cells (z) | Cell Height (nm) |
|-------------|--------|----------------------|-------------|-------------------|------------------|
| **Coarse**  | 0.5×   | 40×200×20 (substrate)| ~326k       | 8                 | ~8.9             |
| **Medium**  | 1.0×   | 80×400×40 (substrate)| ~1.3M       | 16                | ~4.5             |
| **Fine**    | 1.5×   | 120×600×60           | ~2.9M       | 24                | ~3.0             |
| **Very Fine**| 2.0× | 160×800×80           | ~5.2M       | 32                | ~2.2             |

### Domain Geometry

- **X-direction**: 0 to 50 μm (laser spot centered at 25 μm)
- **Y-direction**: 0 to 32.0714 μm (layered structure)
  - Substrate (receiver): 0 to 8 μm (8 μm thick)
  - Air gap: 8 to 20 μm (12 μm gap)
  - **Titanium film**: 20 to 20.0714 μm (**71.4 nm thick** - critical region)
  - Donor glass: 20.0714 to 32.0714 μm (12 μm thick)
- **Z-direction**: 0 to 10 μm

### Refinement Strategy

- **Uniform refinement**: All regions scaled by the same factor
- **Minimum constraint**: Ti film always has ≥8 cells in thickness
- **Refinement ratio**:
  - Coarse→Medium: r = 2.0
  - Medium→Fine: r = 1.5
  - Fine→Very Fine: r ≈ 1.33

## Usage

### Step 1: Generate Mesh Definitions

Generate blockMeshDict files for all refinement levels:

```bash
cd meshConvergenceStudy
python3 generateMeshes.py
```

This creates mesh definitions in the `meshes/` directory.

**Custom levels** (optional):
```bash
python3 generateMeshes.py coarse medium fine  # Only generate 3 levels
```

### Step 2: Setup Case Directories

Create complete case directories for each mesh level:

```bash
python3 setupCases.py
```

This:
- Copies case structure from `LiftTest1`
- Installs appropriate blockMeshDict for each level
- Modifies controlDict for shorter runs (100 ps instead of 2 ns)
- Creates Allrun/Allclean scripts for each case

**Output**: Case directories in `meshStudy/`

### Step 3: Run Simulations

#### Option A: Run individual cases

```bash
cd meshStudy/coarse
./Allrun
```

Monitor progress:
```bash
tail -f log.compInterFoam
```

#### Option B: Run all cases sequentially

```bash
cd meshStudy
./runAll.sh
```

⚠️ **Warning**: This may take **hours to days** depending on your hardware!

**Computational requirements**:
- **Time step**: ~1 femtosecond (10⁻¹⁶ s)
- **Simulation time**: 100 ps (testing) or 2 ns (full)
- **Steps**: ~1-20 million steps
- **Recommended**: HPC cluster with parallel processing

#### Option C: Run in parallel (if MPI available)

Modify Allrun scripts to use `mpirun`:
```bash
# In each case's Allrun script, replace:
compInterFoam > log.compInterFoam 2>&1

# With:
mpirun -np 8 compInterFoam -parallel > log.compInterFoam 2>&1
```

Then decompose mesh before running:
```bash
decomposePar
```

### Step 4: Analyze Convergence

After all simulations complete:

```bash
python3 analyzeConvergence.py
```

**Outputs**:
- `convergenceResults/convergence_plots.png`: Key metrics vs mesh size
- `convergenceResults/gci_values.png`: GCI bar chart
- `convergenceResults/convergence_report.txt`: Detailed text report

## Key Metrics Analyzed

The convergence analysis examines:

1. **Temperature Fields**:
   - `T_peak`: Peak single-temperature field
   - `Te_peak`: Peak electron temperature (two-temperature model)
   - `Tl_peak`: Peak lattice temperature

2. **Fluid Dynamics**:
   - `U_mag_max`: Maximum velocity magnitude
   - `p_max`: Maximum pressure

3. **Multiphase**:
   - `alpha.metal_max`: Maximum metal volume fraction

4. **Mesh Quality**:
   - Total cells
   - Characteristic cell size (h)
   - Non-orthogonality
   - Skewness

## Grid Convergence Index (GCI) Methodology

The analysis uses the GCI method (Roache, 1994):

### 1. Richardson Extrapolation

For three consecutive meshes (coarse, medium, fine), compute:

```
ε₁₂ = φ₂ - φ₁  (coarse-medium difference)
ε₂₃ = φ₃ - φ₂  (medium-fine difference)
```

### 2. Observed Order of Accuracy

```
p = |ln(ε₁₂/ε₂₃)| / ln(r)
```

where `r` is the refinement ratio between meshes.

### 3. Grid Convergence Index

```
GCI = (F_s × |ε₂₃/φ₃|) / (r^p - 1)
```

where:
- `F_s` = 1.25 (safety factor for 3+ grids)
- `φ₃` = solution on fine mesh

### 4. Extrapolated Value

```
φ_ext = φ₃ + ε₂₃/(r^p - 1)
```

### 5. Convergence Assessment

- **GCI < 1%**: Excellent convergence
- **GCI < 3%**: Good convergence
- **GCI < 5%**: Acceptable convergence
- **GCI > 5%**: Further refinement recommended

### 6. Convergence Types

- **Monotonic**: ε₁₂ and ε₂₃ have same sign (expected for good convergence)
- **Oscillatory**: ε₁₂ and ε₂₃ have opposite signs (may indicate numerical issues)

## Interpreting Results

### Example Output

```
T_peak:
  Coarse:  4523.45 K
  Medium:  4587.32 K
  Fine:    4612.18 K
  Order of accuracy: 2.14
  GCI (fine): 0.87%
  Convergence: monotonic
  Extrapolated value: 4623.45 K
```

**Interpretation**:
- ✓ **Order ~2**: Near theoretical 2nd-order accuracy (good spatial schemes)
- ✓ **GCI < 1%**: Excellent mesh independence
- ✓ **Monotonic**: Smooth convergence
- ✓ **Small difference**: Fine mesh differs from extrapolated by <0.3%

### Red Flags

❌ **Order < 1**: Potential discretization issues or insufficient refinement
❌ **GCI > 5%**: Mesh-dependent results, refine further
❌ **Oscillatory**: May indicate numerical instability or under-resolution
❌ **Negative order**: Non-convergent behavior, check simulation setup

## Recommendations

### For Production Simulations

Based on convergence results:

1. **If GCI < 1% for all metrics**: Use **fine** mesh (good accuracy/cost balance)
2. **If GCI 1-3%**: Use **very fine** mesh for critical applications
3. **If GCI > 5%**: Run additional refinement level and repeat study

### Mesh Selection Guidelines

| Application | Recommended Mesh | Reason |
|-------------|-----------------|--------|
| **Preliminary studies** | Coarse | Fast turnaround for parameter exploration |
| **Standard simulations** | Medium | Acceptable accuracy for most purposes |
| **High-accuracy work** | Fine | Well-resolved physics, manageable cost |
| **Publication/validation** | Very fine | Demonstrable mesh independence |

### Critical Region (Ti Film)

The 71.4 nm titanium film requires special attention:
- **Minimum**: 8 cells in thickness (~9 nm/cell)
- **Recommended**: 16+ cells (~4.5 nm/cell or finer)
- **Ideal**: 24-32 cells (~2-3 nm/cell)

Laser penetration depth (~17 nm) should be resolved by **at least 4-6 cells**.

## Computational Cost Estimates

Approximate runtime for **100 ps simulation** on 8-core CPU:

| Mesh Level  | Total Cells | Runtime (est.) | Disk Space |
|-------------|-------------|----------------|------------|
| Coarse      | ~326k       | ~4-6 hours     | ~10 GB     |
| Medium      | ~1.3M       | ~12-18 hours   | ~25 GB     |
| Fine        | ~2.9M       | ~24-36 hours   | ~50 GB     |
| Very Fine   | ~5.2M       | ~48-72 hours   | ~80 GB     |

**Full 2 ns simulation**: Multiply times by ~20×

**Parallelization** can reduce wall-clock time significantly (near-linear scaling up to ~32-64 cores).

## Troubleshooting

### Common Issues

1. **"blockMesh failed"**
   - Check `log.blockMesh` for errors
   - Verify mesh definitions in `system/blockMeshDict`
   - Ensure geometry is valid (no negative volumes)

2. **"Simulation diverges"**
   - Reduce time step (lower `maxCo` in controlDict)
   - Check initial conditions (0.orig/)
   - Verify material properties (constant/)

3. **"Out of memory"**
   - Use coarser mesh
   - Enable parallel decomposition
   - Increase system swap space

4. **"Analysis script fails"**
   - Ensure simulations completed successfully
   - Check for time directories in case folders
   - Verify field files exist in latest time

5. **"No convergence detected"**
   - Run for longer simulation time
   - Ensure transients have settled
   - Check if steady-state is achievable

### Getting Help

- Check OpenFOAM documentation: https://www.openfoam.com/documentation
- OpenFOAM forum: https://www.cfd-online.com/Forums/openfoam/
- Review compInterFoam source code for physics details

## References

1. **Roache, P.J.** (1994). "Perspective: A Method for Uniform Reporting of Grid Refinement Studies." *Journal of Fluids Engineering*, 116(3), 405-413.

2. **ASME V&V 20-2009**: "Standard for Verification and Validation in Computational Fluid Dynamics and Heat Transfer"

3. **Celik, I.B., et al.** (2008). "Procedure for Estimation and Reporting of Uncertainty Due to Discretization in CFD Applications." *Journal of Fluids Engineering*, 130(7), 078001.

4. **OpenFOAM User Guide**: https://www.openfoam.com/documentation/user-guide

## License

This mesh convergence study framework is provided as-is for use with the compInterFoam solver. See main repository for license information.

## Author

Generated for compInterFoam LIFT simulation
Date: 2025-11-07

---

## Quick Start Summary

```bash
# 1. Generate meshes
python3 generateMeshes.py

# 2. Setup cases
python3 setupCases.py

# 3. Run simulations
cd meshStudy
./runAll.sh  # or run individual cases

# 4. Analyze results
cd ..
python3 analyzeConvergence.py

# 5. Review results
cat convergenceResults/convergence_report.txt
```

**Enjoy your mesh convergence study! 🚀**
