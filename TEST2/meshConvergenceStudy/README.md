# Mesh Convergence Study for TEST2 (Graded Mesh)

## Overview

This directory contains a comprehensive mesh convergence study for **TEST2**, which uses a **graded mesh** configuration with refinement in critical regions for the compInterFoam LIFT simulation.

### TEST2 Configuration

- **Mesh type**: Graded (non-uniform with strategic refinement)
- **Simulation time**: 100 ps (1e-10 s)
- **Base mesh**: 80×400×40 (substrate), 80×400×80 (air gap), 80×400×36 (Ti film)
- **Total base cells**: ~4.99M cells
- **Grading strategy**:
  - Substrate: simpleGrading (1 1 0.5) - finer at top
  - Air gap: simpleGrading (1 1 2.0) - finer at bottom
  - Ti film: simpleGrading (1 1 0.67) - finer at top

## Directory Structure

```
meshConvergenceStudy/
├── README.md                 # This file
├── generateMeshes.py         # Generate blockMeshDict for each refinement level
├── setupCases.py             # Setup complete case directories
├── meshes/                   # Generated blockMeshDict files
│   ├── blockMeshDict.coarse
│   ├── blockMeshDict.medium
│   ├── blockMeshDict.fine
│   └── blockMeshDict.very_fine
├── cases/                    # Case directories (created by setupCases.py)
│   ├── coarse/
│   ├── medium/
│   ├── fine/
│   ├── very_fine/
│   └── runAll.sh             # Master script to run all cases
└── results/                  # Analysis outputs (to be created)
```

## Mesh Refinement Levels

The study uses 4 refinement levels with graded mesh refinement:

| Level       | Factor | Substrate      | Air Gap        | Ti Film        | Total Cells |
|-------------|--------|----------------|----------------|----------------|-------------|
| **Coarse**  | 0.5×   | 40×200×20      | 40×200×40      | 40×200×18      | ~0.62M      |
| **Medium**  | 1.0×   | 80×400×40      | 80×400×80      | 80×400×36      | ~4.99M      |
| **Fine**    | 1.5×   | 120×600×60     | 120×600×120    | 120×600×54     | ~16.85M     |
| **Very Fine**| 2.0× | 160×800×80     | 160×800×160    | 160×800×72     | ~39.94M     |

### Domain Geometry

- **X-direction**: 0 to 50 μm (laser spot centered at 25 μm)
- **Y-direction**: 0 to 20.0714 μm (layered structure)
  - Substrate (receiver): 0 to 8 μm (8 μm thick)
  - Air gap: 8 to 20 μm (12 μm gap)
  - **Titanium film**: 20 to 20.0714 μm (**71.4 nm thick** - critical region)
- **Z-direction**: 0 to 10 μm

### Grading Strategy

The graded mesh is designed to optimize resolution in critical regions:

1. **Substrate (grading 0.5)**:
   - Finer cells at top (interface with air gap)
   - Benefits: Better thermal coupling, smooth transition

2. **Air gap (grading 2.0)**:
   - Finer cells at bottom (near Ti film)
   - Benefits: Captures jet formation, resolves recoil pressure gradients

3. **Ti film (grading 0.67)**:
   - Finer cells at top (laser entry surface)
   - Benefits: Better resolution of laser absorption (Beer-Lambert law)
   - Peak absorption occurs at top surface where laser enters

### Physics-Based Refinement Rationale

- **Laser enters from TOP (donor side)** going DOWN
- **Beer-Lambert absorption**: I(y) = I₀ exp(-α(y₃-y))
- **Peak absorption** at top surface → finest cells needed there
- **Penetration depth**: ~17 nm (should be resolved by 4-6 cells)

## Quick Start

### 1. Meshes are already generated

The mesh files have been created in the `meshes/` directory with graded refinement.

### 2. Cases are already setup

The case directories have been created in the `cases/` directory.

### 3. Run Simulations

#### Option A: Run individual cases

```bash
cd cases/coarse
./Allrun
```

Monitor progress:
```bash
tail -f log.compInterFoam
```

#### Option B: Run all cases sequentially

```bash
cd cases
./runAll.sh
```

⚠️ **Warning**: This may take **hours to days** depending on your hardware!

**Computational requirements**:
- **Time step**: ~1 femtosecond (10⁻¹⁵ s)
- **Simulation time**: 100 ps
- **Steps**: ~100,000 steps
- **Recommended**: HPC cluster with parallel processing

### 4. Analyze Convergence

After all simulations complete, you can analyze the results using post-processing tools or custom scripts.

## Computational Cost Estimates

Approximate runtime for **100 ps simulation** on 8-core CPU:

| Mesh Level  | Total Cells | Runtime (est.) | Disk Space |
|-------------|-------------|----------------|------------|
| Coarse      | ~0.62M      | ~1-2 hours     | ~6 GB      |
| Medium      | ~4.99M      | ~10-15 hours   | ~30 GB     |
| Fine        | ~16.85M     | ~30-50 hours   | ~80 GB     |
| Very Fine   | ~39.94M     | ~80-140 hours  | ~180 GB    |

**Parallelization** can reduce wall-clock time significantly (near-linear scaling up to ~32-64 cores).

## Key Metrics to Analyze

The convergence analysis should examine:

1. **Temperature Fields**:
   - `T_peak`: Peak single-temperature field
   - `Te_peak`: Peak electron temperature (two-temperature model)
   - `Tl_peak`: Peak lattice temperature
   - Temperature distribution in Ti film (critical for laser absorption)

2. **Fluid Dynamics**:
   - `U_mag_max`: Maximum velocity magnitude
   - `p_max`: Maximum pressure
   - Jet velocity and trajectory

3. **Multiphase**:
   - `alpha.metal_max`: Maximum metal volume fraction
   - Interface evolution

4. **Mesh Quality**:
   - Total cells
   - Cell size distribution (due to grading)
   - Non-orthogonality
   - Skewness

## Grid Convergence Index (GCI) Methodology

For detailed GCI methodology and interpretation, see the main repository's meshConvergenceStudy/README.md.

Key steps:
1. Richardson Extrapolation
2. Observed Order of Accuracy
3. Grid Convergence Index calculation
4. Extrapolated value estimation
5. Convergence assessment

### Convergence Criteria

- **GCI < 1%**: Excellent convergence
- **GCI < 3%**: Good convergence
- **GCI < 5%**: Acceptable convergence
- **GCI > 5%**: Further refinement recommended

## Comparison with TEST1

**TEST2** uses a **graded mesh** while **TEST1** uses a **uniform mesh**:

| Aspect | TEST1 (Uniform) | TEST2 (Graded) |
|--------|-----------------|----------------|
| Substrate grading | (1 1 1) | (1 1 0.5) |
| Air gap cells | 60 | 80 |
| Air gap grading | (1 1 1) | (1 1 2.0) |
| Ti film grading | (1 1 1) | (1 1 0.67) |
| Base cells | ~4.35M | ~4.99M |
| Strategy | Uniform everywhere | Refined in critical regions |

### Expected Benefits of Grading (TEST2)

1. **Better laser absorption resolution**: Finer cells at top of Ti film
2. **Improved jet formation capture**: Finer cells at bottom of air gap
3. **Efficient use of cells**: Refinement where needed most
4. **Better thermal coupling**: Finer cells at substrate-air interface

### Comparison Study Goals

Compare TEST1 vs TEST2 to assess:
- **Solution accuracy**: Does grading improve key metrics?
- **Computational efficiency**: Better accuracy per cell?
- **Convergence behavior**: Different convergence rates?
- **Physics resolution**: Better capture of critical phenomena?

## Critical Regions and Resolution

### Ti Film (71.4 nm thick)

The titanium film is the most critical region:

- **Laser penetration depth**: ~9.71 nm (Ti @ 800 nm wavelength)
- **Required resolution**: 4-6 cells in penetration depth
- **Medium mesh**: 36 cells → ~2 nm average (graded finer at top)
- **Fine mesh**: 54 cells → ~1.3 nm average
- **Very fine mesh**: 72 cells → ~1.0 nm average

With grading (0.67), top cells are finer than bottom cells, providing better resolution where laser absorption is strongest.

### Air Gap (12 μm)

Critical for jet formation:

- **Expected jet width**: ~1-2 μm
- **Medium mesh**: 80 cells → ~150 nm average (finer at bottom)
- **Fine mesh**: 120 cells → ~100 nm average
- **Goal**: 10-15 cells across jet width

## Troubleshooting

### Common Issues

1. **"blockMesh failed"**
   - Check `log.blockMesh` for errors
   - Verify mesh definitions in `system/blockMeshDict`
   - Ensure geometry is valid (no negative volumes)
   - Check grading values are positive

2. **"Simulation diverges"**
   - Reduce time step (lower `maxCo` in controlDict)
   - Check initial conditions (0/)
   - Verify material properties (constant/)
   - Graded meshes may be more sensitive to time step

3. **"Out of memory"**
   - Use coarser mesh
   - Enable parallel decomposition
   - Increase system swap space

4. **"Poor mesh quality"**
   - Check checkMesh output
   - Verify grading values are reasonable (0.2 to 5.0 typically)
   - Consider adjusting grading factors

## References

1. **Roache, P.J.** (1994). "Perspective: A Method for Uniform Reporting of Grid Refinement Studies." *Journal of Fluids Engineering*, 116(3), 405-413.

2. **ASME V&V 20-2009**: "Standard for Verification and Validation in Computational Fluid Dynamics and Heat Transfer"

3. **Celik, I.B., et al.** (2008). "Procedure for Estimation and Reporting of Uncertainty Due to Discretization in CFD Applications." *Journal of Fluids Engineering*, 130(7), 078001.

4. **OpenFOAM blockMesh Documentation**: https://www.openfoam.com/documentation/user-guide/4-mesh-generation-and-conversion/4.3-mesh-generation-with-the-blockmesh-utility

## Notes

- Generated meshes use the same domain geometry as the original TEST2 case
- All cases use the same physics models and material properties
- Only the mesh resolution varies between cases
- Grading factors are preserved across all refinement levels
- Results should be compared to TEST1 (uniform mesh) to assess mesh grading impact

---

**Date**: 2025-11-12
**Mesh Type**: Graded (strategic refinement)
**Base Configuration**: TEST2
**Grading Philosophy**: Refine where physics demands it most
