# Mesh Convergence Study for TEST1 (Uniform Mesh)

## Overview

This directory contains a comprehensive mesh convergence study for **TEST1**, which uses a **uniform mesh** configuration (no grading) for the compInterFoam LIFT simulation.

### TEST1 Configuration

- **Mesh type**: Uniform (simpleGrading 1 1 1 throughout)
- **Simulation time**: 100 ps (1e-10 s)
- **Base mesh**: 80×400×40 (substrate), 80×400×60 (air gap), 80×400×36 (Ti film)
- **Total base cells**: ~4.35M cells

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

The study uses 4 refinement levels with uniform refinement (no grading):

| Level       | Factor | Substrate      | Air Gap        | Ti Film        | Total Cells |
|-------------|--------|----------------|----------------|----------------|-------------|
| **Coarse**  | 0.5×   | 40×200×20      | 40×200×30      | 40×200×18      | ~0.54M      |
| **Medium**  | 1.0×   | 80×400×40      | 80×400×60      | 80×400×36      | ~4.35M      |
| **Fine**    | 1.5×   | 120×600×60     | 120×600×90     | 120×600×54     | ~14.69M     |
| **Very Fine**| 2.0× | 160×800×80     | 160×800×120    | 160×800×72     | ~34.82M     |

### Domain Geometry

- **X-direction**: 0 to 50 μm (laser spot centered at 25 μm)
- **Y-direction**: 0 to 20.0714 μm (layered structure)
  - Substrate (receiver): 0 to 8 μm (8 μm thick)
  - Air gap: 8 to 20 μm (12 μm gap)
  - **Titanium film**: 20 to 20.0714 μm (**71.4 nm thick** - critical region)
- **Z-direction**: 0 to 10 μm

### Refinement Strategy

- **Uniform refinement**: All cells scaled by the same factor
- **No grading**: simpleGrading (1 1 1) throughout
- **Refinement ratio**:
  - Coarse→Medium: r = 2.0
  - Medium→Fine: r = 1.5
  - Fine→Very Fine: r ≈ 1.33

## Quick Start

### 1. Meshes are already generated

The mesh files have been created in the `meshes/` directory.

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
| Coarse      | ~0.54M      | ~1-2 hours     | ~5 GB      |
| Medium      | ~4.35M      | ~8-12 hours    | ~25 GB     |
| Fine        | ~14.69M     | ~24-48 hours   | ~70 GB     |
| Very Fine   | ~34.82M     | ~72-120 hours  | ~150 GB    |

**Parallelization** can reduce wall-clock time significantly (near-linear scaling up to ~32-64 cores).

## Key Metrics to Analyze

The convergence analysis should examine:

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

## Comparison with TEST2

**TEST1** uses a **uniform mesh** while **TEST2** uses a **graded mesh** with refinement in critical regions:
- TEST1: All regions use simpleGrading (1 1 1)
- TEST2: Graded mesh with finer cells in critical regions

Compare the results between TEST1 and TEST2 to assess the impact of mesh grading on:
- Solution accuracy
- Computational efficiency
- Convergence behavior

## Troubleshooting

### Common Issues

1. **"blockMesh failed"**
   - Check `log.blockMesh` for errors
   - Verify mesh definitions in `system/blockMeshDict`
   - Ensure geometry is valid (no negative volumes)

2. **"Simulation diverges"**
   - Reduce time step (lower `maxCo` in controlDict)
   - Check initial conditions (0/)
   - Verify material properties (constant/)

3. **"Out of memory"**
   - Use coarser mesh
   - Enable parallel decomposition
   - Increase system swap space

4. **"No convergence detected"**
   - Run for longer simulation time
   - Ensure transients have settled
   - Check if steady-state is achievable

## References

1. **Roache, P.J.** (1994). "Perspective: A Method for Uniform Reporting of Grid Refinement Studies." *Journal of Fluids Engineering*, 116(3), 405-413.

2. **ASME V&V 20-2009**: "Standard for Verification and Validation in Computational Fluid Dynamics and Heat Transfer"

3. **Celik, I.B., et al.** (2008). "Procedure for Estimation and Reporting of Uncertainty Due to Discretization in CFD Applications." *Journal of Fluids Engineering*, 130(7), 078001.

## Notes

- Generated meshes use the same domain geometry as the original TEST1 case
- All cases use the same physics models and material properties
- Only the mesh resolution varies between cases
- Results should be compared to TEST2 (graded mesh) to assess mesh grading impact

---

**Date**: 2025-11-12
**Mesh Type**: Uniform (no grading)
**Base Configuration**: TEST1
