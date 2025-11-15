# LIFT Simulation Case - TEST2

This is an OpenFOAM case for simulating Laser-Induced Forward Transfer (LIFT) using the compInterFoam solver.

## Case Description

**Physics:**
- Two-phase compressible flow (metal + air)
- Laser heating and material transfer
- Metal phase: ρ=4515 kg/m³, ν=5.2×10⁻⁷ m²/s
- Air phase: ρ=1.2 kg/m³, ν=1.48×10⁻⁵ m²/s

**Geometry:**
- Domain: 50 × ~24 × 10 μm (X × Y × Z)
- Substrate (receiver): Y ∈ [0, 8] μm
- Air gap: Y ∈ [8, 20] μm
- Donor film: Y ∈ [20, 20.0714] μm (thickness ~71.4 nm)

**Features:**
- Lift process tracking enabled
- Mid-plane visualization (VTK output)
- Optional mesh refinement in laser focus region
- Parallel execution support

## Running the Case

### Option 1: Complete Run (Preprocessing + Solver)

```bash
./Allrun
```

This script will:
1. Generate mesh (`blockMesh`)
2. Create cell zones (`topoSet`)
3. Optionally refine mesh (`refineMesh`)
4. Initialize fields (`setFields`)
5. Decompose for parallel (`decomposePar`)
6. Run solver in parallel (`compInterFoam`)
7. Reconstruct results (`reconstructPar`)

### Option 2: Preprocessing Only

```bash
./Allrun.pre
```

Then manually run the solver:
```bash
mpirun -np 4 compInterFoam -parallel | tee log.lift
```

### Option 3: Manual Steps

```bash
# Clean previous results
./Allclean

# Generate mesh
blockMesh

# Create cell zones (required for fvOptions)
topoSet

# Initialize fields (alpha.metal, temperatures)
setFields

# Decompose for parallel run
decomposePar -copyZero

# Run solver
mpirun -np 4 compInterFoam -parallel | tee log.lift

# Reconstruct
reconstructPar
```

## Important Notes

**Before Running:**
- The solver **requires** proper preprocessing:
  - ✓ Mesh must be generated (`blockMesh`)
  - ✓ Cell zones must be created (`topoSet`)
  - ✓ Fields must be initialized (`setFields`)
  - ✓ Case must be decomposed (`decomposePar`) for parallel runs

**Common Errors:**
- `alpha.metal is effectively empty` → Run `setFields` first
- `No cells selected` for substrate → Run `topoSet` first
- Parallel run fails → Run `decomposePar` first

## Visualization

```bash
# Serial case
paraFoam

# Parallel case
paraFoam -builtin

# Or reconstruct first, then view
reconstructPar
paraFoam
```

## Monitoring Progress

Check solver log:
```bash
tail -f log.lift
```

Monitor specific metrics:
```bash
grep -i "courant" log.lift
grep -i "PIMPLE" log.lift
```

## Case Structure

```
TEST2/
├── 0.orig/          # Initial and boundary conditions
├── constant/        # Physical properties, mesh, fvOptions
├── system/          # Numerical schemes, solver settings, mesh generation
├── Allrun           # Complete run script
├── Allrun.pre       # Preprocessing only
└── Allclean         # Clean case script
```
