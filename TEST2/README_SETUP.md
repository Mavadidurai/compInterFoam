# TEST2 Case Setup Instructions

## Issue: Missing Mesh Files

If you encounter this error when running `compInterFoam`:

```
--> FOAM FATAL ERROR: (openfoam-2406)
Cannot find file "points" in directory "TEST2/polyMesh" in times "0" down to constant
```

This means the mesh has not been generated yet.

## Solution

The TEST2 case requires mesh generation before running the solver. Follow these steps:

### Quick Start

1. **Generate mesh and run simulation:**
   ```bash
   cd TEST2
   ./Allrun
   ```

2. **Or, generate mesh only (without running solver):**
   ```bash
   cd TEST2
   ./Allrun.pre
   ```

3. **To clean and start fresh:**
   ```bash
   cd TEST2
   ./Allclean
   ```

### Manual Setup (Alternative)

If you prefer to run commands manually:

```bash
cd TEST2

# 1. Clean previous results (optional)
./Allclean

# 2. Generate the mesh
blockMesh

# 3. Check mesh quality
checkMesh

# 4. Copy initial conditions (if using non-zero startTime)
cp -r 0.orig 0

# 5. Set initial fields
setFields

# 6. Run the solver
compInterFoam
```

### Prerequisites

- OpenFOAM v2406 must be installed and sourced
- If OpenFOAM is not in your environment, source it first:
  ```bash
  source $HOME/OpenFOAM/OpenFOAM-v2406/etc/bashrc
  # or
  source /opt/openfoam/etc/bashrc
  ```

### Mesh Configuration

The mesh is configured in `system/blockMeshDict` with:
- **Total cells:** ~1.78 million (80 × 556 × 40)
- **Substrate:** 40 cells (8 μm thick)
- **Air gap:** 80 cells (12 μm thick) - critical for jet formation
- **Ti film:** 36 cells (71.4 nm thick) - ablation zone
- **Grading:** Optimized for laser absorption and jet dynamics

### Troubleshooting

1. **Command not found errors:** Make sure OpenFOAM is sourced
2. **blockMesh fails:** Check `log.blockMesh` for errors
3. **Mesh quality issues:** Review `log.checkMesh` for warnings

### Directory Structure After Setup

```
TEST2/
├── 0/                  # Initial fields (created from 0.orig)
├── 0.orig/            # Original initial conditions
├── constant/
│   ├── polyMesh/      # Generated mesh (created by blockMesh)
│   └── ...
├── system/
│   ├── blockMeshDict  # Mesh specification
│   └── ...
├── Allrun             # Full setup and run script
├── Allrun.pre         # Mesh generation only
└── Allclean           # Clean script
```

## Notes

- Mesh files (in `constant/polyMesh/`) are not committed to git as they are large binary files
- They must be regenerated each time you clone or clean the repository
- The `Allrun` and `Allrun.pre` scripts handle this automatically
