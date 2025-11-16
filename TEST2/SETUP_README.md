# OpenFOAM LIFT Simulation Setup Guide

## Issues Identified and Fixed

### 1. Missing Mesh
**Problem**: No mesh was generated before running the simulation.
**Solution**: Run `blockMesh` to generate the computational mesh.

### 2. Uninitialized Fields
**Problem**: The `alpha.metal` field was empty (max = 0), preventing proper simulation.
**Solution**:
- Copy `0.orig` to `0`
- Run `setFields` to initialize phase fractions according to geometry

### 3. Cell Zone Warnings
**Problem**: fvOptions referenced a "substrate" cellZone that doesn't exist.
**Solution**: Commented out unused constraints in `constant/fvOptions`.

### 4. Laser Beam Not Found (Decomposition Issue)
**Problem**: "No cells found in laser beam path" - laser focus at x=25μm was on processor boundary.
**Solution**: Changed decomposition from (2 2 1) to (1 4 1):
- **Old**: Split in X at 25μm → laser on boundary!
- **New**: Split only in Y-direction → laser fully contained in one X-slice
- Focus position: (25e-6, 20.0357e-6, 5e-6) - centered in Ti film ✓
- Film bounds: y ∈ [20.0000e-6, 20.0714e-6] m (71.4 nm thick)
- Laser penetration depth: ~9.7 nm (absorption coefficient: 1.03e8 m⁻¹)

## Geometry Summary

```
Domain: 50 μm × 20.0714 μm × 10 μm (X × Y × Z)

Layers (Y-direction, gravity in -Y):
├─ Substrate (Receiver):  0.0000 - 8.0000 μm  (8 μm,    Ti metal)
├─ Air Gap:              8.0000 - 20.0000 μm  (12 μm,   air)
└─ Donor Film (Ti):     20.0000 - 20.0714 μm  (71.4 nm, Ti metal)

Laser:
- Focus: (25, 20.0357, 5) μm (center of domain, middle of donor film)
- Direction: (0, -1, 0) - shooting downward into film
- Spot size: 6 μm diameter (3 μm 1/e² radius)
- Pulse energy: 60 nJ (0.2 J/cm² fluence threshold)
- Pulse width: 200 fs FWHM
- Wavelength: 343 nm (Ti absorption peak)
```

## Mesh Details

- **Total cells**: ~4.99M cells (80 × 400 × [40+80+38])
- **Grading**: Optimized for laser absorption and jet formation
  - Substrate: Finer at top (interface with air)
  - Air gap: Finer at bottom (jet formation zone) - 80 cells
  - Ti film: Finer at top (laser entry surface) - 38 cells
  - Finest cells in Ti film: ~1.5 nm at laser entry surface

## Setup Instructions

### Method 1: Using Allrun Script (Recommended)

```bash
cd ~/OpenFOAM/mavadi-v2406/run/TEST2
chmod +x Allrun
./Allrun
```

This will automatically:
1. Clean previous results
2. Copy initial fields
3. Generate mesh with blockMesh
4. Check mesh quality
5. Initialize fields with setFields
6. Decompose for parallel run (4 processors)

### Method 2: Manual Setup

```bash
cd ~/OpenFOAM/mavadi-v2406/run/TEST2

# Clean
./Allclean

# Copy initial fields
cp -r 0.orig 0

# Generate mesh
blockMesh

# Optional: Check mesh quality
checkMesh

# Initialize phase fields
setFields

# Decompose for parallel
decomposePar -force
```

## Running the Simulation

After setup is complete:

```bash
# Parallel run (recommended, 4 cores)
mpirun -np 4 compInterFoam -parallel | tee log.lift

# OR serial run (slower)
compInterFoam | tee log.lift
```

## Expected Behavior

After proper setup, you should see:
- ✓ Mesh cells properly initialized with metal and air phases
- ✓ Laser energy deposited in the donor film
- ✓ Phase change (melting/vaporization) in the Ti film
- ✓ Jet formation and material transfer

The simulation will run from t=0 to t=200 ps with adaptive time stepping.

## Verification Checklist

Before running, verify:
- [ ] `constant/polyMesh/` directory exists (mesh generated)
- [ ] `0/` directory exists with all field files
- [ ] `max(alpha.metal)` > 0 in initialized fields (check with `foamListTimes`)
- [ ] `processor0/`, `processor1/`, etc. exist (for parallel runs)

## Post-Processing

Results are saved in:
- `postProcessing/` - Field min/max values, probe data
- `VTK/` or time directories - Full field data
- `log.lift` - Simulation log with convergence info

## Troubleshooting

### "alpha.metal is effectively empty"
→ Run `setFields` to initialize phase fractions

### "No cells found in laser beam path"
**Most common cause**: Laser on processor boundary
→ Re-decompose with corrected decomposeParDict:
```bash
rm -rf processor*
decomposePar -force
```
The new decomposition (1 4 1) avoids splitting at the laser focus

### Parallel run fails immediately
→ Run `decomposePar` before parallel execution

### Already ran with old decomposition?
If you already decomposed and ran with the old (2 2 1) settings:
```bash
./Allclean          # Clean everything
./Allrun            # Re-setup with correct decomposition
```

### Very slow convergence
→ Check `log.lift` for PIMPLE iteration counts and residuals
→ May need to adjust time step or solver tolerances in `system/controlDict`
