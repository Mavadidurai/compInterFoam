# Laser Energy Deposition Diagnostic Analysis

## Problem Summary

The femtosecond laser model is reporting:
1. "Temporal envelope expects X J this step, but spatial integration returned zero power"
2. "No cells found in laser beam path at time 0.018 ps!"

## Configuration Review

### Laser Parameters (from TEST2/constant/laserProperties)
- **Focus position**: (25e-6, 20.0357e-6, 5e-6) m = (25, 20.0357, 5) μm
- **Direction**: (0, -1, 0) - pointing downward in Y
- **Spot size**: 6.0e-6 m = 6 μm diameter (3 μm radius)
- **Film Y bounds**: [20.0, 20.0714] μm (71.4 nm thick)
- **Pulse energy**: 60 nJ
- **Pulse width**: 200 fs

### Mesh Parameters (from TEST2/system/blockMeshDict)
- **Domain**: X ∈ [0, 50] μm, Y ∈ [0, 20.0714] μm, Z ∈ [0, 10] μm
- **Resolution**:
  - X: 80 cells (625 nm/cell)
  - Y: 156 cells total (variable with grading)
  - Z: 40 cells (250 nm/cell)
- **Ti film**: 38 cells in Y direction for 71.4 nm thickness (~1.9 nm/cell)

## Theoretical Analysis

### Beam Geometry
The `isInBeam()` function checks if a cell center is within:
1. **Radial distance**: R ≤ 3 × beamRadius = 9 μm from beam axis
2. **Axial distance**: |z| ≤ axialHalfLength ≈ 3 μm along beam direction

For the configuration:
- Focus at Y = 20.0357 μm (middle of Ti film)
- Beam pointing down (-Y direction)
- Cells should be found if: 17.0357 μm ≤ Y ≤ 23.0357 μm
- Ti film (20.0 to 20.0714 μm) is well within this range ✓

### Expected Behavior
With the focus in the middle of the domain and film, hundreds of cells should be found within the beam path.

## Diagnostic Enhancements Added

### 1. Search Box Diagnostics (Line ~1320)
Added output showing:
- Focus position and beam direction
- Spot size and beam radius
- Search box bounds
- Mesh bounds
- Total cell count

### 2. Cell Processing Statistics (Line ~1648)
Added tracking of:
- Total mesh cells
- Cells in search box
- Cells passing isInBeam() test

### 3. Enhanced Error Messages (Line ~2215)
When no cells are found, now displays:
- Focus position
- Beam direction
- Spot size
- Film Y bounds
- Mesh bounds
- Total cell count

## Potential Root Causes

### 1. Mesh Not Generated (MOST LIKELY)
- No `polyMesh` directory found in constant/
- No time directories (0, 0.001, etc.)
- User needs to run `blockMesh` first

### 2. Coordinate System Mismatch
- blockMeshDict uses `scales 1e-6` (coordinates in μm)
- laserProperties uses SI units (m)
- Should be compatible, but worth verifying with diagnostics

### 3. Numerical Precision Issues
- Very thin film (71.4 nm) with coarse mesh in X/Z
- Cell centers might not align well with beam axis

## Next Steps

### To Run Diagnostics:
```bash
cd TEST2
blockMesh                    # Generate mesh
cp -r 0.orig 0              # Copy initial conditions
<solver_name>               # Run solver (will show diagnostic output)
```

### What to Look For:
1. Check if "Total mesh cells" shows > 0
2. Check if "Cells in search box" shows > 0
3. Compare search box bounds with mesh bounds
4. Verify focus position is inside mesh

### If Mesh is Generated and Issue Persists:
- Check units in mesh vs. laser properties
- Verify film Y bounds match blockMeshDict
- Consider increasing spot size temporarily to test
- Check if scanVelocity is moving focus outside domain

## Code Changes Summary

**File**: `femtosecondLaserModel.C`
- Added diagnostic output in `applySpatialWeighting()` function
- Enhanced warning message in `calculateSource()` function
- Added cell counting variables to track search box efficiency
- All changes are non-intrusive (diagnostic output only)
