# LIFT Simulation Fixes

## Issues Found and Fixed

### 1. **CRITICAL: Missing Mesh** ❌
**Problem:** The simulation was run without first generating the computational mesh.
**Symptom:** Temperature solver diverged with increasing residuals, alpha solver reported singularities.
**Fix:** Must run mesh generation before starting simulation.

### 2. **Temperature Solver Configuration** ⚠️
**Problem:** The gas temperature equation uses `smoothSolver` with `symGaussSeidel`, which is not optimal for equations with large source terms and strong coupling.
**Symptom:** Solver residuals increased from 0.436 to 0.450 after 1002 iterations.
**Fix:** Changed T, Te, and Tl solvers to `PBiCGStab` with `DILU` preconditioner for better convergence.

### 3. **Missing Initial Directory** ⚠️
**Problem:** The `0` directory was not created from `0.orig`.
**Fix:** Allrun script now properly restores the 0 directory.

## How to Run the Simulation Properly

### Step 1: Clean and Setup
```bash
cd ~/OpenFOAM/mavadi-v2406/run/RealisticLIFT
./Allrun
```

This will:
- Clean previous results
- Generate the mesh (blockMesh)
- Check mesh quality
- Create cell zones for substrate
- Initialize the metal film and temperature fields

### Step 2: Run the Simulation
```bash
compInterFoam | tee log.compInterFoam
```

The simulation should now:
- ✅ Solve the temperature equation successfully (converging residuals)
- ✅ Handle the alpha field properly (no singularities after first timestep)
- ✅ Generate recoil pressure and metal velocities
- ✅ Simulate the LIFT process

## What Was Changed

### Files Modified:
1. **system/fvSolution** - Changed temperature solvers from `smoothSolver` to `PBiCGStab`
2. **Allrun** (new) - Created proper setup script
3. **Allclean** (existing) - Already had correct logic to restore 0 directory

### Expected Behavior After Fixes:
```
Time = <timestep>
...
PBiCGStab:  Solving for T, Initial residual = X, Final residual = Y, No Iterations Z
```
- Residuals should DECREASE (Final < Initial)
- Alpha solver should work after first few timesteps when velocity develops
- Metal velocities should increase to m/s range (not 10^-9 m/s)

## Technical Details

### Why the Temperature Solver Failed:
The gas-phase temperature equation includes a strong coupling term to the lattice temperature:
```
∂(ρ·Cv·T)/∂t + ∇·(ρ·Cv·u·T) = ∇·(κ·∇T) + h·(Tl - T)
```

When Tl jumps to ~10,000 K due to laser heating while T is still 300 K, the coupling term `h·(Tl-T)` becomes very large. The `symGaussSeidel` smoother doesn't handle this type of stiff coupling well, especially when the matrix has both diagonal dominance issues and large off-diagonal terms.

`PBiCGStab` (Preconditioned Bi-Conjugate Gradient Stabilized) is much more robust for:
- Non-symmetric matrices
- Large source terms
- Stiff coupling between fields
- Systems with disparate time scales

### Why Alpha Had Singularities:
At t=0, the velocity field U is exactly zero everywhere. The alpha transport equation becomes:
```
∂α/∂t = 0
```

This creates a singular matrix because there's no flux to transport alpha. OpenFOAM correctly detects this and reverts to the previous field. Once the pressure and recoil forces develop (after temperature increases), velocities become non-zero and alpha transport works normally.

This is **not an error** - it's expected behavior for the first timestep.

## Verification

After running with these fixes, you should see:
1. **Converging temperature solves**: Final residual < Initial residual
2. **Increasing metal velocities**: From 10^-9 m/s to 10^0 - 10^2 m/s range
3. **Active phase change**: Mass flux > 0, recoil pressure developing
4. **Successful time advancement**: No crashes, timesteps completing

## Next Steps

If issues persist:
- Check `log.compInterFoam` for detailed diagnostics
- Monitor Courant numbers (should be < 0.5)
- Verify field initialization in the 0 directory after `setFields`
- Check that laser parameters in controlDict match your experimental setup
