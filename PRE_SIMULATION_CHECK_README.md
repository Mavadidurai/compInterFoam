# Pre-Simulation Validation Script

## Overview

The `preSimulationCheck.py` script performs comprehensive validation of compInterFoam simulation cases before running, helping to catch configuration errors early and ensure simulation stability.

## Usage

```bash
# Check current directory
python3 preSimulationCheck.py

# Check specific case directory
python3 preSimulationCheck.py /path/to/case

# Example: Check TEST1
python3 preSimulationCheck.py ./TEST1
```

## What It Checks

### 1. File Structure Validation
- **Required directories**: `system/`, `constant/`, `0/` (or `0.orig/`)
- **System files**: `controlDict`, `fvSchemes`, `fvSolution`, `blockMeshDict`
- **Constant files**: `thermophysicalProperties`, `transportProperties`, `g`, `laserProperties`
- **Initial fields**: `p`, `U`, `T`, `alpha.air`, `alpha.metal`, etc.

### 2. Timing Parameters (Critical for Femtosecond Simulations)

#### Basic Time Settings
- **startTime**: Must be non-negative
- **endTime**: Must be positive and greater than startTime
- **deltaT**: Initial timestep must be positive
- **Simulation duration**: Validates reasonable duration for femtosecond simulations

#### Timestep Validation
- **Number of timesteps**: Estimates total steps and warns if too many (>1M) or too few (<10)
- **Timestep resolution**: Checks if pulse can be resolved (>5 timesteps per laser pulse)

#### Adaptive Timestepping
- **adjustTimeStep**: If enabled, validates:
  - **maxCo** (velocity Courant): Warns if >1.0 (instability) or <0.01 (too conservative)
  - **maxAlphaCo** (interface Courant): Errors if >1.0 (interface smearing), warns if <0.01
  - **maxDeltaT**: Upper timestep limit validation
  - **maxThermalCourant**: Thermal diffusion stability check
  - **Timestep range**: Warns if maxDeltaT/deltaT > 1000 (abrupt changes possible)

#### Write Control
- **writeControl**: Validates against standard options (timeStep, runTime, adjustableRunTime, cpuTime, clockTime)
- **writeInterval**:
  - For timeStep mode: Must be ≥1, warns if =1 (writes every step)
  - Estimates number of output times and warns if >1000 or <5
- **purgeWrite**: Notes if enabled (keeps only last N time directories)
- **runTimeModifiable**: Whether parameters can be changed during simulation
- **timePrecision**: Warns if <6 (naming issues for small timesteps)
- **writePrecision**: Warns if <8 (accuracy loss in output)

### 3. Laser Properties (Femtosecond LIFT Specific)

- **pulseEnergy**:
  - Must be positive
  - Warns if >1 mJ (very high) or <1 nJ (too low for material removal)

- **pulseWidth**:
  - Must be positive
  - Warns if >1 ps (not femtosecond) or <0.1 fs (unrealistic)
  - Validates temporal resolution vs simulation timestep

- **wavelength**:
  - Must be positive
  - Warns if outside 200-2000 nm range

- **spotSize/beamRadius**:
  - Must be positive
  - Warns if >100 μm (large for LIFT)

- **focusPoint**: Validates (x,y,z) format

- **absorptionCoefficient**:
  - Must be positive
  - Calculates penetration depth
  - Compares with film thickness

- **reflectivity**:
  - Must be in range [0,1]
  - Warns if >90% (little energy absorbed)

- **filmThickness**:
  - Validates ratio to penetration depth
  - Warns if film << penetration depth (energy passes through)

### 4. Thermophysical Properties

#### For Each Phase (air, metal):
- **Density (rho)**: Must be positive
- **Specific heat (Cp)**: Must be positive
- **Thermal conductivity (kappa/k)**: Cannot be negative
- **Viscosity (mu/nu)**: Cannot be negative

#### Transport Properties:
- **Surface tension (sigma)**:
  - Cannot be negative
  - Warns if >2.0 N/m (very high)

### 5. Numerical Settings

#### fvSchemes:
- **Time discretization (ddtSchemes)**:
  - Validates default scheme exists
  - Notes if using Euler (1st order) vs higher-order schemes

- **Gradient schemes (gradSchemes)**: Validates default exists

- **Divergence schemes (divSchemes)**:
  - Checks for alpha divergence scheme (critical for interface tracking)
  - Warns if missing (interface may smear)

#### fvSolution:
- **PIMPLE settings**:
  - `nOuterCorrectors`: Warns if <3 (may not converge)
  - `nCorrectors`: Pressure corrector iterations
  - `nNonOrthogonalCorrectors`: For mesh non-orthogonality

- **Relaxation factors**:
  - Validates all factors are in range (0,1]
  - Warns if outside valid range

### 6. Boundary Conditions

- **Field files**: Checks that p, U, T files exist and are not empty
- **Alpha fields**:
  - Requires at least 2 alpha fields for two-phase flow
  - Warns if only 1 or none found

### 7. Advanced Physics Parameters

- **Two-temperature model**:
  - `electronHeatCapacityCoeff`: Electron heat capacity coefficient
  - `electronPhononCoupling`: Electron-phonon coupling strength

- **Phase change**:
  - `TSat`: Saturation temperature validation

- **Gravity**:
  - Validates gravity vector format (x,y,z)
  - Calculates magnitude
  - Notes if near-zero or standard Earth gravity

- **Dynamic mesh**:
  - Checks if mesh motion is enabled
  - Notes solver type

## Output Format

The script provides color-coded output:

- **🔴 ERRORS** (Red): Critical issues that will prevent simulation from running correctly
- **🟡 WARNINGS** (Yellow): Issues that may cause problems or inefficiencies
- **🔵 INFORMATION** (Cyan): Informational messages about configuration
- **🟢 PASSED** (Green): Successful validations

### Exit Codes
- **0**: Validation passed (with or without warnings)
- **1**: Validation failed (errors found)

## Example Output

```
================================================================================
CompInterFoam Pre-Simulation Validation
================================================================================
Case directory: /home/user/compInterFoam/TEST1

[1/8] Checking file structure...
[2/8] Loading dictionary files...
[3/8] Validating timing parameters...
[4/8] Validating laser properties...
[5/8] Validating thermophysical properties...
[6/8] Validating numerical settings...
[7/8] Validating boundary conditions...
[8/8] Validating advanced physics parameters...

================================================================================
VALIDATION SUMMARY
================================================================================

WARNINGS (9):
  [W01] writeControl 'writeTime' may not be standard
  [W02] focusPoint not defined - laser position unclear
  ...

INFORMATION (3):
  [I01] Adaptive timestepping enabled - estimated ~500,000 steps
  ...

PASSED CHECKS: 53

RESULT: VALIDATION PASSED WITH WARNINGS - Review warnings before proceeding
```

## Best Practices

### Before Running Simulation:
1. Always run this validation script
2. **Address all ERRORS** - these will cause simulation failures
3. **Review all WARNINGS** - these indicate potential problems
4. Pay special attention to:
   - Timestep settings (especially for femtosecond simulations)
   - Courant number limits
   - Write interval (disk space considerations)
   - Laser energy and pulse parameters

### Common Issues:

#### Too Many Timesteps
```
WARNING: Very large number of timesteps (2,000,000) - simulation may take long time
```
**Solution**: Increase `deltaT` or reduce `endTime`, or ensure `adjustTimeStep` is enabled

#### Interface Smearing
```
WARNING: No alpha divergence scheme found - interface may smear
```
**Solution**: Add proper alpha compression schemes in `fvSchemes/divSchemes`

#### Unstable Courant Numbers
```
ERROR: maxAlphaCo (1.5) > 1.0 will cause interface smearing
```
**Solution**: Reduce `maxAlphaCo` to ≤1.0 in `controlDict`

#### Too Many Output Files
```
WARNING: Very large number of outputs (5000) - consider increasing writeInterval
```
**Solution**: Increase `writeInterval` to reduce disk usage

#### Poor Temporal Resolution
```
WARNING: Only 3.2 timesteps per pulse - increase resolution (reduce deltaT)
```
**Solution**: Reduce `deltaT` to have at least 5-10 timesteps per laser pulse

## Integration with Workflow

### Recommended Workflow:
```bash
# 1. Set up case
./Allrun.pre

# 2. Run validation
python3 preSimulationCheck.py .

# 3. If errors, fix and re-validate
# Edit system/controlDict or other files
python3 preSimulationCheck.py .

# 4. Run simulation
./Allrun
```

### In Allrun Script:
```bash
#!/bin/bash
cd ${0%/*} || exit 1

# Run pre-simulation checks
python3 ../preSimulationCheck.py . || {
    echo "Pre-simulation validation failed!"
    exit 1
}

# Proceed with simulation...
runApplication blockMesh
runApplication decomposePar
runParallel $(getApplication)
```

## Customization

The script can be extended by modifying:

- `check_timing_parameters()`: Add custom timing checks
- `check_laser_properties()`: Add laser-specific validations
- `check_thermophysical_properties()`: Add material property checks
- `check_numerical_settings()`: Add scheme/solver validations

## Limitations

- Simple dictionary parser (may not handle complex nested structures)
- Does not validate mesh quality (use `checkMesh` for that)
- Does not check for typos in field names
- Does not validate parallelization settings
- Cannot predict actual simulation stability (runtime issue)

## See Also

- OpenFOAM User Guide: https://www.openfoam.com/documentation/user-guide
- `checkMesh`: For mesh quality validation
- `foamDictionary`: For querying dictionary values
- `IMPLEMENTATION_GUIDE.md`: For case setup instructions
