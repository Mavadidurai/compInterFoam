# compInterFoam

## Project purpose and key features

compInterFoam is an OpenFOAM solver for two compressible, non-isothermal
immiscible fluids using a volume-of-fluid (VOF) phase-fraction based
interface capturing approach. It extends the standard `compressibleInterFoam`
solver with femtosecond laser modeling, two-temperature physics, and
advanced interface capturing geared toward Laser-Induced Forward Transfer
(LIFT) studies.

The two-temperature model enforces energy conservation using an
`energyTol`/`energyTolerance` parameter specified in case dictionaries.

## Required OpenFOAM version and dependencies

* OpenFOAM v10 or later
* Standard OpenFOAM build dependencies (C++ compiler, wmake, MPI)
* Libraries bundled with this repository: `femtosecondLaserModel`,
  `twoTemperatureModel`, `advancedInterfaceCapturing`,
  `twoPhaseMixtureThermo`, and `DimensionValidator`

Ensure that your OpenFOAM environment variables (e.g. `WM_PROJECT_DIR`,
`FOAM_USER_APPBIN`) are configured.

## Build instructions

Inside a configured OpenFOAM environment:

```bash
# 1. Source your OpenFOAM installation
source /path/to/OpenFOAM-<version>/etc/bashrc

# 2. Confirm the build tool is available
command -v wmake

# 3. Inspect the user application directory derived from the environment
echo $FOAM_USER_APPBIN

# 4. Move to the application directory
cd $FOAM_USER_APPBIN
# Clone this repository if not already present
# git clone <repository-url> compInterFoam
cd compInterFoam

# 5. Build bundled libraries (required before compiling the solver)
wmake libso DimensionValidator
wmake libso twoPhaseMixtureThermo
wmake libso femtosecondLaserModel
wmake libso twoTemperatureModel
wmake libso advancedInterfaceCapturing

# 6. Build the solver
wmake
# or build from the OpenFOAM tree
wmake applications/solvers/compInterFoam
```

The executable will be placed in `$FOAM_USER_APPBIN`.

## Key dictionary entries

Some solver features rely on additional entries in case dictionaries:

* `pressureScale` (`controlDict` → `advancedInterfaceCapturing`):
  nondimensional scaling for evaporation-induced recoil pressure.
* `recoilMax` (`controlDict` → `advancedInterfaceCapturing`): cap on recoil
  pressure, default `5e6` Pa.
* `recoilUpdateInterval` (`controlDict` → `advancedInterfaceCapturing`):
  number of time steps between recoil field updates, default `5`.
* `useAdvancedInterfaceCapturing` (`controlDict`): enable the recoil-based
  interface capturing model.
* `laserProperties` (`constant/laserProperties`): configure the femtosecond
  laser model. Common entries include `laserModel`, `pulseEnergy`,
  `pulseWidth`, `spotSize`, `focus`, `direction`, `wavelength`, optional
  scanning entries (`scanVelocity`, `pulseFrequency`, `pulseDutyCycle`) and
  `continuousLaser`.
* `twoTemperatureProperties` (`system/controlDict`): coefficients for the
  two‑temperature model—electron heat capacity `Ce`, lattice heat capacity
  `Cl`, coupling factor `G`, and the energy conservation tolerance
  `energyTol`/`energyTolerance`.

## Relevant source files

* `compInterFoam.C` – solver main.
* `alphaEqnSubCycle.H`, `UEqn.H`, `pEqn.H`, `TEqn.H` – phase, momentum,
  pressure and temperature equations.
* `advancedInterfaceCapturing.C/H` – VOF interface capturing with recoil
  pressure; defines `recoilMax`.
* `femtosecondLaserModel.C/H` – laser energy deposition.
* `twoTemperatureModel.C/H` – electron/phonon energy coupling.
* `twoPhaseMixtureThermo.C/H` – thermophysical properties.
* `DimensionValidator.C/H` – run-time dimension checks.

## Example usage

Run the solver on a case directory:

```bash
compInterFoam -case path/to/case
```

For parallel runs:

```bash
decomposePar -case path/to/case
mpirun -np <N> compInterFoam -parallel -case path/to/case
```

### Provided `TestCase`

The repository includes a minimal LIFT setup in `TestCase`. Run it with:

1. **Clean the case** *(optional)*:
   ```bash
   cd TestCase
   ./Allclean
   ```
2. **Generate the mesh**:
   ```bash
   blockMesh
   ```
3. **Set up initial fields**:
   ```bash
   cp -r 0.orig 0   # or use: setFields
   ```
4. **Decompose for parallel runs** *(optional)*:
   ```bash
   decomposePar
   ```
5. **Launch the simulation**:
   ```bash
   compInterFoam                        # serial
   mpirun -np <N> compInterFoam -parallel  # parallel
   ```
6. **Post-process** *(optional)*:
   ```bash
   reconstructPar
   foamToVTK
   paraFoam
   ```

Log files can be captured by redirecting command output, e.g. `blockMesh > log.blockMesh` and `compInterFoam | tee log.compInterFoam`. Output fields are written to time directories (e.g. `0`, `0.001`, …) in the case folder, or to `processor*/` subdirectories during parallel runs. After `reconstructPar`, fields are merged back into the main case directory for visualization.

## References

* Brown, M. S., & Arnold, C. B. (2010). *Fundamentals of Laser-Induced Forward
  Transfer*. Journal of Laser Micro/Nanoengineering, 5(3), 236–258.
  https://doi.org/10.2961/jlmn.2010.03.0001
* Anisimov, S. I., Kapeliovich, B. L., & Perel'man, T. L. (1974). *Electron
  emission from metal surfaces by ultrashort laser pulses*. Soviet Physics JETP,
  39(2), 375–377.


Refer to the standard OpenFOAM multiphase tutorials for sample cases and
adapt them for compressible and laser-induced physics scenarios.
