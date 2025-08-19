# compInterFoam

## Project purpose and key features

compInterFoam is an OpenFOAM solver for two compressible, non-isothermal
immiscible fluids using a volume-of-fluid (VOF) phase-fraction based
interface capturing approach. It extends the standard `compressibleInterFoam`
solver with femtosecond laser modeling, two-temperature physics, and
advanced interface capturing geared toward Laser-Induced Forward Transfer
(LIFT) studies.

The two-temperature model enforces energy conservation using an
`energyTolerance` parameter specified in case dictionaries.

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
cd $FOAM_USER_APPBIN
# Clone this repository if not already present
# git clone <repository-url> compInterFoam
cd compInterFoam
wmake
```

The executable will be placed in `$FOAM_USER_APPBIN`.

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

Refer to the standard OpenFOAM multiphase tutorials for sample cases and
adapt them for compressible and laser-induced physics scenarios.
