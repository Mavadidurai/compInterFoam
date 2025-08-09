# compInterFoam

compInterFoam is an OpenFOAM solver for two compressible, non-isothermal immiscible fluids using a volume-of-fluid (VOF) phase-fraction based interface-capturing approach. It solves a single momentum equation for the mixture and supports either mixture or two-phase transport modelling. This fork adds femtosecond laser modelling, two-temperature physics, and advanced interface capturing for Laser-Induced Forward Transfer (LIFT).

## Build

This code requires a working [OpenFOAM](https://www.openfoam.com/) installation. After sourcing the OpenFOAM environment, compile the solver with:

```bash
source /path/to/OpenFOAM/etc/bashrc
wmake
```

The resulting executable will be placed in `$FOAM_USER_APPBIN`.

## Run

Create or load an OpenFOAM case for the solver and execute:

```bash
compInterFoam
```

## Citation

If you use this solver in academic work, please cite the repository. A suggested citation is:

```
compInterFoam developers (2024). compInterFoam: solver for compressible immiscible flows with LIFT modelling. GitHub repository. https://github.com/OWNER/compInterFoam
```

