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
  * `alphaBounds` (`controlDict` → `advancedInterfaceCapturing` → `alphaBounds`):
  bounds the phase fraction range contributing to recoil pressure. Defaults:
  `alphaMin 0.001; alphaMax 0.999;`.
  * `verbose` (`controlDict`): set to `true` to enable additional runtime
  logging for debugging; default `false`.
* `laserProperties` (`constant/laserProperties`): configure the femtosecond
  laser model. Common entries include `laserModel`, `pulseEnergy`,
  `pulseWidth`, `spotSize`, `focus`, `direction`, `wavelength`, optional
  scanning entries (`scanVelocity`, `pulseFrequency`, `pulseDutyCycle`) and
  `continuousLaser`.
* `twoTemperatureProperties` (`system/controlDict`): coefficients for the
  two‑temperature model—electron heat capacity `Ce`, lattice heat capacity
  `Cl`, coupling factor `G`, and the energy conservation tolerance
  `energyTol`/`energyTolerance`.
### Sample `controlDict` excerpt

```cpp
advancedInterfaceCapturing
{
    meltingTemperature   1941;    // [K]
    vaporTemperature     3200;    // [K]
    pressureScale        2e4;     // 1e4–1e5 [-], cf. Brown & Arnold (2010)
    recoilMax            5e7;     // 1e6–1e8 Pa, cf. Brown & Arnold (2010)
    recoilUpdateInterval 1;       // 1–10 steps
    alphaBounds
    {
        alphaMin 0.001;           // default 0.001
        alphaMax 0.999;           // default 0.999
    }
}
```

### Sample `laserProperties` snippet

```cpp
laserModel      gaussian;
pulseEnergy     8e-7;      // 10⁻⁸–10⁻⁶ J, cf. Brown & Arnold (2010)
pulseWidth      150e-15;   // 100–300 fs, cf. Anisimov et al. (1974)
spotSize        6e-6;      // 5–50 µm, cf. Brown & Arnold (2010)
wavelength      1.03e-6;   // 0.8–1.1 µm, cf. Brown & Arnold (2010)
absorptionCoeff 5e7;       // 10⁷–10⁸ m⁻¹, cf. Anisimov et al. (1974)
```

Typical ranges above are drawn from the LIFT and ultrafast laser literature
([Brown & Arnold, 2010](https://doi.org/10.2961/jlmn.2010.03.0001);
[Anisimov et al., 1974](https://doi.org/10.1134/1.1478536)).
The `femtosecondLaserModel` validates these entries at run time. If the
computed peak intensity

\[
I_0 = \frac{2\,E}{\pi r^2 \tau}
\]

exceeds `1e16` W m⁻² (or other bounds) the solver aborts with a
`FatalError`. For example, choosing `spotSize = 10e-6` m and
`pulseWidth = 200e-15` s limits `pulseEnergy` to roughly `3e-7` J to stay
within the threshold (`pulseEnergy = 2e-7` J yields
`I_0 ≈ 6.4e15` W m⁻²). When limits are violated, `FatalError` messages in
the log identify the offending field.

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

* Kaganov, I. M., Lifshitz, I. M., & Tanatarov, L. V. (1957). *Relaxation between
  electrons and the crystalline lattice*. Soviet Physics JETP, 4, 173–178.
* Anisimov, S. I., Kapeliovich, B. L., & Perel'man, T. L. (1974). *Electron
  emission from metal surfaces by ultrashort laser pulses*. Soviet Physics JETP,
  39(2), 375–377.
* Anisimov, S. I., Bäuerle, D., & Luk'yanchuk, B. S. (1997). *Gas-dynamics and
  recoil pressure in pulsed laser ablation of solids*. Applied Surface Science,
  96–98, 24–32. https://doi.org/10.1016/S0169-4332(96)00564-0
* Brown, M. S., & Arnold, C. B. (2010). *Fundamentals of Laser-Induced Forward
  Transfer*. Journal of Laser Micro/Nanoengineering, 5(3), 236–258.
  https://doi.org/10.2961/jlmn.2010.03.0001
* Piqué, A., & Serra, P. (2013). *Laser-Induced Forward Transfer: Fundamentals
  and Applications*. Applied Surface Science, 274, 52–58.
  https://doi.org/10.1016/j.apsusc.2012.07.029
* Rethfeld, B., Sokolowski-Tinten, K., von der Linde, D., & Anisimov, S. I.
  (2017). *Modelling ultrafast laser ablation*. Progress in Surface Science,
  92(1), 273–325. https://doi.org/10.1016/j.progsurf.2017.04.001


Refer to the standard OpenFOAM multiphase tutorials for sample cases and
adapt them for compressible and laser-induced physics scenarios.
