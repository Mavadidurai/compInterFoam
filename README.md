# compInterFoam LIFT configuration guide

The extended `compInterFoam` solver introduces several optional dictionaries and
switches to control the laser, two-temperature, and phase-change models. The
following tables summarize the key entries, their purpose, and in-code defaults.
## Building the solver

The code relies on OpenFOAM's `wmake` build system. To avoid the
`bash: command not found: wmake` error, first source the OpenFOAM
environment (replace the path with the one that matches your installation):

```bash
source /opt/openfoam/etc/bashrc
```

After sourcing the environment you can compile the solver with the usual
command:

```bash
wmake
```

If you are unsure whether the environment is set up correctly you can also run
the helper script in the repository root:

```bash
./wmake
```

It will delegate to the system `wmake` when available or print a clear
instruction when the OpenFOAM environment is missing.
## `system/controlDict`
The reference case in [`Latest/system/controlDict`](Latest/system/controlDict)
keeps the base time step fixed with `adjustTimeStep no`, so the Courant-number
limits (`maxCo`, `maxAlphaCo`, `maxThermalCourant`) act as documentation only
unless you explicitly set `adjustTimeStep yes`. Output is controlled with
`writeControl runTime`, where `writeInterval` is specified in seconds of
simulated time (`5e-13` corresponds to 0.5 ps in the default case).

### `verbose`
* **Type:** `Switch`
* **Default:** `false`
* **Effect:** Enables detailed diagnostics for laser, two-temperature, and
  interface-capturing models. Messages are written by the master process only.

### `twoTemperatureProperties`
Controls the two-temperature model and related mixture properties. Frequently
used entries include:

| Entry | Default | Description |
| --- | --- | --- |
| `Ce` | `210 J/m^3/K` | Electron heat capacity (value or `Function1`). |
| `Cl` | `2.3e6 J/m^3/K` | Lattice heat capacity. Also cached by `twoPhaseMixtureThermo` via `setClTTM`. |
| `G` | `5e17 W/m^3/K` | Electron-phonon coupling. |
| `De` | `1e-4 m^2/s` | Electron thermal diffusivity. |
| `gasMetalExchangeCoeff` | `5e17 W/m^3/K` | Gas/metal interfacial coupling. |
| `energyTolerance` (`energyTol`) | `0.1` | Relative energy-change tolerance checked after each time step. |
| `minTe` / `maxTe` | `300 K` / `3500 K` | Bounds applied to `Te` and `Tl`. |
| `metalFractionFloor` | `1e-6` | Minimum metal fraction considered for phase-change coupling and temperature bounds. |
| `metalFractionCutoff` | `metalFractionFloor` | Cells with metal fraction below this value are excluded from phase-change coupling. |
| `metalAmbientBlendWidth` | `1e-3` | Width of blending region for ambient fallback temperature. |
| `electronSubCycles` | `1` | Planned electron sub-iterations per time step.
| `maxElectronDeltaT` | `Δt` | Optional cap on electron sub-cycle time-step.
| `minElectronSubCycles` / `maxElectronSubCycles` | `1` / unlimited | Bounds on
sub-cycle count.
| `temperatureDiagnostics` | `verbose` | Print temperature statistics when true.
| `energyDiagnostics` | `false` | Enables extra energy reporting.
| `energyAudit` | `verbose` | Enables detailed energy conservation log.

### `advancedInterfaceCapturing`
Defines parameters for the recoil-pressure helper used when
`useAdvancedInterfaceCapturing` (default `true`) is enabled.

| Entry | Default | Description |
| --- | --- | --- |
| `meltingTemperature` | mixture `T_melt` | Lower activation bound. |
| `vaporTemperature` | mixture `T_vapor` | Upper activation bound. |
| `phaseChangeTempOffset` | `0 K` | Shifts temperature threshold. |
| `recoilTempOffset` | `0 K` | Extra offset applied before recoil activation. Must be non-negative and less than `vaporTemperature`. |
| `pressureScale` | `1 Pa·s` (`[1 -1 -1 0 0 0 0]`) | Scaling factor for recoil response; supply as a Pa·s dimensioned scalar. |
| `recoilMax` | `5e6 Pa` | Absolute recoil-pressure cap. |
| `clampRecoil` | `true` | Enables capping of recoil pressure. |
| `scaleRecoilMax` | `false` | Scales `recoilMax` with interface weighting. |
| `recoilRelax` | `1.0` | Temporal relaxation factor (`relaxFactor` supported for backward compatibility). |
| `boltzmannConstant` | `1.38e-23 J/K` (`[1 2 -2 -1 0 0 0]`) | Boltzmann constant used in recoil pressure; must be positive. |
| `vaporParticleMass` | `7.95e-26 kg` (`[1 0 0 0 0 0 0]`) | Mass of the evaporated species for the recoil model; must be positive. |
| `momentumAccommodationCoeff` | `0.18` | Momentum accommodation coefficient applied to the recoil pressure; constrained to `[0, 1]`. |
| `alphaMin` / `alphaMax` | `0.001` / `0.999` | Valid phase-fraction window for recoil evaluation.

## `system/controlDict`

### `phaseChangeCoeffs`
Governs metal evaporation/condensation coupling supplied to the two-temperature
model. Moving this dictionary into `controlDict` keeps timing adjustments close
to the main run controls.

| Entry | Default | Description |
| --- | --- | --- |
| `Tvapor` | mixture `T_vapor` | Vaporization temperature. |
| `windowWidth` | `0 K` | Half-width of smooth transition window around `Tvapor`. |
| `dtFloor` | `1e-12 s` | Minimum implicit time-step; clamps implicit relaxation to `1/dtFloor`. |
| `relaxationRate` | – | Required source relaxation rate `[1/s]`. |
| `relaxationTime` | – | Alternative to `relaxationRate`; converted internally. |
| `maxSource` (`minCoefficient`) | See note | Maximum allowable source term. Uses `transportProperties.phaseChangeMaxSourceDefault` if unspecified (fallback `1e7`). |
| `phaseChangeMetalCutoff` | `metalFractionCutoff` | Minimum metal fraction for source activation (falls back to `twoTemperatureProperties`). |
| `onlyAboveVapor` | `false` | Disable phase-change source below `Tvapor`. |
| `activationTime` | – | List of `(start stop)` time windows when the source is active. |
| `gasConstant` | `173.7 J/(kg·K)` | Hertz–Knudsen gas constant. Supply as a dimensioned scalar with units `[0 2 -2 -1 0 0 0]`. |
| `evaporationCoeff` | `0.18` | Evaporation/condensation accommodation coefficient (dimensionless, > 0). |
| `evapRelaxationTime` | `1e-12 s` | Sets the implicit relaxation time applied to the net evaporation source (`> 0`). |
| `alphaMin` / `alphaMax` | `0.01` / `0.99` | Bounding phase-fraction window for evaluating evaporation (`0 ≤ alphaMin < alphaMax ≤ 1`). |

### `massTransferCoeffs`
Optional bulk mass-transfer limiter for the gas phase. Relocating this
dictionary into `controlDict` puts its activation windows alongside other
timing settings.

| Entry | Default | Description |
| --- | --- | --- |
| `rateMax` | `-1` | Maximum `dg/dt` magnitude. Negative disables limiting. |
| `tStart`, `tEnd` | – | Parallel lists defining activation windows. Missing or empty lists keep the limiter active at all times. |

### `laserStartTime` / `laserEndTime`
Top-level entries that define when the femtosecond laser source is active.
Entries in `controlDict` override values supplied in `constant/laserProperties`.
Single-pulse runs may additionally supply `pulseCenterTime` in
`constant/laserProperties` to delay the Gaussian peak; values outside the
window are clamped to the configured start/end times, and omitting the entry
retains the legacy near-start placement. When the keyword is absent the solver
centres the Gaussian at `laserStartTime + min(0.5·(laserEndTime -
laserStartTime), 3σ)` so the energy lands within a few pulse widths of the
trigger. Providing an explicit `pulseCenterTime` simply shifts this centre—it
does not change the deposited energy, only when the peak arrives during the
`[laserStartTime, laserEndTime]` window.
## Notes
* All diagnostics controlled by `verbose` are now emitted by the master MPI
  rank only to prevent duplicate messages in parallel runs.
* The femtosecond-laser pulse-energy check scales its expectation by the
  Gaussian beam area that intersects the mesh, so clipping the spot no longer
  triggers spurious mismatch warnings.
* Defaults above reflect the values compiled into the solver; supply explicit
  entries in your dictionaries to override them as needed.
* The femtosecond laser accepts `spatialIntegrationMode` to control axial
  chord estimation: `exact` (default) scans cell vertices, while `centerline`
  approximates the span from `V^(1/3)` for faster but less precise deposition.
## `system/fvSolution:compInterFoamCoeffs`
The solver re-reads the `compInterFoamCoeffs` sub-dictionary from
`fvSolution` whenever it needs clamp limits or solver fallbacks. The active
entries and their in-code defaults are:

| Entry | Default | Used in |
| --- | --- | --- |
| `maxPressureGradient` | `GREAT` | Caps `∇p_rgh` used in the momentum predictor and zeroes non-finite components. |
| `maxUEqnVelocity` | `GREAT` | Limits `U` before computing `phi = fvc::flux(U)` to prevent runaway mass flux. |
| `maxVelocity` | `500` | Bounds the velocity magnitude after the pressure correction. |
| `minUEqnDiag` | `1e-9` | Floors the diagonal before inverting `UEqn` to form `rAU`. |
| `enableRAUClamp` | `false` | Enables explicit `rAU`/`rAUf` bounding. |
| `minRAU`, `maxRAU` | `1e-10`, `GREAT` | Cell-wise `rAU` clamp (only when `enableRAUClamp` is `true`). |
| `minRAUf`, `maxRAUf` | `1e-14`, `GREAT` | Face-wise `rAUf` clamp (only when `enableRAUClamp` is `true`). |
| `pressureClamp` | `false` | Enables bounding of `p` and `p_rgh`. |
| `minPressure`, `maxPressure` | `-GREAT`, `GREAT` | Limits applied when `pressureClamp` is `true`. |

Any other keys placed in `compInterFoamCoeffs` are ignored by the current
code, so make sure to set the values above if you need runtime clamping.
