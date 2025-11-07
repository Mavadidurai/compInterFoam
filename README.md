# compInterFoam LIFT configuration guide

The extended `compInterFoam` solver introduces several optional dictionaries and
switches to control the laser, two-temperature, and phase-change models. The
following tables summarize the key entries, their purpose, and in-code defaults.
## Building the solver
## RealisticLIFT restart policy
The `RealisticLIFT/system/controlDict` case disables `purgeWrite` so every
checkpoint directory is preserved until post-processing. This keeps restart
points consistent with the archived log. If you must re-enable purging for
storage reasons, also change `startFrom` to `startTime` so the solver always
restarts from a clean state rather than a potentially pruned intermediate time.
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
| `gasConstant` | required | Hertz-Knudsen gas constant. Supply a material-appropriate dimensioned scalar with units `[0 2 -2 -1 0 0 0]`. |
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
| `maxKineticEnergyDensity` | `1e12` | Aborts the run when `0.5*rho*|U|^2` exceeds the ceiling before solving `UEqn`. |
| `maxUEqnVelocity` | `GREAT` | Limits `U` before computing `phi = fvc::flux(U)` to prevent runaway mass flux. |
| `maxVelocity` | `500` | Bounds the velocity magnitude after the pressure correction. |
| `velocityAlphaThreshold` | `0.01` | Fraction of phase-1 (`alpha1`) used when reporting the post-correction velocity diagnostic. Speeds in cells below this volume fraction are treated as gas/plume motion. |
| `minUEqnDiag` | `1e-9` | Floors the diagonal before inverting `UEqn` to form `rAU`. |
| `enableRAUClamp` | `false` | Enables explicit `rAU`/`rAUf` bounding. |
| `minRAU`, `maxRAU` | `1e-10`, `GREAT` | Cell-wise `rAU` clamp (only when `enableRAUClamp` is `true`). |
| `minRAUf`, `maxRAUf` | `1e-14`, `GREAT` | Face-wise `rAUf` clamp (only when `enableRAUClamp` is `true`). |
| `pressureClamp` | `false` | Enables bounding of `p` and `p_rgh`. |
| `minPressure`, `maxPressure` | `-GREAT`, `GREAT` | Limits applied when `pressureClamp` is `true`. |

Any other keys placed in `compInterFoamCoeffs` are ignored by the current
code, so make sure to set the values above if you need runtime clamping.

## Solver field symbols
The solver allocates the following named fields and helper scalars during
startup. Units follow the dimensions supplied when the objects are
constructed in `createFields.H` and the equation includes; the short
descriptions summarise their role in the femtosecond LIFT workflow.

| Symbol | Units | Short description |
| --- | --- | --- |
| `trDeltaT` | s⁻¹ | Inverse local time-step used when LTS is active (`fv::localEulerDdt::rDeltaTName`). 【F:createFields.H†L4-L34】 |
| `trSubDeltaT` | s⁻¹ | Inverse sub-cycled time-step companion to `trDeltaT` for LTS momentum updates. 【F:createFields.H†L4-L41】 |
| `p_rgh` | Pa | Pressure relative to the hydrostatic column, read from the latest time directory. 【F:createFields.H†L46-L68】 |
| `U` | m s⁻¹ | Mixture velocity field solved by the momentum predictor. 【F:createFields.H†L70-L93】 |
| `phi` | m³ s⁻¹ | Face volumetric flux obtained from `fvc::flux(U)` for use in transport terms. 【F:createFields.H†L95-L118】 |
| `recoilTractionf` | N m⁻³ | Face-based recoil traction initialised to zero and later populated by the recoil model. 【F:createFields.H†L120-L143】 |
| `T` | K | Single-temperature field (initialised to 300 K) that seeds metal/gas temperatures. 【F:createFields.H†L145-L170】 |
| `gasMetalHeatFlux` | W m⁻³ | Volumetric accumulator for Kapitza-style heat exchange between gas and metal. 【F:createFields.H†L172-L195】 |
| `alpha1`, `alpha2` | – | Phase-fraction fields for metal (`alpha.metal`) and gas (`alpha.gas`) retrieved from the mixture. 【F:createFields.H†L324-L341】 |
| `rho` | kg m⁻³ | Mixture density assembled from `alpha1*rho1 + alpha2*rho2`. 【F:createFields.H†L353-L369】 |
| `g` | m s⁻² | Uniform gravitational acceleration vector read from `constant/g`. 【F:createFields.H†L381-L401】 |
| `hRef` | m | Reference height used when reconstructing hydrostatic pressure. 【F:createFields.H†L403-L412】 |
| `ghRef` | m² s⁻² | Reference gravitational potential (`g·hRef`) preserving the gravity direction. 【F:createFields.H†L414-L427】 |
| `gh`, `ghf` | m² s⁻² | Cell- and face-based gravitational potentials used when forming `p` and flux corrections. 【F:createFields.H†L429-L436】 |
| `p` | Pa | Absolute pressure field defined as `p_rgh + rho*gh` for output and diagnostics. 【F:createFields.H†L438-L454】 |
| `pRefCell`, `pRefValue` | –, Pa | Reference cell index and pressure used to stabilise the pressure solve. 【F:createFields.H†L456-L463】 |
| `rhoPhi` | kg s⁻¹ | Mass flux computed from interpolated density and volumetric flux. 【F:createFields.H†L465-L478】 |
| `dgdt` | kg m⁻³ s⁻¹ | Volumetric mass-transfer rate supplied by the phase-change thermo model. 【F:createFields.H†L480-L488】 |
| `alphaPhi10` | m³ s⁻¹ | Face flux of the metal phase used in the compressive VOF transport. 【F:createFields.H†L490-L520】 |
| `nAlphaSubCycles`, `nAlphaCorr` | count | Controls for the alpha equation sub-cycling and correction passes read from `fvSolution`. 【F:createFields.H†L522-L539】 |
| `MULESCorr`, `alphaApplyPrevCorr` | Boolean | Switches governing the compressive VOF algorithm behaviour. 【F:createFields.H†L541-L551】 |
| `icAlpha`, `scAlpha` | – | Interface-compression and smoothing coefficients applied in the alpha transport. 【F:createFields.H†L553-L561】 |
| `recoilPressure` | Pa | Legacy recoil-pressure field allocated when advanced interface capturing is disabled. 【F:createFields.H†L642-L688】 |
| `K` | m² s⁻² | Specific kinetic energy (`0.5·|U|²`) updated after the momentum correction for diagnostics. 【F:createFields.H†L571-L579】【F:UEqn.H†L602-L610】 |
| `rAU` | aP⁻¹ | Reciprocal diagonal of the momentum matrix (`1/A`) with dimensions inherited from `UEqn.A()`. 【F:UEqn.H†L612-L624】 |
| `rhoPhi` (post-momentum) | kg s⁻¹ | Updated mass flux after the velocity-limiting pass before the pressure solve. 【F:UEqn.H†L646-L656】 |
| `p_rgh` (corrected) | Pa | Pressure-relative field repeatedly refreshed in the pressure equation include. 【F:pEqn.H†L172-L207】 |
| `magUAll`, `magULiquid` | m s⁻¹ | Velocity magnitude diagnostics for mixture and metal-dominated regions. 【F:pEqn.H†L378-L405】 |
| `phaseChangeSource`, `phaseChangeRelaxCoeff` | K s⁻¹, s⁻¹ | Temperature-source term and its implicit relaxation coefficient imported from the thermo model. 【F:compInterFoam.C†L724-L734】【F:twoPhaseMixtureThermo.C†L70-L101】 |
