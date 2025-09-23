# compInterFoam LIFT configuration guide

The extended `compInterFoam` solver introduces several optional dictionaries and
switches to control the laser, two-temperature, and phase-change models. The
following tables summarize the key entries, their purpose, and in-code defaults.

## `system/controlDict`

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
| `pressureScale` | `1 Pa·s` | Scaling factor for recoil response. |
| `recoilMax` | `5e6 Pa` | Absolute recoil-pressure cap. |
| `clampRecoil` | `true` | Enables capping of recoil pressure. |
| `scaleRecoilMax` | `false` | Scales `recoilMax` with interface weighting. |
| `recoilRelax` | `1.0` | Temporal relaxation factor (`relaxFactor` supported for backward compatibility). |
| `alphaMin` / `alphaMax` | `0.001` / `0.999` | Valid phase-fraction window for recoil evaluation.

## `constant/transportProperties`

### `phaseChangeCoeffs`
Governs metal evaporation/condensation coupling supplied to the two-temperature
model.

| Entry | Default | Description |
| --- | --- | --- |
| `Tvapor` | mixture `T_vapor` | Vaporization temperature. |
| `windowWidth` | `0 K` | Half-width of smooth transition window around `Tvapor`. |
| `dtFloor` | `1e-12 s` | Minimum implicit time-step used for source limiting. |
| `relaxationRate` | – | Required source relaxation rate `[1/s]`. |
| `relaxationTime` | – | Alternative to `relaxationRate`; converted internally. |
| `maxSource` (`minCoefficient`) | See note | Maximum allowable source term. Uses `transportProperties.phaseChangeMaxSourceDefault` if unspecified (fallback `1e7`). |
| `phaseChangeMetalCutoff` | `metalFractionCutoff` | Minimum metal fraction for source activation (falls back to `twoTemperatureProperties`). |
| `onlyAboveVapor` | `false` | Disable phase-change source below `Tvapor`. |
| `activationTime` | – | List of `(start stop)` time windows when the source is active. |

### `massTransferCoeffs`
Optional bulk mass-transfer limiter for the gas phase.

| Entry | Default | Description |
| --- | --- | --- |
| `rateMax` | `-1` | Maximum `dg/dt` magnitude. Negative disables limiting. |
| `tStart`, `tEnd` | – | Parallel lists defining activation windows. Missing or empty lists keep the limiter active at all times.

## Notes
* All diagnostics controlled by `verbose` are now emitted by the master MPI
  rank only to prevent duplicate messages in parallel runs.
* Defaults above reflect the values compiled into the solver; supply explicit
  entries in your dictionaries to override them as needed.
