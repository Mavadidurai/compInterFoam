# LIFT Physics & Numerics Fixes

**Date:** 2026-07-08
**Branch:** `claude/lift-code-review-debug-h6hfzr`
**Scope:** solver sources + all four cases (TEST1, TEST2, TEST3, TestCase)

This change set repairs the defects found in the full code review so that the
LIFT chain — fs-laser absorption → two-temperature heating → Hertz–Knudsen
evaporation → recoil pressure → film acceleration/ejection → vapour-driven
expansion — is internally consistent. Fixes are ordered by impact; each entry
states what was wrong, what changed, and what downstream outputs it affects.

---

## 1. Pressure–velocity coupling (`pEqn.H`) — CRITICAL

**Was:** the matrix was assembled as `laplacian(rAUf,p_rgh) == div(phiHbyA)+S`
while the updates kept the stock idiom `phi = phiHbyA + pEqn.flux()`. For this
arrangement `flux()` has the opposite sign, so every corrector pushed flux *up*
the pressure gradient. This is the direct cause of the GAMG stall, the 5 GPa
clamp hits, and the velocity/recoil inconsistency documented in
`CONVERGENCE_ISSUE_DIAGNOSIS.md`. In addition `pEqn.relax(0.2)` decoupled
`flux()` from the solved pressure (mass-conservation error every corrector).

**Now:** stock `compressibleInterFoam` arrangement
(`ddt + div(phiHbyA) - laplacian == (rho1/rho2 - 1)*dgdt`), no matrix
relaxation, and a mixture acoustic-compliance term
`(alpha_i*psi_i/rho_i)*ddt(p_rgh)` so the gas/vapour phase responds
compressibly (Boussinesq metal has psi = 0).

**Affects:** everything downstream — pressure level, jet velocity, time-step
evolution. Old runs' pressure/velocity fields are not comparable.

## 2. Gas temperature no longer overwritten (`compInterFoam.C`)

**Was:** `T = ttm.Tl()` after the thermal loop replaced the entire solved gas
temperature with the lattice field each outer iteration; in gas cells Tl is
ambient-clamped, so the plume stayed at 300 K, the gas EOS never expanded, and
the 1400-line `TEqn.H` had no lasting effect.

**Now:** removed. `TEqn.H` already projects `T → Tl` in metal-dominated cells
(alpha1 > 0.9); interface/gas cells keep their solved temperature.

**Affects:** vapour plume dynamics, gas density (perfectGas: rho = p/RT),
pressure in the gap, receiver-side heat flux.

## 3. Interface-localized evaporation (`twoPhaseMixtureThermo.C/.H`)

**Was:** every cell in the alpha window (which includes bulk alpha = 1 cells,
since `alphaMax = 1` disables the upper gate) evaporated through its full
volume, using a cell-thickness heuristic — i.e. volumetric boiling of the
whole hot film instead of surface evaporation.

**Now:** cells need a resolved interface (`|grad(alpha)|*cellSize >=
interfaceGradientCutoff`, default 1e-3, new optional entry in
`phaseChangeCoeffs`), and the volumetric source is
`q_vol = j_net*L*|grad(alpha)|` (interface-area-density form, Hardt & Wondra
2008). `maxSource` now caps the volumetric source directly.

**Affects:** evaporative cooling magnitude and location, `dgdt`
(alpha/pressure source), `phaseChangeMassFlux` (recoil input) — recoil is now
confined to the free surface as physics requires.

## 4. Recoil closure calibration + deterministic damping (`advancedInterfaceCapturing`)

**Was:** `p_recoil = ((2-beta_m)/(2*alpha_e))*j_net*sqrt(2*pi*R*T)` ≈
0.91·p_sat — ~70 % above the accepted Knudsen-layer result. The temporal
relaxation reference updated on *every* call, so damping depended on how many
times recoil was recomputed per step (outer iters × correctors × sub-cycles).

**Now:** `p_recoil = (recoilCoefficient/alpha_e)*j_net*sqrt(2*pi*R*T)` =
`recoilCoefficient*p_sat`, with `recoilCoefficient` a dictionary entry
(default 0.54, Anisimov 1968 / Knight 1979). The relaxation/rate-limit
reference is frozen once per time step (`relaxTimeIndex_`).

**Affects:** peak recoil (×0.59 vs before at equal T), jet velocity scale
(∝ sqrt(p), ≈ −23 %), and reproducibility w.r.t. PIMPLE iteration counts.

## 5. Lattice conductivity and Te/Tl advection (`twoTemperatureModel.C/.H`)

**Was:** `kl = Cl*De` used the *electron* diffusivity for the lattice →
~250 W/m/K vs a physical Ti phonon conductivity of ~5 W/m/K (melt zone smeared,
cooling far too fast). Te/Tl had no advection at all: ejected material left
its temperature behind (film moves a thickness-comparable distance during
ejection).

**Now:** optional dimensioned `kl` entry in `twoTemperatureProperties` (cases
set 4.5 W/m/K for Ti; falls back to Cl*De with a warning), and advective-form
transport `C*(div(phi,T*) − T*·div(phi))` for both temperatures (switch
`advectTemperatures`, default true; uses the existing `div(phi,Te)/div(phi,Tl)`
scheme entries).

**Affects:** melt-pool width, cooling curves, droplet thermal history,
solidification timing.

## 6. Bounded alpha source terms (`alphaSuSp.H`)

**Was:** `Sp` stored *positive* ("strengthens matrix" — it does the opposite
in the MULES convention) plus a spurious `2*Sp*alpha1` contribution.

**Now:** upstream Weller form: evaporation → `Sp -= dgdt/max(alpha1,1e-4)`
(bounded sink), condensation → bounded `Su/Sp` pair. Exact `−dgdt` net rate,
shuts off smoothly at alpha = 0 and 1.

## 7. Predictor/corrector force consistency (`UEqn.H`, `pEqn.H`)

* Explosive and plasma pressures now also enter `phig` in `pEqn.H`; before,
  they existed only in the momentum predictor and the corrector cancelled
  most of their effect.
* `recoilPressurePtr->relax(0.5)` removed — it was a **compile error**
  (non-const call through a `const` pointer), proving the pushed tree did not
  match the tested binary. Recoil smoothing lives in
  `advancedInterfaceCapturing::recoilRelax`.
* `UEqn.relax(URelaxationFactor)` is skipped on the final PIMPLE iteration so
  the converged state is unbiased.

## 8. Object lifetime and run-loop hygiene (`createFields.H`, `compInterFoam.C`)

* `autoPtr::ptr()` → `get()` (twice): `ptr()` released ownership, leaking the
  legacy recoil field and disabling its refresh/write path permanently.
* Removed the manual `functionObjects().start()/execute()` calls (2–3×
  duplicated sampling/probe output per step); `Time` handles them.
* Alpha sub-cycling decoupled from the electron micro-step count
  (`compressibleAlphaEqnSubCycle.H`) — they are unrelated stability controls.

## 9. Diagnostics corrected (`compInterFoam.C`)

* Electron energy audit uses `∫Ce dT` (`electronInternalEnergy`, now public)
  instead of `Ce(Te)*Te` (factor-2 error for the linear Ce model).
* Energy-balance warning compares the per-step change against the per-step
  boundary loss (was: cumulative total — warning became meaningless noise).
* Tracker kinetic energy uses metal *mass* (alpha·rho·V, was volume), and the
  laser→kinetic efficiency uses time-integrated absorbed energy (was
  instantaneous power × elapsed time).

## 10. Enhanced-physics integrity (`enhancedLIFTPhysics`)

* Breakup melting gate reads `mixture.T_melt()` (a **gold** value, 1337 K, was
  hardcoded in this titanium solver).
* The non-conservative direct alpha redistribution in `applyBreakup()` is now
  opt-in (`applyAlphaRedistribution`, default **off**); breakup fields remain
  as diagnostics. It created/destroyed metal volume with no source term.
* `explosiveMassSource` remains diagnostic-only by design: the phase-explosion
  feedback path is the `explosionIndicator` multiplier on the Hertz–Knudsen
  flux (coupling both would double-count).

## 11. Case configuration (all four cases)

| Change | Reason |
|---|---|
| `evaporationCoeff 0.18` | literature value for Ti (own diagnosis doc, Priority 2 — was 0.05/0.3) |
| `activationTime` → `tStart/tEnd` | `activationTime` was never read by the solver |
| `kl 4.5 W/m/K` added | Ti phonon conductivity (k_total − Wiedemann–Franz electronic part) |
| `advectTemperatures true` | makes the new Te/Tl transport explicit |
| `recoilCoefficient 0.54` | documents the new recoil calibration |
| `pRelaxationFactor` removed (fvSolution) | key no longer used (see fix 1) |
| `div(phi,T/Te/Tl)` → `Gauss limitedLinear 1` (TEST1/TEST3) | bare `Gauss linearUpwind` lacks its mandatory gradient argument and now these keys are actually exercised |

---

## What to expect on the first re-run (validation checklist)

1. **Build first** (`wmake` with OpenFOAM v2406 sourced) — this tree previously
   did not compile (fix 7), so verify the binary is rebuilt from it.
2. GAMG on `p_rgh` should converge in O(10) iterations; no
   `pressureClamp: enforcing bounds` messages in normal operation.
3. `Max |j_net|` of O(10²–10³) kg/m²/s at peak heating and recoil of
   O(10–10²) MPa (0.54·p_sat at the surface temperature), *only* in interface
   cells.
4. Film velocity of O(10²) m/s consistent with sqrt(2·p_recoil/rho_metal);
   the "velocity is N× the recoil-driven expectation" warning should stay
   near N ≈ 1.
5. The vapour plume now heats up (check `max(T)` in gas cells > 300 K during
   and after the pulse) and expands.
6. Total metal volume is conserved to the evaporated fraction (breakup no
   longer modifies alpha directly).
7. Energy-audit warnings should be rare; if the clamp-correction energy
   (`Temperature clamp energy correction`) is a significant fraction of the
   pulse energy, raise `maxTe/maxTl`.

## Known remaining modelling limitations (unchanged, by design)

* Metal EOS is Boussinesq (no acoustic/stress-confinement wave inside the
  film); gas-phase compressibility is now included, film-side is not.
* Condensation is inactive until a `pMetalVapor` field is provided; metal
  vapour inherits "air" gas properties.
* Fixed optical properties (reflectivity, absorption depth); no
  Te-dependence, no ballistic-electron correction.
* No latent heat of melting, constant surface tension (no Marangoni).
* The receiver substrate is initialized as liquid metal (no solid mechanics).
