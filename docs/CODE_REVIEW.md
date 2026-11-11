# Code review of compInterFoam-fsLIFT

Reviewed: commit `5844a69` (every file in `src/`, `tools/`, `cases/` and `docs/`). Each finding was checked against the code and, for the MULES item, against the OpenFOAM v2406 sources. The corrected solver was compiled with OpenFOAM v2406 and tested (§8).

**Status:** ✅ fixed · ❌ not fixed (reason given) · ℹ️ the review's claim was corrected.

**Effect on the old run:** what the issue did to the reference results in `figures/` and the docs, which were produced before the fixes. The reference case is being rerun with the corrected solver.

---

## 1. Defects that change results

| # | Issue | Status | Reason and fix | Effect on the old run |
|---|---|---|---|---|
| 1 | The lattice equation had an extra term C_l(T_l^old − T_l)/τ in every evaporating cell, with τ = 1e-12 s. The case key `relaxationTime` was not the one read (`evapRelaxationTime`) | ✅ | It multiplied the heat capacity of those cells by 1 + Δt/τ: 2× at Δt = 1 ps, 3.5× at 2.5 ps, and the energy audit did not count it. Removed; the evaporation sink is now explicit. That is stable: its stiffness time is ≥ 33 ps at the 10 kK cap, against Δt ≤ 2.5 ps | **Active.** The surface cells that set the vapour pressure heated and cooled 2–3.5× too slowly |
| 2 | The trapped-vapour pocket received only 25–50 % of the mass that evaporation removed from the film | ✅ | The film sink uses the geometric interface area. The pocket used j·\|∇α\|·V with a Gauss gradient, which gives half that area at a sharp step and less at a smeared one. The pocket inflow is now ρ_m·dgdt·V, exactly the mass the α equation removes | **Active.** The ideal-gas pocket pressure was too low wherever it, not p_sat, set the pressure |
| 3 | `MULES::correct` was given the full source, so evaporation was applied again in every correction pass | ✅ | Stock v2406 `VoF/alphaEqn.H` passes −Sp·α here. Now the same | **Not active**: item 4 kept the correction off |
| 4 | `MULESCorr` was switched off at the first femtosecond step and stayed off for the whole run, although `fvSolution` says `true`; a restart could switch it on again | ✅ | `alphaEqn.H` overwrote the run-long variable. It now uses a per-step copy (the same for `icAlpha`, `scAlpha`) | **Active.** The old run used plain explicit MULES throughout. The corrected run uses the MULES correction whenever Δt > 1 ps (the flight stage) |
| 5 | The laser energy was normalised over ∫Q dV, but the two-temperature model deposits α·Q | ✅ | Partly filled cells absorbed less than the pulse energy. Now normalised over ∫α·Q dV | Small while the film was sharp: the deposited energy was within 0.2 % |
| 6 | Recoil started at 3560 K, but evaporation removes mass and latent heat from 1941 K | ✅ | Recoil now acts in the same cells as evaporation | Negligible: below 3560 K the recoil is ≤ 0.055 MPa |
| 7 | Two specific gas constants: 174 in the evaporation flux, 173.6 (k_B/m) in the recoil and the pocket | ✅ | One value (174 J/(kg K)) everywhere | Recoil 0.1 % low, pocket pressure 0.2 % low |

---

## 2. Other logic problems

| Issue | Status | Reason and fix | Effect on the old run |
|---|---|---|---|
| The compression-limiter override came after the limiter had been applied | ✅ | It did nothing. Removed; the effective compression coefficient is min(c_α, λ) = 0.5, as documented | None |
| Gas velocity clipped to 2000 m/s by default in every cell with α < 0.1, after the flux was built | ✅ | It does not conserve momentum and made U and φ disagree. Now off by default and applied before the flux | None: gas speeds stayed near 200 m/s |
| `pressureClamp false` still applied an "emergency clamp" at ±7.5 GPa; the default was `true` | ✅ | "Off" now means off: p is never modified, a non-finite p stops the run. The default is `false` | Very probably none: the bounds are 14× the largest pressure (524 MPa) |
| `alphaMin`, temperature offsets and smoothing were silently overwritten for femtosecond cases | ✅ | The femtosecond defaults now apply only to keys missing from the dictionary, and are reported | Active: `alphaMin` 0.01 was replaced by 0 |
| The global energy check left out the laser input and the latent heat | ✅ | Its warning fired on every heating step. Both sources are now included | None (diagnostic) |
| The time step was limited only by `maxCo` and `maxDeltaT`, although the docs listed `maxAlphaCo` and `maxThermalCourant` | ✅ | `maxAlphaCo` now limits Δt over interface cells (stock VoF form). `maxThermalCourant` bounds only the gas–metal exchange; the docs now say so | None: with maxAlphaCo 0.5 > maxCo 0.2 the new limit does not bind |
| `gasMetalExchangeCoeff` Function1 looked up inside its own sub-dictionary | ✅ | Now looked up in the parent dictionary, like `Ce` and `G` | None: no case uses that form |
| Phase-explosion and plasma sub-models defaulted to **on** when their switch was missing | ✅ | Default is now off (both are unverified) | None: the cases switch them off |
| Saha prefactor 2^1.5 too large; critical density 1.1e27 m⁻³ (the value for 1 µm light) | ✅ | Saha: 2(g_i/g_0)λ⁻³. Critical density at 343 nm: 1.115e27/0.343² = 9.5e27 m⁻³ | None: plasma model off |
| Different defaults for the same key: `maxTe` (4000 / 3500 / 20000), `metalFractionFloor`, `metalAmbientBlendWidth`, `clampRecoil` | ✅ | One default each | None: the cases set these keys |
| ℹ️ The review said the Kapitza gas–metal exchange is "effectively off" | corrected | It is off at sharp interfaces, in cells with α < 0.5 and above 3560 K, but **active in partly filled interface cells (0.5 ≤ α < 1) below 3560 K** | Active in those cells |

---

## 3. Messages and comments that were wrong

| Issue | Status | Fix |
|---|---|---|
| Pulse-width check printed "10 fs" but tested 0.1 fs, and only warns | ✅ | Message corrected; `laserProperties` comments now say the checks only warn |
| `laserProperties` said a peakIntensity mismatch aborts the run | ✅ | The code rewrites `pulseEnergy` and warns; the comment now says so |
| Debug print every 10 steps inside a per-cell loop, claiming "~90 kg/m²/s at 5000 K" (the inputs give 537) | ✅ | Removed |
| Stale comments: missing dgdt feedback, "SLM: argon … powder" headers, references to files not in the repository, "the predictor uses grad(p_rgh + recoil)", "stock solver double-counts surface tension", "opt-in" next to "default on", "maxSource set to zero" | ✅ | Corrected |
| Garbled characters in two messages; π as 3.14159; k_B as 1.38e-23 | ✅ | Corrected |

---

## 4. Documentation and reported numbers

| Issue | Status | Fix |
|---|---|---|
| "Metal volume changed by −0.0135 %" included the 8 µm metal receiver | ✅ | For the film alone it is about **−1.5 %** (0.0135 / 0.885 % film share). README, RESULTS, MODEL and VALIDATION corrected |
| MODEL.md: electron diffusion across the film 0.05 ps | ✅ | (71.4 nm)² / 10⁻⁴ m²/s = **51 ps** |
| MODEL.md: absorbed fraction 0.5717; electron sub-step "at most 1 fs" | ✅ | 0.5713; the sub-step is ≤ 1 fs only in stage 1, then Δt/20 (50–125 fs) |
| VERIFICATION.md: the old absorption coefficient 1.07e7 1/m "from 2k/λ" | ✅ | 2k/λ = 1.65e7; the origin of 1.07e7 is unknown, now stated |
| Overclaims: "the code is verified", "converging" (two meshes), the Laplace error "settles" (it does not decrease: −2.75, −3.80, −4.09 %), a 14 % gap called "consistent", "four inputs" against "three" | ✅ | Reworded; the untested paths (Drude k_e, gas–metal exchange, pocket mass balance) are named |
| 4515 kg/m³ called "solid titanium at room temperature" (that is 4506); "the liquid is 7.7 % less dense" (7.1 %; 7.7 % is the expansion) | ✅ | Corrected |
| The mesh-study case has different inputs and a different radial mesh, so it is not an error bar for the reference result | ✅ | USER_GUIDE and VALIDATION say so |
| The mass counted in the detached droplet grows by 3.3 % between 13.4 and 19.8 ns | ✅ | Now stated in RESULTS as a sign of α smearing |

---

## 5. Case inputs

| Issue | Status | Fix |
|---|---|---|
| Keys the solver never reads: 15 in `fvSolution`, `relaxationTime` and `logRecoilSuppression` in `controlDict`, `temporalShape` and `spatialShape` in `laserProperties`, `Tsol`, `Tliq`, `Tvap`, `hf`, `kappa` in `thermophysicalProperties.metal`, `sigma`, `icAlpha`, `scAlpha` in `transportProperties` | ✅ | Removed from all four cases (each checked against `src/`) |
| Unread `metal{}` and `air{}` blocks in `thermophysicalProperties`; Ti density 4500 in two places | ✅ | Removed; one density per case |
| Unused files: `dynamicMeshDict`, `refineMeshDict`, two inactive `fvOptions` entries, an empty topoSet `focusBox` | ✅ | Removed |
| Gas density 1.2 in `transportProperties` against 1.1614 in the thermo file | ✅ | 1.1614 (the value only enters a switched-off term) |
| Kapitza model used argon values for a gas called air (383 / 1.7 = 225 m/s; argon's sound speed is about 322 m/s) | ✅ | Air: 1.1614 kg/m³ × 347 m/s |
| ρ·C_p = 2.93e6 against C_l = 2.5e6, and three metal conductivities | ❌ | The metal C_p and conductivity enter only the gas-temperature equation, which does not act back on the flow or the two-temperature model. Documented |
| `donor` declared `wall` in `blockMeshDict` but open in `0/` | ❌ | No effect with laminar flow; changing it would change the reference mesh |

---

## 6. Duplicated and dead code

| Issue | Status |
|---|---|
| 13 copy-pasted mass-weighted property functions; γ and W mass-weighted, which is physically wrong | ✅ Three helpers; γ = C_p/C_v, C_p/C_pv and 1/W = ΣY_i/W_i of the mixture |
| Two ~90-line α sub-cycle branches differing by one call | ✅ One code path |
| Repeated blocks: laser α test, pulse centre, pulse-train windows, recoil clamp, exchange cut-offs, vaporisation-temperature aliases, `simulationType` parsing | ✅ Each in one place |
| Dead code: unreachable LTS block, unused recoil clamp, no-op δ round trip, discarded C_e update, `explosionFluxMultiplier`, commented-out dynamic-mesh code and its libraries | ✅ Removed |
| Registry lookups inside per-cell loops; a message printed on every call | ✅ Hoisted; printed once |
| Process tracker: hard-coded 1941/3560 K, "PLASMA" and "SOLIDIFICATION" phases the model cannot produce, plume detected from T_l (held at ambient in gas), a "turbulent KE" that was not a fluctuation energy, `maxSource` labelled kg/(m³ s), header re-written after every restart | ✅ Corrected (column now `fluctuation_KE_nJ`) |
| Analysis tools hard-coded ρ = 4515 and σ = 1.64 | ✅ Read from the case |
| `postProcess.H` is the openfoam.org 2016–2018 version | ❌ `compInterFoam.C` relies on a macro this copy defines; it affects only the `-postProcess` option |
| `incidenceAngle` was read but unused; a Fresnel fallback with hard-coded constants (n = 2.0, k = 2.8, light from air) was used when `reflectivity` was missing | ✅ Both removed. `reflectivity` (or `transmission`) is now required; every case sets it, so no number changes. **Not yet compiled** |

---

## 7. Not fixed

| Issue | Reason |
|---|---|
| Gas–metal exchange energy: the gas gains the exchange heat per full cell volume, the lattice loses α times it | The exchange is active in partly filled cells, so a fix changes results, and the reference rerun was started with this code. Open; the energy involved is not yet quantified |
| `postProcess.H`, `donor` patch type, metal C_p and conductivity | No effect on any run (reasons in §5 and §6) |
| The `if (LTS)` branch of the time loop in `compInterFoam.C` | Unreachable: `createFields.H` stops the run when `localEuler` is selected. Harmless; left because removing it also touches `alphaEqn.H` and needs a compile |

---

## 8. Tests of the corrected solver (output)

| Test | Output |
|---|---|
| Compilation, OpenFOAM v2406 (`wclean; wmake`) | no errors, no warnings |
| Galilean invariance: film at −170 m/s, 0–100 ps | velocity exactly −170 m/s; pressure within 5 × 10⁻⁶ of 10⁵ Pa |
| 1-D laser and two-temperature column against the independent finite-difference code | peak electron temperature 6640 K against 6672 K (−0.5 %); electrons and lattice meet at 0.500 ps against 0.496 ps (+0.8 %); unchanged by the fixes |
| Clean-up (§3, §6) against the binary used for the rerun | all written fields (U, p, p_rgh, α, T_e, T_l, T) bit-identical on the Galilean test and the 1-D column: the clean-up changes no number |
| Static drop, Laplace pressure jump (MULES correction now active) | −2.75 % at R/Δx = 10, −3.80 % at R/Δx = 20: the same as before the fixes; drop volume within 1×10⁻⁴ % |
| Oscillating drop against Lamb, R/Δx = 10 | frequency −12.9 %, damping −52 %: the same as before (the weak damping is a known open point); volume within 3×10⁻³ % |
| Viscous decay of a velocity field | decay rate 0.9988 of the exact value at 5 ns (0.12 %), the same as before |
| Force on the film in the corrected reference run | change of momentum against the recoil and vapour impulse: 1.7 % (2–102 ps) and 4 % (102–202 ps); before the fixes the worst interval was 6 % |
| Vapour-pocket mass balance in the corrected reference run (fix 2) | pocket mass / metal lost by the film = 0.98–1.02 from 2 to 16 ps; afterwards lower, as designed (vapour above p_sat condenses; the gap side evaporates into the open gas). Before the fix the pocket received about half |
| Corrected reference run against the old one, launch | trapped-vapour pressure on the axis 961 MPa at 26 ps (old: 524 MPa, the factor ≈ 2 expected from fix 2); film speed on the axis about 10–30 % higher (205 against 185 m/s at 0.45–0.65 ns); no NaN, no fatal error, no energy-balance warning |
| Reference case rerun to 20 ns | running; results replace those in `figures/` and the docs when it finishes |
