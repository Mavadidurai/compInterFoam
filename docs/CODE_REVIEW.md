# Line-by-line review of compInterFoam-fsLIFT

Review date: 2026-11-07. Commit reviewed: `5844a69`.

**Status of the fixes (same branch).** A1, A2, A3 and A4 have been corrected in `src/`:
- A1: the relaxation pair was removed from `solveLatticeEquation`;
- A2: the pocket inflow now uses ρ_m·dgdt·V, the mass the α equation removes;
- A3: `MULES::correct` now gets `Su = −Sp·α`;
- A4: per-call copies of `MULESCorr`, `icAlpha` and `scAlpha`.

The medium-severity items have also been corrected:
- A5: the ineffective override was removed, so the code matches the documented effective compression of 0.5;
- A6: `maxVaporVelocity` is off by default and acts before the flux is built;
- A7: femtosecond defaults apply only to entries not in the dictionary and are reported, and the pocket's glass-side test no longer depends on the recoil α window;
- A8: `pressureClamp false` never modifies p (warning only; a non-finite p is fatal), and the default is now `false`;
- A9: the energy check includes the laser and latent-heat sources of the step;
- B2, B3: MODEL.md numbers, sub-step and time-step text corrected; density statements corrected;
- B4: USER_GUIDE and VALIDATION state that the mesh-study case is not an error bar for the reference result;
- C1: the unread `metal{}`/`air{}` blocks were removed from `constant/thermophysicalProperties`, `rho_metal` was set to 4515 in the reference and column cases, and the unread `defaultMetalDensity`/`defaultAirDensity` were removed.

A10 (Kapitza gradient, low–medium) is not changed: fixing it would switch on a gas–metal exchange that has never been active in any run.

These fixes were written without a compiler. They have since been compiled and tested, the remaining items fixed, and every finding checked against the code and the OpenFOAM v2406 sources: see **§10, verification of this review**. The reference results in `figures/` and the docs were produced with the old code; the reference case is being rerun with the fixed solver.

## 0. Scope, method and what this review could *not* do

| Part | Lines | Read |
|---|---|---|
| `src/` (C++, 24 files) | 16,693 | every line |
| `tools/` (Python and shell, 25 files) | 2,075 | every line |
| `cases/` (4 cases, 108 dictionaries) | ≈ 4,700 | the reference case in full; the other three as full `diff`s against it |
| `docs/`, `README.md`, results table | ≈ 2,000 | every line; numbers recomputed where possible |
| `literature/*.csv` | 120 rows | headers and counts only; the bibliographic entries were **not** checked against their DOIs |
| `figures/**/*.png` | – | not inspected (images) |
| `figures/**/data/*.csv` | – | used to recompute the numbers in the docs |

**Limits.** OpenFOAM is not installed in this environment, so nothing was compiled or run. Every finding below comes from reading the code plus small hand or Python calculations. Where a finding depends on OpenFOAM library internals that I am recalling rather than reading here, it is marked [UNCERTAIN].

Spatial layout used below (reference case, y = beam axis, r = radius, 5° wedge):

| Label | y range [µm] | Content | α (metal fraction) at t = 0 |
|---|---|---|---|
| L0 | 0 – 8.0 | receiver (metal) | 1 |
| L1 | 8.0 – 20.0 | gap (gas) | 0 |
| L2 | 20.0 – 20.0714 | Ti film | 1 |
| L3 | 20.0714 – 32.0714 | gas standing in for the glass | 0 |
| P_top | y = 32.07 | patch `donor` | open boundary (totalPressure) |
| P_r | r = 8 µm | patch `right` | open boundary |

---

## 1. Summary, most severe first

| # | Severity | Finding | Where |
|---|---|---|---|
| A1 | **High** | Artificial time-step-dependent damping of the lattice temperature in every evaporating cell | `twoTemperatureModel.C:784,790`, `twoPhaseMixtureThermo.C:1473-1477` |
| A2 | **High** | The trapped-vapour pocket receives only about 25–50 % of the mass that evaporation removes from the film | `advancedInterfaceCapturing.C:1152-1155,1248-1252` vs `twoPhaseMixtureThermo.C:1243-1263` |
| A3 | **High** (when active) | `MULES::correct` is given the full source `Su`, so phase change and compressibility are applied again in every correction pass | `alphaEqn.H:353-363,441-451` |
| A4 | **High** | `MULESCorr` is switched off permanently after the first femtosecond step; on a restart with Δt > 1 ps it comes back on, so the algorithm differs between runs and stages | `alphaEqn.H:393-405` |
| B1 | **High** (reporting) | The headline "metal volume changed by only −0.0135 %" includes the 8 µm metal receiver; for the film alone it is about **−1.5 %** | `compInterFoam.C:523-529`; README:130,280; MODEL:231; RESULTS:67; VALIDATION:112 |
| A5 | Medium | The femtosecond override of the compression limiter comes after the limiter has been applied, so the override does nothing | `alphaEqn.H:100-126` |
| A6 | Medium | Every cell with α < 0.1 (air included) is clipped to 2000 m/s by default, after φ has been computed | `UEqn.H:461-516` |
| A7 | Medium | `femtosecondCase_` silently overrides the user's `alphaMin`, temperature offsets and smoothing | `advancedInterfaceCapturing.C:130-137,434-451,480-489` |
| A8 | Medium | `pressureClamp false` does not switch clamping off: an "emergency clamp" still acts | `pEqn.H:807-859` |
| A9 | Medium | The global energy check leaves out the laser input, and the latent term is always 0 in this case, so its warning means nothing | `compInterFoam.C:1854-1970` |
| B2 | Medium | MODEL.md: electron diffusion time across the film is **51 ps, not 0.05 ps** (off by 1000×) | `docs/MODEL.md:180` (and the argument at :126) |
| B3 | Medium | MODEL.md says Δt is also limited by `maxAlphaCo` and `maxThermalCourant`; the code limits it only by `maxCo` and `maxDeltaT` | `docs/MODEL.md:449-452`, `compInterFoam.C:1662-1687` |
| B4 | Medium | The "mesh study" case has different inputs **and** a different radial mesh from the reference case, so it cannot give the reference result an error bar | `cases/fsLIFT_332mJ_Ti_meshStudy`, USER_GUIDE §2.4 |
| C1 | Medium | The reference case uses three titanium densities (4515 / 4500 / 4506) and has contradictory material blocks that are never read | see §4 |
| A10 | Low–Medium | The Kapitza gas–metal exchange (switched on in the case) uses a cell-limited ∇α, which is ≈ 0 at a sharp interface, so the exchange is effectively off | `twoTemperatureModel.C:2040` |

Short summary of §1: four code defects (A1–A4) can change the physical result, and one reported number (B1) is misleading. Most of the remaining items are wrong documentation, inconsistent inputs, or dead and duplicated code.

---

## 2. Defects in the code (A)

### A1. Artificial relaxation term in the lattice equation

`twoTemperatureModel.C:784` puts `Sp(α·Cl·relax, Tl)` on the left-hand side, and `:790` puts `α·Cl·relax·TlOld` on the right-hand side. Together they add

$$ C_l\,\frac{T_l^{old} - T_l}{\tau}, \qquad \tau = \max(\texttt{evapRelaxationTime},\ \texttt{dtFloor}) $$

to every cell where phase change is active (`twoPhaseMixtureThermo.C:1473-1477`).

- The case sets `relaxationTime 5e-11`, but the solver reads `evapRelaxationTime` (`twoPhaseMixtureThermo.C:682`). The default τ = 1e-12 s is therefore used silently. The case comment ("read and logged only") admits the key is unused, but does not say which value takes its place.
- Over one implicit step this multiplies the heat capacity by (1 + Δt/τ), as a calculation of mine:
  - Δt = 1e-14 (stage 1): 1.01, so negligible;
  - Δt = 1 ps (stage 2): 2.0, so the surface temperature responds half as fast;
  - Δt = 2.5 ps (stage 3): 3.5.
- It slows both evaporative cooling and conductive heating of exactly the surface cells that set p_sat and therefore the vapour drive.
- The result depends on the time step, and the energy audit does not count this term (`phaseChangeEnergy` includes only `Cl·phaseChangeSource`).
- Size of the effect on the 2-D result: [UNCERTAIN]. It depends on how often stage 2 actually reaches Δt ≈ 1 ps rather than a smaller Courant-limited step.

### A2. The vapour pocket gets only part of the evaporated mass

- **Mass removed from the film.** `twoPhaseMixtureThermo.C:1243-1263` uses the *geometric* area of the face that crosses α = 0.5. It assigns that area to the metal-side cell, so the removed mass is j·A.
- **Mass added to the pocket.** `advancedInterfaceCapturing.C:1248-1252` adds j·|∇α|·V, with ∇α from `fvc::grad(alpha, "grad(alphaConfinedVapour)")`. That gradient falls back to the default `Gauss linear` scheme, and j is non-zero only in the evaporating (metal-side) cell.
- A 1-D replica of both formulas gives a ratio of pocket mass to removed mass of:

| Interface profile (cell α values) | Ratio |
|---|---|
| sharp: 1, 1, 1, 0, 0 | **0.50** |
| one smeared cell: 1, 1, 0.6, 0, 0 | 0.50 |
| smeared over two cells | 0.40 |
| smeared over three cells | 0.25 |

- Consequence: p_ideal = m·R·T/V is under-predicted by 2× or more wherever the ideal-gas branch of p_b = min(p_sat, p_ideal) is the binding one. Whether that branch binds at the 524 MPa peak cannot be told from the shipped data [UNCERTAIN].
- This contradicts MODEL.md §6.2, which says the same |∇α| is used for the film and the pocket. The force-balance test checks momentum, not pocket mass, so it cannot catch this.

### A3. Source applied twice in `MULES::correct`

- `alphaEqn.H:353-363` and `:441-451` pass `Sp(), Su()` to `MULES::correct`.
- `MULES::correct` updates ψ = (ψ·rΔt + Su − ∇·φ_corr)/(rΔt − Sp). The sources are therefore applied again in every correction pass, on top of the predictor, which already contains them.
- Stock `compressibleInterFoam` passes `Su = −Sp·α` at this point for exactly that reason. [UNCERTAIN: this is from memory of the OpenFOAM `MULESTemplates.C` and the stock `alphaEqn.H`; please check against the v2406 source.]
- The defect acts only when `MULESCorr` is on (see A4).

### A4. `MULESCorr` changes with the run history

- `alphaEqn.H:403` writes `MULESCorr = false` into the variable created in `createFields.H:840`. That variable lives for the whole run.
- The first step has Δt = 1e-15, so MULESCorr is off for the rest of an uninterrupted run, although `fvSolution` says `MULESCorr true`.
- Each new solver process (`tools/run_case.sh` stages, `RESUME=1`) re-reads `true`. Whether it then switches off again depends on whether the restart Δt is ≤ 1e-12.
- Example from my calculation, not from a run: in `tools/run_mesh_study.sh`, level 1 has `MAXDT2` = 1.68e-12. Its stage 3 (1 → 20 ns) can therefore run with MULESCorr **on**, including defect A3, while levels 2 and 3 run with it off. That would make the 20 ns comparison between levels algorithm-dependent. The comparison at 0.802 ns is not affected.
- `icAlpha` and `scAlpha` are overwritten the same way (`alphaEqn.H:116-117`); they are 0 in the case anyway.

### A5. The compression-limiter override is in the wrong place

- `alphaEqn.H:105` applies `phic = min(cAlpha·|u|, 0.5·|u|)` first.
- Only after that does `:113-118` set `compressionLimiter = 1` for femtosecond time steps, which has no effect any more.
- The effective coefficient is therefore always 0.5. The docs state this, but the code clearly intends something else.

### A6. Default velocity clip on gas cells

- `UEqn.H:461-516` (`maxVaporVelocity`, default 2000 m/s, not set in the case) rescales U in every cell with α < 0.1, including plain air.
- This happens after `phi = fvc::flux(U)` (line 457), so U and φ disagree until the pressure equation.
- It does not conserve momentum. It is probably rarely active, because the case's gas speeds are about 200 m/s.

### A7. Silent override in `advancedInterfaceCapturing`

`femtosecondCase_` defaults to true whenever the initial Δt ≤ 5 ps, which includes this case. It then:
- overwrites `alphaMin` (case value 0.01 → 0), `alphaMax`, `phaseChangeTempOffset` and `recoilTempOffset`;
- turns ramping, smoothing and rate limiting off;
- does all this without printing anything.

### A8. `pressureClamp false` still clamps

With `pressureClamp false`, `pEqn.H:807-859` still clamps p to fallback bounds and prints a warning.
- In this case the bounds work out to ±7.5 GPa, widened by 1.5× the observed recoil pressure, so the bounds move during the run.
- They are wide enough that they probably never act, but "false" does not mean "off".
- The default when the key is missing is `true` (`pEqn.H:429,648`), which is the behaviour the case comment says broke continuity.

### A9. The global energy check does not work

- `compInterFoam.C:1854-1970` compares E_tot between steps but leaves out the laser input.
- `El = cumulativePhaseChangeEnergy` is accumulated only when `energyDiagnostics` is on (`twoTemperatureModel.C:1004-1019`), and it is off in the case, so El is always 0.
- During heating the "energy change exceeds tolerance" warning must therefore fire on every step, so it tells the user nothing.

### A10. Other code defects (lower impact)

- **Kapitza gas–metal exchange.** `twoTemperatureModel.C:2040` uses `fvc::grad(alpha.metal)`, which is `cellLimited Gauss linear 1` in `fvSchemes`. This is the scheme the code's own comment (`twoPhaseMixtureThermo.C:1208-1213`) measured as 3.7e-17 at a sharp step. The coefficient is then 0 (`:2078`). In addition, `minMetalFractionForExchange 0.5` and `gasMetalCutoffTemperature 3560` disable it. So it is configured but in practice off.
- **Function1 lookup.** `twoTemperatureModel.C:425-431` calls `Function1::New("gasMetalExchangeCoeff", gasMetalDict)`, which looks the key up inside its own sub-dictionary. `Ce` and `G` correctly pass `dict_`. This only matters for a non-Kapitza dictionary form, which no case uses.
- **Laser vs TTM weighting.** The laser model normalises ∫Q dV without α (`femtosecondLaserModel.C:2336-2362`), but the TTM applies α·Q (`twoTemperatureModel.C:864`). Partly filled cells therefore absorb less energy than intended. This is negligible while the film is sharp.
- **Gas–metal energy bookkeeping.** The gas gains exchange·(Tl − T) per full cell volume (`TEqn.H:1196-1197`). The lattice loses α·gasMetalHeatFlux one step later (`twoTemperatureModel.C:791`). The two differ by a factor α.
- **Recoil threshold.** Recoil acts only at T ≥ T_vap = 3560 K (`advancedInterfaceCapturing.C:780-789`), while evaporation removes mass and latent heat from T_melt = 1941 K. Between those temperatures there is mass loss without recoil. MODEL.md §5.3 does not mention this.
- **Saha factor (sub-model off).** `enhancedLIFTPhysics.C:445` uses (2/λ²)^1.5 = 2^1.5/λ³. Saha has 2·g_i/g_0/λ³, so the expression is off by about 2^1.5.
- **Plasma model (off).** The critical density 1.1e27 m⁻³ is the value for about 1 µm light; at 343 nm it is about 9.4e27.
- **Two explosion multipliers.** The mixture uses a hard-coded 1 + 99·ξ (`twoPhaseMixtureThermo.C:1417`). `explosionFluxMultiplier()` uses (M − 1)·ξ² and is never called (dead code).
- **Defaults for unverified sub-models.** Phase explosion, plasma and breakup default to **enabled** when their coefficient dictionary exists but the switch is missing.
- **Inconsistent defaults for the same key:**
  - `maxTe`: 4000 (`twoTemperatureModel.C:170`), 3500 (`:1378`) and 20000 (tracker);
  - `metalFractionFloor`: 1e-6 (TTM) and 1e-8 (`TEqn.H:221`);
  - `metalAmbientBlendWidth`: 1e-3 and 5e-3;
  - `clampRecoil`: false (AIC) and true (`UEqn.H:44`).

Short summary of §2: A1 and A2 act directly on the vapour drive, which is the main launch mechanism. A3 and A4 make the α transport depend on the run history. The remaining items are smaller or affect switched-off code.

---

## 3. Documentation and reported numbers (B)

### B1. Metal-volume "conservation" number

- The tracker sums α·V over the **whole** domain (`compInterFoam.C:523-529`). The receiver L0 is metal (α = 1).
- Film share of the metal (all layers span the same radius): 0.0714 / 8.0714 = 0.885 %.
- So −0.0135 % of the total corresponds to −0.0135 / 0.00885 ≈ **−1.5 % of the film**.
- The docs present it as a conservation check (VALIDATION.md:112 is titled "Conservation").

### B2. Errors in MODEL.md

- **:180.** h²/D_e with a cold lattice is (71.4e-9)² / 1e-4 = 5.1e-11 s = **51 ps**, not 0.05 ps. The statement at :126 ("electrons carry part of the energy deeper during the pulse") leans on that number; the actual diffusion length in 0.3 ps is √(D·t) ≈ 5.5 nm.
- **:122.** (1 − R)(1 − e^−ah) = 0.57168 × 0.99938 = **0.5713**, not 0.5717.
- **:192.** "Electron sub-step at most 1 fs" holds only in stage 1. `maxElectronSubCycles 20` caps the sub-step at Δt/20, which is 50 fs at 1 ps and 125 fs at 2.5 ps.
- **:449-452.** The time step is set by the finiteVolume `setDeltaT.H`, which uses `maxCo` and `maxDeltaT` only:
  - there is no local `setDeltaT.H` in `src/`;
  - `alphaCoNum` (`compInterFoam.C:1666-1678`) is computed over all cells, only logged, and never used;
  - `maxThermalCourant` only scales the gas–metal exchange coefficient (`TEqn.H:979-1054`).

  The claim in VALIDATION.md:111 ("below every stability limit … interface Courant, thermal limit") is wrong for the same reason.
- **§6.2.** Says the pocket inflow uses the same |∇α| as the film sink; it does not (A2).

### B3. Density and expansion statements

- MODEL.md:497, RESULTS.md:109 and VALIDATION.md:137 call 4515 kg/m³ "solid titanium at room temperature". That value is **4506** (the repository's own mesh-study case says so).
- "The liquid is 7.7 % less dense": 4193/4515 is **7.1 %** less dense. The 7.7 % is the volume expansion.

### B4. Explanation of the old absorption coefficient

VERIFICATION.md:123 says the old value 1.07e7 1/m came from 2k/λ. With the case's k = 2.8242 and λ = 343 nm:
- 2k/λ = 1.65e7;
- k/λ = 8.2e6;
- 2πk/λ = 5.2e7.

None of these is 1.07e7, so the explanation does not check out.

### B5. Laser-input comments that the code contradicts

| Comment | What the code does |
|---|---|
| `laserProperties:34-35`: checks "enforce pulseWidth 10 fs–10 ps … peakIntensity ≤ 6e16" | `femtosecondLaserModel.C:551-577` only **warns**, uses 1e-16 s (0.1 fs) as the lower bound, and warns on intensity above 1e17 |
| `laserProperties:86-87`: a peakIntensity mismatch "aborts" | `:285-310` overwrites `pulseEnergy` and prints a warning |
| VALIDATION.md:269: "the solver accepts pulse durations from 10 fs to 10 ps" | there is no hard limit, only a warning |

Also: `incidenceAngle` is read but never used, `temporalShape` and `spatialShape` are never read, and `checkPhysicalBounds` prints "10fs" while testing against 0.1 fs.

### B6. Overclaims

- README.md:78, VERIFICATION.md:15 and VALIDATION.md:35,251 say "every added equation verified; the code is verified". Not covered by any test:
  - the Drude k_e path (the column test uses k_e = C_e·D_e and an absorption coefficient 10× smaller than the reference);
  - A1, A2 and A3;
  - gas–metal exchange;
  - pocket mass conservation.
- VERIFICATION.md:137 itself admits that the logged force budget's "pressure" entry (60–75 %) is unexplained.
- "Frequency converging" is based on two meshes only.
- The Laplace-jump error grows under refinement (−2.75 → −3.80 → −4.09 %). Calling that "settles" is generous; that standard CSF behaves this way is [UNCERTAIN].
- VERIFICATION.md:41 calls 1.6 against an expected 1.4 (a 14 % gap) "consistent".
- VALIDATION.md says "failed for 3 inputs" (:31, :176) but "four corrections" (:35, :254, :267).

### B7. Stale comments

| File | Stale text |
|---|---|
| `thermophysicalProperties.air` | says the solver lacks the dgdt compressibility feedback; `dgdtComp` has since been added |
| `0.orig/alpha.*` | "SLM: Argon … metal powder" (copied from a selective-laser-melting case) |
| `twoTemperatureModel.C:1968` | refers to `RealisticLIFT/system/controlDict`, which is not in the repository |
| `pEqn.H:280` | refers to `CONVERGENCE_ISSUE_DIAGNOSIS.md`, which is not in the repository |
| `pEqn.H:26-28` | says the predictor uses grad(p_rgh + recoil); it does not |
| `UEqn.H:199-210` | says stock compressibleInterFoam "double-counts" surface tension and gravity; it does not, because HbyA = rAU·UEqn.H() does not contain the predictor's right-hand side (`pEqn.H:24`) |
| `alphaEqn.H:264` vs `:266` | "opt-in until re-verified" next to "default on" |
| `controlDict:71-73` | "set to zero so the cap is disabled", while the value is 1e25 |
| `twoPhaseMixtureThermo.C:1442` | debug line "Expected at 5000K: ~90 kg/m²/s"; the case inputs give **537** |

### B8. Docs numbers checked and found **correct**

All of the following match the shipped data or hand calculation:
- README and RESULTS: 191.7 m/s at 0.536 ns; the 0.54–10 ns table; boiling at 2.86 ns; 524.0 MPa at 26 ps; below 1 MPa from 0.586 ns; centre of mass 101.7 → 97.7 m/s; fit 100.06 m/s; 0.254 µm³ ≈ 0.79 µm; We = 7.2; KE/2σ = 1.81; drag estimates.
- MODEL §3 and §8: I₀, Q₀, R, absorption depth, E, 19,700 K, τ_eq, Re, Oh, Taylor–Culick speed, Bo, capillary Δt, the p_sat/j table.
- LITERATURE §3.3: all three sensitivity cases reproduce from their CSVs.

Not mentioned in the docs: the mass of the detached droplet *grows* by 3.3 % between 13.4 and 19.8 ns (`droplet_motion.csv`). A detached body should not gain mass; this is a sign of α-threshold or smearing effects.

Short summary of §3: the results numbers are reproducible, but several explanatory numbers in MODEL.md are wrong, and the verification and conservation claims are stronger than the evidence.

---

## 4. Case inputs (C)

- **Three titanium densities in one case:**
  - 4515: `thermophysicalProperties*`, `transportProperties`, `confinedVapourPocket.metalDensity`;
  - 4500: `controlDict:184` (Kapitza `rho_metal`) and `fvSolution:159` (`defaultMetalDensity`, never read);
  - 4506: the mesh-study case.
- **Contradictory, unread material blocks.** `constant/thermophysicalProperties` contains `metal{}` and `air{}` blocks (air as perfectGas R = 287; metal β = 7.6e-5). The solver reads `thermophysicalProperties.metal` and `.air` instead (air rhoConst, metal β = 0). [UNCERTAIN: this is the stock behaviour of `rhoThermo::New(mesh, phaseName)`.] The key scan agrees: `R`, `beta` and `T0` in that file are never referenced.
- **Unread keys in `thermophysicalProperties.metal`:** `Tsol`, `Tliq`, `Tvap`, `hf`, `kappa`.
- **Unused or overridden keys in `transportProperties`:** `icAlpha`, `scAlpha` and `phaseChangeMaxSourceDefault` are unused or overridden, and `sigma` is repeated in two files (the one in `thermophysicalProperties` is the one used).
- **Inconsistent lattice heat capacity:** ρ·Cp = 4515 × 650 = 2.93e6 against the TTM's Cl = 2.5e6 J/(m³ K).
- **Inconsistent metal conductivity:** μ·Cp/Pr = 47.7 W/(m K) (used in the gas temperature equation), `kappa 17.2` (unused), and kl 4.5 + ke0 18.9 in the TTM.
- **Kapitza "argon" values for a gas called air:** Z_gas/ρ_gas = 383/1.7 = 225 m/s, while argon's sound speed is about 320 m/s [UNCERTAIN, from memory].
- **Gas density** 1.2 (`transportProperties`) against 1.1614 (thermo).
- **Unused case files:** `donor` is declared `wall` in `blockMeshDict` but is an open boundary in `0/`; `dynamicMeshDict` and `refineMeshDict` are unused (the solver never builds a dynamic mesh); `fvOptions` entries are inactive; the topoSet `focusBox` is empty.
- **G table.** By the case's own comment it exceeds the published titanium values from 2000 K up (5–40× at 5–20 kK), and Ce = 630·Te is about 2× published. This is documented, but the reference result rests on it.
- **Mesh-study case.** Its level 2 has 35 nm radial cells for r < 3.5 µm against 50 nm in the reference, plus different G, Ce, ρ and σ. It is a new case, not a mesh study of the reference.
- **Unread keys that look meaningful:**
  - `fvSolution compInterFoamCoeffs`: `defaultMetalDensity`, `defaultAirDensity`, `maxInitialResidual`, `enableAlphaClamp`, `minAlpha`, `maxAlpha`, `densityClamp`, `minDensity`, `maxDensity`, `enableRAUClamp`, `minRAU`, `maxRAU`, `minRAUf`, `maxRAUf`, `targetResU`, `targetResP`;
  - `recoilSmoothingFactor` (its comment promises smoothing that does not exist);
  - `controlDict`: `relaxationTime`, `model`, `logRecoilSuppression`;
  - `laserProperties`: `temporalShape`, `spatialShape`.

---

## 5. Dead code and duplication (D)

### Duplicated code

| Duplication | Where |
|---|---|
| 13 copy-pasted mass-weighted property functions (one template would do); `gamma` and `W` are not mass-additive quantities | `twoPhaseMixtureThermo.C:1663-2326` |
| two ~90-line branches that differ by one call | `compressibleAlphaEqnSubCycle.H` |
| `inFilm` and `depositInMetal` computed twice | `femtosecondLaserModel.C:1616-1654` |
| pulse-centre logic | `:1310-1332` and `:2139-2163` |
| pulse-train window loop | `:1258` and `:2216` |
| clamp block | `advancedInterfaceCapturing.C:839-862` and `:865-888` |
| cutoff checks | `twoTemperatureModel.C:2068-2099` and `:2122-2149` |
| `updateElectronHeatCapacityField` duplicates `electronHeatCapacity()`, and its result is discarded | `twoTemperatureModel.C:1519-1537,1616` |
| no-op round trip: coefficient × δ then ÷ δ | `TEqn.H:433-438` |
| alias lookups, partly repeated | `twoPhaseMixtureThermo.C:163-348` |
| `simulationType` parsing | `createFields.H:603-681` and `compressibleInterPhaseTransportModel.C` |

### Dead code

- `createFields.H:25-62`: code after a FatalError, never reached.
- `UEqn.H:67-120`: the clamped recoil is computed and never used; the traction at `:234` uses the unclamped field.
- The dynamic-mesh blocks are commented out (`compInterFoam.C:1479-1492,1691-1729`), but `Make/options` still links `dynamicFvMesh` and `topoChangerFvMesh`.
- `try/catch(IOerror)` around `rhoThermo::New` (`twoPhaseMixtureThermo.C:893-915`): these exceptions are not thrown unless `throwExceptions` is set.
- `explosionFluxMultiplier()` is never called.
- The 1.6k-line ASCII-art `liftProcessTracker` (enabled in the case):
  - hard-codes the titanium temperatures 1941 and 3560 K instead of reading T_melt and T_vapor from the model;
  - labels a "PLASMA FORMATION" phase although the plasma model is off, and a "SOLIDIFICATION" phase although the model has no solid phase;
  - detects the vapour plume from Tl in gas cells, which the TTM pins to ambient;
  - uses a "turbulent KE" formula that is not a fluctuation energy;
  - re-appends the CSV header after every restart;
  - labels `maxSource` "kgm3s" though it is W/m³;
  - logs `relaxationTime` = 5e-11 although 1e-12 is used.
- `tools/benchmarks/drop_analyze.py:21`: a dead `if False` expression.

### Minor

- Hard-coded π as 3.14159 (`twoPhaseMixtureThermo.C:1401`) and k_B as 1.38e-23.
- Mojibake from a wrong text encoding: `TEqn.H:126` (`â‰ˆ`), `pEqn.H:958` (`Ã—`).
- A debug print every 10 steps inside a per-cell loop (`twoPhaseMixtureThermo.C:1427-1443`).
- `foundObject` called inside per-cell loops (`:1409`, `femtosecondLaserModel.C:1745`).
- An Info message printed on every call (`twoPhaseMixtureThermo.C:1538`).
- The density 4515 is hard-coded in `droplet_motion.py`, `metal_momentum.py`, `results_package.py` and `plot_droplet_motion.py`, so they give wrong masses for the 4506 mesh-study case. The laser end time is hard-coded in `compare_ttm_1d.py:44`.
- `postProcess.H` is the openfoam.org 2016–2018 version, not the v2406 one.

---

## 6. What I checked and found correct

- **Laser:** spatial and temporal normalisation, the erf window, the wedge sector coverage, exact cell-average Beer–Lambert deposition; all laser input arithmetic.
- **Phase change:** `alphaSuSp.H` Sp/Su algebra (net source = −rate in both branches); Hertz–Knudsen and Clausius–Clapeyron as implemented; the geometric interface-area construction and its processor-patch logic.
- **Recoil:** the 0.54·p_sat closure and the sign-aware localMax face value.
- **Flow:** pEqn dilatation terms, flux reconstruction and the stock HbyA form.
- **Finite-difference TTM** (`tools/ttm_fd_1d.py`): FWHM pulse, energy-conserving midpoint C_e, block-Thomas solver.
- **Grid-convergence script** (`mesh_study_gci.py`): Celik 2008 procedure, including the iteration for p, sign handling and GCI = 1.25·e_a/(r^p − 1).
- **Benchmarks:** Lamb ω₀² = 8σ/(ρR³), β = 5ν/R², the J₀ third zero.
- **Mesh generator** (`make_mesh_level.py`): cell counts 62,010 / 124,200 / 248,252 reproduce.
- **Files:** `LICENSE` is the unmodified GPL-3.0 text (md5 matches).

---

## 7. The opposing view (steelman)

- **"These are research-code rough edges, not errors."** Partly true. Most of D (duplication, dead keys, banners) does not change any number. A1–A4, however, sit on the main physics path: the vapour drive, the α transport and restarts.
- **"A1 is a deliberate semi-implicit linearisation of evaporative cooling."** That reading is plausible: S(T) ≈ S(T⁰) − (C_l/τ)(T − T⁰). But a linearisation should use the actual ∂S/∂T, which here is a steep exponential. It should not use a fixed 1/τ that also brakes heating and is not counted in the energy audit.
- **"A2 hardly matters, because p_b is capped by p_sat."** If the saturation branch binds throughout, the pocket mass deficit has no effect on the pressure. The shipped data do not record p_sat per bin, so this cannot be decided [UNCERTAIN]. It needs one logged run.
- **"The −0.0135 % is labelled 'total metal volume', so it is not wrong."** The label is literally correct. Presenting it as a conservation result for a 71 nm film next to an 8 µm metal receiver is still misleading.
- **"The docs are unusually candid."** True, and to the author's credit: many limitations are stated openly (no validation, separation limited by the mesh, inputs off the literature values, unexplained force budget). The problems above are mostly places where code and docs have drifted apart.

---

## 8. Suggested order of work

1. **A2.** Use one interface-area definition for the film sink and the pocket source (for example, publish the mixture's `interfaceAreaDensity` field and use it in the pocket). Add a test that the pocket mass equals the mass removed from the film.
2. **A1.** Remove the relaxation term, or base it on the real ∂S/∂T. Fix the `relaxationTime` / `evapRelaxationTime` key mismatch.
3. **A3 and A4.** Pass `−Sp·α` to `MULES::correct`. Never write the configuration variables `MULESCorr`, `icAlpha` and `scAlpha` inside `alphaEqn.H`; use local copies.
4. **B1.** Report the film volume (metal above y = 8.05 µm), as `mesh_study_fields.py` already does.
5. **B2 and B3.** Correct MODEL.md (51 ps, time-step limits, sub-step size, 4506 / 7.1 %).
6. **Rerun.** After 1–3, repeat the launch phase to 1 ns. Only then is a mesh study of the reference inputs meaningful.
7. **Clean-up.** Remove the dead keys and files, deduplicate the code in §5, and make the shared defaults consistent.

---

## 9. Confidence check

| Claim | Confidence | Why |
|---|---|---|
| A1 relaxation term exists and scales as (1 + Δt/τ), τ = 1e-12 | High | read directly in the code; key-mismatch confirmed by searching the source |
| A1 size of the effect on the 2-D result | Low | depends on the actual stage-2 Δt; no logs available |
| A2 pocket receives 25–50 % of the evaporated mass | Medium–High | both formulas read in the code and reproduced in 1-D; assumes the default `Gauss linear` for `grad(alphaConfinedVapour)` (confirmed in `fvSchemes`) |
| A2 effect on the peak pressure | Low | depends on whether p_sat or p_ideal binds |
| A3 double-counted source | Medium | depends on my recollection of `MULES::correct` and the stock pattern [UNCERTAIN] |
| A4 MULESCorr stays off after the first step | High | plain C++ variable scope |
| A4 mesh-study level 1 runs with MULESCorr on | Medium | requires the stage-2 Δt to actually reach > 1e-12 |
| B1 the −0.0135 % is diluted by the receiver | High | tracker code and `setFieldsDict` read; arithmetic simple |
| B2 51 ps, not 0.05 ps | High | arithmetic |
| B3 only `maxCo` limits Δt | Medium–High | no local `setDeltaT.H`; relies on the standard finiteVolume header |
| Docs results numbers reproduce | High | recomputed from the shipped CSVs |
| C1 unused contradictory `thermophysicalProperties` blocks | Medium–High | key scan plus stock behaviour of `rhoThermo::New` |
| Kapitza / argon sound-speed value | Low–Medium | from memory |
| Nothing was compiled or run | — | no OpenFOAM in this environment; all runtime statements are inferred |

---

## 10. Verification of this review (compiled and run, 2026-11-08)

The review above was written from the source alone. Every finding was then checked against the code at `5844a69`, and A3 against the OpenFOAM v2406 sources installed with the solver. The fixes were compiled (no warnings) and tested. Verdicts: **confirmed** (the finding is right), **corrected** (right, with a correction below), **wrong**.

### 10.1 Verdicts

| Item | Verdict | Evidence | Status |
|---|---|---|---|
| A1 relaxation term | confirmed | `Sp(α·C_l/τ, T_l)` with `α·C_l·T_l^old/τ` on the right; τ = max(1e-12, dtFloor 1e-12) | fixed (term removed). The evaporation sink is now explicit; it is stable because its stiffness time, 1/(∂S/∂T)/C_l, is ≥ 33 ps at the 10 kK cap against Δt ≤ 2.5 ps |
| A2 pocket mass | confirmed | Gauss-linear |∇α| in the metal-side cell of a sharp step is 0.5/Δx; the film sink uses the full face area | fixed (pocket inflow = ρ_m·dgdt·V, the mass the α equation removes) |
| A3 MULES source | **corrected** | not uncertain: v2406 `applications/solvers/multiphase/VoF/alphaEqn.H` passes `(-Sp*alpha1)()` to `MULES::correct`, and `CMULESTemplates.C` updates ψ = (ρψ/Δt + Su − ∇·φ_corr)/(ρ/Δt − Sp). **Correction:** A3 never acted in the reference run, because A4 kept `MULESCorr` off in all three stages | fixed |
| A4 MULESCorr | confirmed | each stage of `run_case.sh` restarts with Δt ≤ 1 ps, the femtosecond test switches `MULESCorr` off at the first step, and the variable stays off. So the old reference run used plain explicit MULES throughout | fixed. With the fix, the MULES correction is **on** whenever Δt > 1 ps, i.e. in the flight stage: the corrected run uses a different α algorithm there |
| A5 limiter override | confirmed | the override came after the limiter | fixed (removed) |
| A6 gas speed clip | confirmed | **addition:** it could not act in the reference run: the largest gas speed was about 200 m/s against the 2000 m/s limit | fixed (off by default; applied before φ) |
| A7 silent overrides | confirmed | | fixed (reported, only for missing keys) |
| A8 emergency clamp | confirmed | its bounds are ±1.5 × `maxRecoilPressure` = ±7.5 GPa; they widen only if the recoil exceeds 5 GPa, which it never did, so they did not move. They are more than 10× the largest pressure, 524 MPa, so the clamp very probably never acted [inferred, the logs do not record it] | fixed (off means off; non-finite p is fatal) |
| A9 energy check | confirmed | | fixed; compiles |
| A10 Kapitza exchange | **partly wrong** | the coefficient is zero where α < 0.5 (`minMetalFractionForExchange`), where the cell has no gas, where the cell-limited |∇α| is ~0 (a sharp step) and where T_l ≥ 3560 K. It is **not** zero in partly filled interface cells (0.5 ≤ α < 1) cooler than 3560 K, so it was active in those cells in every run, not "in practice off" | not changed: the gradient scheme decides where it acts; changing it changes results. The case comment now states where it acts |
| A10 Function1 lookup | confirmed | | fixed (looked up in the parent dictionary) |
| A10 laser vs TTM α | confirmed | | fixed: the laser normalises ∫αQ dV, as the TTM deposits αQ |
| A10 gas–metal bookkeeping | confirmed | the gas gains the exchange heat per full cell volume, the lattice loses α times it | **not changed**: since the exchange is active in some cells (row above), fixing it changes results, and the reference rerun was started with this code. Open item; the energy involved is not quantified |
| A10 recoil threshold | confirmed | | fixed: recoil from T_melt, where evaporation starts (below T_vapor the recoil is ≤ 0.055 MPa) |
| A10 Saha, critical density | confirmed | Saha: 2(g_i/g_0)λ⁻³; n_c = 1.115e27/(λ/µm)² = 9.5e27 m⁻³ at 343 nm | fixed (sub-model off) |
| A10 explosion multipliers | confirmed | | unused `explosionFluxMultiplier` removed |
| A10 defaults | confirmed | | phase explosion and plasma default to off; `maxTe`, `metalFractionFloor`, `metalAmbientBlendWidth`, `clampRecoil` defaults made consistent. The breakup model still defaults to on: it is a diagnostic and changes α only with `applyAlphaRedistribution` |
| B1 metal volume | confirmed | −0.0135 % of all metal ÷ 0.885 % film share = −1.5 % of the film | docs corrected |
| B2 MODEL numbers | confirmed | (71.4 nm)²/10⁻⁴ m²/s = 51 ps | corrected |
| B3 time-step limits | confirmed | the standard `setDeltaT.H` uses `maxCo` and `maxDeltaT` only | **code changed:** `maxAlphaCo` now limits Δt over interface cells (stock VoF form). With maxAlphaCo 0.5 > maxCo 0.2 it does not bind in these cases. `maxThermalCourant` still only bounds the exchange; the docs say so |
| B4 (§3) absorption coefficient | confirmed | 2k/λ = 1.65e7, k/λ = 8.2e6, 2πk/λ = 5.2e7 | corrected. Note: the review uses "B4" twice (summary table: mesh-study case; §3: this item) |
| B4 (summary) mesh-study case | confirmed | | docs say so |
| B5, B6, B7 | confirmed | | corrected |
| B8 droplet mass | confirmed | 1.5396e-17 → 1.5904e-17 kg (+3.3 %) between 13.4 and 19.8 ns | stated in RESULTS |
| C inputs | confirmed | key scan of `src/` for every listed key | removed from all four cases: the unread keys of `fvSolution`, `controlDict` (`relaxationTime`, `logRecoilSuppression`) and `laserProperties`; `Tsol`, `Tliq`, `Tvap`, `hf`, `kappa` from `thermophysicalProperties.metal` (read only from `controlDict`); `sigma`, `icAlpha`, `scAlpha` from `transportProperties` (read from `thermophysicalProperties` and `fvSolution`); `dynamicMeshDict`, `refineMeshDict`, the two inactive `fvOptions` entries and the empty topoSet `focusBox`. Gas density 1.1614 in `transportProperties`; Kapitza values for air instead of argon (argon: 383/1.7 = 225 m/s against a real 322 m/s). `phaseChangeMaxSourceDefault` **is** read, but `maxSource` overrides it. Kept and documented: `donor` as `wall` in `blockMeshDict` (no effect with laminar flow; changing it changes the reference mesh); ρ·C_p = 2.93e6 against C_l = 2.5e6 and μC_p/Pr = 47.7 W/(m K): the metal C_p and conductivity enter only the gas-temperature equation, which does not act back on the flow or the two-temperature model |
| D dynamic mesh | confirmed | `createMesh.H` builds a static `fvMesh`; the dynamic-mesh code was commented out | dead code and the `dynamicFvMesh`/`topoChangerFvMesh` libraries removed; `constant/dynamicMeshDict` is not read |
| D duplication | confirmed | | done: the 13 mass-weighted functions now use three helpers, with γ = C_p/C_v, C_p/C_pv and 1/W = ΣY_i/W_i of the mixture instead of mass-weighted γ and W; the two α sub-cycle branches merged; the repeated laser α test; the discarded C_e update; the no-op δ round trip; dead LTS code; the unused recoil clamp; the per-cell debug print; π and k_B to full precision; mojibake; tracker labels, thresholds and CSV header. Also done since: the laser pulse-centre and pulse-train duplicates (one helper each), the two recoil clamp blocks (one), the duplicated exchange cut-off checks (applied once), the vaporisation-temperature alias lookups (a loop), the `simulationType` parsing (only in the transport model now), `foundObject` lookups hoisted out of per-cell loops, the per-call Info message (once), and the tracker: plume detection from the gas temperature instead of T_l (held at ambient in gas), the "turbulent KE" now the kinetic energy of the deviation from the mean velocity vector (`fluctuation_KE_nJ`), `maxSource` labelled W/m³. **Not done:** replacing `postProcess.H` (openfoam.org 2016–2018 text) with the v2406 one. The local file defines `CREATE_FIELDS_DONE`, which `compInterFoam.C` relies on, so the swap needs a change to how the fields are created; it affects only the `-postProcess` option, not a run |

One finding was partly wrong: the Kapitza gas–metal exchange (A10) is not "in practice off"; it acts in partly filled interface cells below 3560 K. The other corrections: A3 is certain, not [UNCERTAIN], and did not act in the old run; A6 and A8 almost certainly did not act; A8's bounds did not move during the run; `phaseChangeMaxSourceDefault` is read (but overridden); and the review numbers two different items "B4".

### 10.2 Changes beyond the review

These also change results slightly and are in the corrected run:
- laser normalisation weighted by α (A10);
- recoil from T_melt (A10);
- `maxAlphaCo` limits Δt (does not bind here);
- one specific gas constant (174 J/(kg K)) for the evaporation flux, the recoil and the pocket; the recoil used k_B/m = 173.6.

### 10.3 Why the old solver still produced results

The defects fall into three groups, and none of them makes the equations unsolvable:
1. **Inactive in the reference run.** A3 (switched off by A4), A6 (limit never reached), A8 (bounds never reached), A5 (a no-op), A9 (a diagnostic) and the switched-off sub-models. They could not change the numbers. (The Kapitza exchange was active in partly filled interface cells below 3560 K; its α-bookkeeping error is still open.)
2. **Biases in a bounded, stable scheme.** A1 made the surface cells respond 2–3.5× more slowly to heating and cooling, and A2 gave the pocket 25–50 % of the evaporated mass. Both change *how much* drive the film gets, not whether the solver converges. The pocket pressure is also capped at p_sat, so where saturation binds, A2 has no effect. Implicit time stepping, bounded MULES transport and the pressure equation still conserve mass and momentum, which is why the conservation, Galilean and force-balance tests passed.
3. **Reporting and documentation.** B1 to B8 change what the docs say, not what the solver computes.

So the old results were plausible and internally consistent, but the size of the vapour drive, and therefore the launch speed, carried an unknown bias from A1 and A2. The rerun with the corrected solver measures it.

### 10.4 Tests of the corrected solver

| Test | Result |
|---|---|
| Compilation (OpenFOAM v2406, `wclean; wmake`) | no errors, no warnings |
| Galilean invariance (film at −170 m/s, 100 ps) | velocity exactly −170 m/s; pressure within 5×10⁻⁶ of 10⁵ Pa |
| 1-D laser and two-temperature column against the finite-difference code | peak electron temperature 6640 K against 6672 K (−0.5 %); electrons and lattice meet at 0.500 ps against 0.496 ps (+0.8 %): unchanged from before the fixes |
| Refactored code (all de-duplication and dead-code removal, two rounds) against the binary used for the rerun | Galilean test (0–100 ps) and the 1-D laser column (0–3 ps): every written field (U, p, p_rgh, α, T_e, T_l, T) bit-identical, so the refactoring changes no number |
| Reference case rerun | running |
