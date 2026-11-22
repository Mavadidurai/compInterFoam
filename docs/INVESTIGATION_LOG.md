# Investigation log of the mesh-study case (November 2026)

This is a chronological record of the work on `cases/fsLIFT_332mJ_Ti_meshStudy`: first mesh-converged film momentum and speed, then the energy balance, the flight on the coarse mesh and a run refined along the flight path. It includes the hypotheses that turned out wrong, so the reasoning can be checked. The defects themselves and their fixes are listed in [CODE_REVIEW.md](CODE_REVIEW.md) (items 8 and 9 of §1); this file explains how they were found.

Units: momentum in 10⁻¹³ kg m/s (full circle), speeds in m/s, times in ps unless stated. Meshes L1 / L2 / L3 have 62,010 / 124,200 / 248,252 cells (refinement ratio √2). All runs on 4 processes, OpenFOAM v2406.

---

## 1. First mesh study to 1 ns (PR #3, merged 17 Nov)

- `tools/mesh_study_gci.py` labelled convergence the wrong way round: the ratio e32/e21 is above 1 for a converging sequence, and the tool called that case "diverging". Fixed in PR #3.
- Result at 1 ns: the peak temperatures converged, but momentum (5.607 / 5.322 / 5.013 on L3 / L2 / L1, order p = 0.23), kinetic energy (p = 0.16) and the speed within 1 µm of the axis (234.9 / 221.9 / 209.5, diverging) did not.

## 2. Phase A: time step against mesh (to 315 ps)

**Hypothesis (wrong):** the time-step error dominates. The coarse levels used larger time steps.

- L2 with maxDeltaT 0.595 / 1.0 / 1.68 ps: momentum 5.317 / 5.324 / 5.337, p = 1.15, GCI 0.2 %. The time-step error is small.
- All three meshes at maxDeltaT 0.595 ps: momentum 5.607 / 5.317 / 5.010 (L3 / L2 / L1), p = 0.16. Axis speed 236.1 / 222.6 / 210.1. The mesh error stays at about 6 % per level.
- **Conclusion:** the time step was not the cause. The non-convergence is a mesh effect.

## 3. Where the meshes part: force budget and temperatures

- The `momentumBudget` log lines, integrated over time, gave the same impulse on all meshes up to 30 ps. The meshes parted after that, in the recoil and pressure terms.
- The integrated budget was about 1.65 × the film's actual momentum change on every mesh. It was then used only to compare meshes; the cause was found later (§8).
- The hottest metal cooled faster on coarse meshes. Peak lattice temperature at 60 ps: 9503 / 9924 / 10284 K (L1 / L2 / L3); recoil maximum 547 / 610 / 703 MPa. A hotter surface gives more recoil for longer, hence more momentum on fine meshes.

**Hypotheses examined for the cooling:**
- *Heat leaking into the gas* (inferred from an "E_gas" measure): **retracted.** That measure mostly tracked the volume of the interface cells, not heat.
- *Gas–metal Kapitza exchange*: **ruled out.** It switches off above `gasMetalCutoffTemperature` (3560 K), and the cooling is in cells far hotter than that.

## 4. Phase B: temperature transport switched off (to 110 ps)

- L1 and L2 rerun with `advectTemperatures false`. Gap between L1 and L2 at 110 ps, before → after:
  - momentum +4.9 % → +0.7 %
  - kinetic energy +10.9 % → +0.8 %
  - speed within 1 µm +4.2 % → +0.07 %
- **Conclusion:** the transport of T_e and T_l causes the mesh dependence. It used the mixture flux φ, which at a receding surface carries gas cells (reset to ambient temperature) into the hottest metal cells, at a rate set by the cell size.

## 5. Fix 8: transport with the metal flux (commit `3b42c29`)

- The energy equations now use the metal flux of the α equation (`alphaPhi10`, summed over the sub-cycles) in advective form, `C[∇·(φ_α T) − T ∇·φ_α]`. Switch `temperatureAdvectionFlux metal|mixture`.
- Galilean tests (film at −170 m/s, 90 ps), heat change:

| Case | Heat change |
|---|---|
| Uniform 1500 K, mixture flux | −3.0 % |
| Uniform 1500 K, metal flux | −0.064 % |
| Step 1000 / 1800 K, at rest (reference) | +1.38 % |
| Step, metal flux | +1.17 % |
| Step, mixture flux | −2.47 % |

- The +1.38 % at rest is physical: electrons pass heat to the lattice. The exact mechanism is not confirmed. The metal flux tracks the reference within 0.2 % and creates no heat.
- 1-D laser column against the finite-difference code: peak T_e 6649 K against 6672 K, coupling time 0.500 ps against 0.496 ps. This is unchanged, but the film barely moves there, so it does not test the transport.

## 6. Phase C: mesh study with fix 8 (to 315 ps, 17–18 Nov)

- The container restarted during L2. L2 resumed from 210.6 ps, its last complete write, which was checked in every processor directory.
- At 315 ps (L1 / L2 / L3):

| Quantity | L1 | L2 | L3 | Spread |
|---|---|---|---|---|
| Momentum | 7.388 | 7.376 | 7.348 | 0.55 % (11.9 % before) |
| Kinetic energy (10⁻¹¹ J) | 8.331 | 8.223 | 8.121 | 2.6 %, not converging |
| Speed within 1 µm | 316.7 | 310.9 | 310.9 | converged |
| Speed within 25 µm | 118.6 | 121.2 | 118.5 | oscillatory |

- Momentum was 47 % higher than with the old transport on L1, and still rising at 315 ps.
- Data committed as `figures/fsLIFT_332mJ_Ti_meshStudy/data/metalFlux/`; PR #4 opened. These data were later superseded by §8.

## 7. Cross-checks asked for

- The 1-D laser column was rerun (§5).
- The budget factor of 1.65 was the one unexplained inconsistency, so it was investigated next.

## 8. Fix 9: the vapour pocket also pushed the receiver (commit `ec3a64d`, 19 Nov)

**Steps:**
1. Window by window, the budget was 1.6–1.8 × the film's momentum change from the first window on, on all three meshes. That pointed to a fixed factor, not a slow leak.
2. The momentum extraction counts only y > 8.05 µm. The budget counts every metal cell, including the receiver (α = 1, 1606 µm³). Film plus receiver change, −12.285, matched the budget, −11.985, within 2.4 %.
3. The receiver held −5.0 momentum at 315 ps (0.68 × the film). It already had −1.8 at 28 ps, before any signal could cross the 12 µm gas gap (about 35 ns at the speed of sound in the gas).
4. The receiver's top cells, at 300 K and not evaporating, carried up to 1.13 GPa of "recoil" pressure at 28 ps. That equalled the film maximum.
5. In the code, `glassSide()` selects every surface with metal below and gas above. `filmMinHeight` was used only for the film velocity average. So the receiver top got the pocket pressure, and its area was added to each pocket ring, roughly doubling the ring volume.

**Alternative explanation considered:** the incompressible gas in the gap carries the film's push to the receiver. Rejected: the gas next to the receiver was below ambient pressure (0.6 bar), and the 840 MPa jump at the receiver surface needs an interface force.

**Test (L1, 0–5 ps, same inputs with and without the fix):**
- Pocket pressure on the receiver: 485 MPa → 0.
- Receiver momentum: −0.0814 → −0.00037.
- Budget per 1 ps window: matched film + receiver within 0.2 % before the fix, and the film alone within 0.4 % after.
- The film was identical up to 3 ps.

**L1 at 305 ps, before → after fix 9:**
- momentum 7.372 → 9.101 (+23.5 %)
- kinetic energy +56 %
- speed within 1 µm 315.4 → 395.0 (+25 %)
- peak temperature unchanged
- The budget is now within 4.1 % of the film's momentum change (0–305 ps). The residual is attributed to the α weighting of interface cells; this is not confirmed.

**Consequences:**
- The Phase C data and the reference-case results in `RESULTS.md` are too low.
- Fixes 8 and 9 together raised the L1 momentum at about 0.3 ns by roughly 80 % and the axis speed by roughly 90 %.

## 9. Phase D: mesh study with fixes 8 and 9 (20 Nov)

- L1 is the fixed on/off-test run continued to 315 ps. L2 and L3 were rerun from 0, because the bug acts from about 3 ps. Driver: `phaseD/run.sh`, restartable (L3 needed three 2-hour runs).
- At 315 ps (L1 / L2 / L3):

| Quantity | L1 | L2 | L3 | Result |
|---|---|---|---|---|
| Momentum | 9.116 | 9.104 | 9.061 | spread 0.6 %; the change grows on the finer pair, so no order or GCI |
| Kinetic energy (10⁻¹⁰ J) | 1.300 | 1.285 | 1.267 | spread 2.6 %; same pattern, apparent p = 0.33 |
| Speed within 1 µm | 396.3 | 389.6 | 389.5 | converged (L2–L3 0.03 %) |
| Speed within 25 µm | 146.4 | 149.6 | 146.2 | oscillatory, ±1.2 % |
| Thickness on the axis (nm) | 68.84 | 68.88 | 68.90 | converged, GCI 0.02 % |

- Spread = (max − min)/min: momentum (9.116 − 9.061)/9.061 = 0.61 %.
- Momentum budget against the film's momentum change over 2–315 ps: 0.969 (L2), 0.977 (L3); 0.959 on L1 over 0–305 ps. The residual shrinks with the cell size, which fits the α weighting of interface cells, but this is not proven.
- Same pattern as Phase C (fix 8 only: 7.388 / 7.376 / 7.348, spread 0.55 %). Fix 9 raised momentum by 23–24 % on every level without making it mesh-dependent.
- Not settled: the film is still accelerating at 315 ps (L1 9.056 at 280 ps → 9.116 at 315 ps), so this is not the final launch momentum.
- Data: `figures/fsLIFT_332mJ_Ti_meshStudy/data/fixes8and9/`.

## 10. Known limits not addressed by these fixes

- **The gas is incompressible** (`rhoConst`, a documented verification choice). The gap reacts instantly: at 28 ps the gap shows +58 bar under the film and −58 bar above it. A real gas cannot have negative pressure, and a real gap would not feel the film this early.
- **Peak temperatures exceed the closure range.** The lattice reaches about 20,000 K at 3.3 ps, far above the critical point of titanium. The evaporation and recoil closures use at most `maxPhysicalTemperature` (10,000 K), so they are outside their range there.
- **No comparison with experiment.** Mesh convergence shows that the equations are solved consistently, not that the model is right.
- **A small notch in the film surface on the axis** shows in the animation frames and in the unprocessed cell data, so it is in the solution, not a rendering artefact. Cause not investigated.

## 11. Energy and mass balance (L1 with fixes 8 and 9, 22 Nov)

- Run: L1 from 0 with `energyDiagnostics true` (`meshstudy/energyTest`), 0–315 ps, then continued towards 20 ns. It matches the Phase D L1 run within 0.04 % at all common times (momentum 9.112 against 9.116 × 10⁻¹³ kg m/s at 315 ps).
- Balance per step from the solver's energy diagnostics: first pass = energy before, last pass = energy after, laser and phase-change terms once per step. "Solve" = change inside the temperature solve minus the sources; "between steps" = change from one step's end to the next step's start (transport of the metal fraction).

| Up to | E − E₀ | Laser | Evaporation | Solve | Between steps | Net unaccounted |
|---|---|---|---|---|---|---|
| 0.51 ps | 321.90 | 324.51 | −0.00 | −2.61 | −0.00 | −2.61 (−0.8 %) |
| 30 ps | 315.41 | 326.04 | −3.06 | +8.92 | −16.50 | −7.58 (−2.3 %) |
| 315 ps | 306.33 | 326.04 | −6.64 | +132.92 | −146.00 | −13.08 (−4.0 %) |
| 2.16 ns | 294.23 | 326.04 | −8.77 | +132.48 | −155.52 | −23.04 (−7.1 %) |

All values in pJ for the 5° wedge; percentages relative to the laser energy.

- **Laser deposition is correct:** 326.04 pJ against 4.1089e-8 J × 0.57168 / 72 = 326.24 pJ (−0.06 %).
- **Evaporation heat is consistent with the mass lost:** −3.06 pJ at 30 ps is 92 % of the latent heat of the metal the film lost (3.32 pJ).
- **Mass:** the film lost 0.0140 µm³ (0.098 %) over 0–315 ps, so evaporation cannot change the momentum by more than about 0.1 %.
- **The two large terms nearly cancel.** Each is about 40 % of the laser energy by 315 ps. The solve carries the temperature with the metal flux (fix 8) and the metal fraction moves between steps; in exact arithmetic these cancel. The net −4.0 % at 315 ps is taken to be the discretisation mismatch between them, not proven: part could be heat leaving the metal that the diagnostics do not count.
- The net loss keeps growing after the launch (−7.1 % at 2.16 ns), while most of the momentum is gained in the first tens of picoseconds (−2.3 % by 30 ps). Its effect on the momentum has **not** been measured; a loss that is similar on every mesh would not show in the mesh study.
- Not done: the same balance on L2 and L3, and an energy-conserving form of fix 8.

## 12. L1 continued after 315 ps: dome, tear and flight (21–22 Nov)

- The energy-diagnostics L1 run (§11) was continued: to 1 ns at maxDeltaT 1 ps (written every 20 ps), then at 2.5 ps (every 0.1 ns). Background jobs were limited to 2 hours, so the driver `energyTest/flight.sh` writes and stops at its time budget and is restarted from the latest write. Data: `figures/fsLIFT_332mJ_Ti_meshStudy/data/flightL1/`.
- **The film keeps accelerating after 315 ps:** speed within 1 µm of the axis 396.1 m/s at 315 ps, 425.8 at 1 ns, highest 432.9 at 1.9 ns, then 413.8 at 6.5 ns. Momentum 9.11 → 9.44 (1 ns) → 9.66 (4 ns). So the 315 ps mesh study converges the early launch, not the final launch momentum.
- **The dome wall breaks between the writes at 3.96 and 4.26 ns.** Its thickness there is about 43 nm by a volume estimate, against 49 nm radial cells: about one cell. The fragments along the wall in the pictures are thinner than a cell; they are taken to be mesh artefacts.
- **The central cap flies free** at about 410–420 m/s, 62–68 nm thick on the axis.
- **From about 7 ns the cap is not resolved.** It enters the graded cells (y 16.5–17 µm) and then the 71 nm cells below. `y_bot_axis` is NaN from 7.195 ns, and from about 9 ns no cell for r < 1 µm below y = 18.5 µm has α above 0.78. The mean speed of that metal falls from 355 m/s (7 ns) to 266 m/s (18.7 ns); real and numerical slowing cannot be separated on one mesh.
- **Momentum rises after the tear** (9.66 at 4 ns → 10.39 at 7.6 ns) while the kinetic energy stays at 1.41–1.44 × 10⁻¹⁰ J. Not analysed.
- At 18.7 ns the lowest metal (α > 0.5) is at 10.52 µm, 2.52 µm above the receiver. Contact is expected at about 24–28 ns; the run end was moved from 20 to 30 ns (run A).
- The animation (`LIFT_launch_and_tear_0-6.5ns`) stops at 6.5 ns, before the cap reaches coarse cells.

## 13. Run C: refined flight path (22 Nov)

**Aim:** resolve the cap through the flight and the impact without rerunning the verified launch.

**Mesh:** L1 with the gap block changed from `((8.5 120 1) (0.5 17 0.12) (3.0 354 1))` to `((8.5 567 1) (0.5 43 0.567) (3.0 354 1))`, i.e. 15 nm cells for y = 8–16.5 µm (71 nm before), graded to the unchanged 8.5 nm cells from y = 17 µm. 106 × 1058 = 112,148 cells; `checkMesh` OK. The receiver, film, glass-side and radial cells are unchanged.

**Start state:** L1 at 2.9625 ns, the last write before the tear. At that time all of the donor metal is above y = 17 µm, where the two meshes have identical cells.

**Mapping:**
1. `mapFields -consistent` of v2406 offers only `mapNearest`, `interpolate` and `cellPointInterpolate`, and maps a field only if the target already has it. With an empty target time directory and `-parallelSource` it mapped nothing.
2. Therefore: the source time was reconstructed in a scratch copy (`reconstructPar`, then `foamFormatConvert` to ASCII). For every volume field a target file was written with the source boundary conditions and a uniform internal field. `mapFields -consistent -mapMethod mapNearest` then mapped all 31 volume fields. The surface fields `phi` and `alphaPhi10` are not mapped; the solver rebuilds them (`READ_IF_PRESENT`).
3. `uniform/vapourPocket` (the pocket state per 0.1 µm ring; mesh-independent) was copied from the source, then the case was decomposed onto 4 processes.

**Checks:**
- Source against target, full circle: metal volume 1620.771214 µm³, momentum −9.623834e-13 kg m/s, Σ α V T_l 4.930897e-13, gas mass 7.522672e-15 kg, total mass 7.308792e-12 kg, p range and peak T_e: identical in every printed digit.
- Boundary conditions identical, except the stored `value` of the mixed conditions and of `rho` on the receiver bottom (y = 0), which the solver recomputes on the first step.
- 5 steps on 4 processes without error; the pocket state was read (m_v = 3.81e-19 kg, depth on the axis 1.47 µm, T_v = 4562 K).

**Run:** driver `flightC/flightC.sh`, to 35 ns at maxDeltaT 2.5 ps, written every 0.2 ns (about 5 GB). It starts when run A has finished. The speed was not measured: the smoke test shared the 4 cores with run A. The estimate is about 30–35 h from 1.8 × the cells at the same time step [UNCERTAIN].

**What C does not fix:** the radial cells are unchanged (49 nm), so the tear of the dome wall after 2.96 ns is as mesh-limited as on L1.

## 14. Why the film momentum keeps rising after the tear (22 Nov)

**Observation (§12):** on L1 the film momentum rises from 9.66 (3.96 ns) to 10.39 (7.6 ns) × 10⁻¹³ kg m/s, while the kinetic energy stays at 1.41–1.44 × 10⁻¹⁰ J and the cap slows down.

**Step 1: the logged budget does not close after about 1 ns.** Integrated over windows, the α-weighted `momentumBudget` total gives −0.30 against a measured −0.33 over 0.32–1 ns, but the opposite sign from 2 to 6 ns (+0.06 to +0.10 per ns against −0.05 to −0.18). It weights each cell's force by α. In a partly filled cell, though, almost all the mass is metal (ρ_m/ρ_g ≈ 3900), so the metal feels almost all of the cell force. After the tear much of the metal is in such cells, so the logged budget cannot be used.

**Step 2: where the momentum is gained** (`data/flightL1/momentum_by_alpha_and_radius.csv`, units 10⁻¹³ kg m/s):

| | 2.96 ns | 7.90 ns |
|---|---|---|
| Film | −9.61 | −10.47 |
| in cells with α > 0.9 | −3.78 | −1.05 |
| in cells with 0.1 < α < 0.5 | −0.69 | −5.70 |
| metal volume in cells with 0.01 < α < 0.5 | 0.24 µm³ | 1.04 µm³ |
| mass-weighted speed of that metal | 71 m/s | 143 m/s |

The metal is spread more and more over partly filled cells, and that metal speeds up.

**Step 3: which force** (`tools/force_decomposition.py`; the force fields written at every 0.1 ns, weighted by the metal mass fraction α ρ_m/ρ and integrated over time; units 10⁻¹³ kg m/s):

| Window | Measured change | Surface tension | Pocket pressure | Gas pressure | Sum | Surface tension summed over all cells |
|---|---|---|---|---|---|---|
| 0.42–3.86 ns (before the tear) | −0.435 | +0.123 | −0.588 | −0.073 | −0.537 | **+0.006** |
| 3.96–7.59 ns (after the tear) | −0.726 | −0.324 | −0.128 | −0.096 | −0.549 | **−0.570** |

Split of the surface tension after the tear (mass-weighted): +0.319 on the cap (r < 1.5 µm, y < 19.5 µm), −0.478 on partly filled metal (α < 0.5) elsewhere, −0.165 on full cells elsewhere.

**Conclusion:**
1. **About 45 % of the rise (−0.32 of −0.73) is a numerical net force of the surface tension.** The surface-tension force (continuum surface force, σκ∇α) summed over all cells is zero before the tear (+0.006, as it must be for closed interfaces) and −0.57 after it. The force sits mainly on smeared metal thinner than a cell: the broken wall and its fragments. The continuum-surface-force form does not conserve momentum when the curvature is computed from the gradient of a smeared α. Why the net force points towards the receiver has not been established [UNCERTAIN].
2. **About 18 % (−0.13) is the trapped-vapour pocket**, which keeps pushing the descending metal after the dome has opened. The model has no venting: once the wall breaks, the vapour should escape through the opening, but the pocket rings keep their mass and push every glass-side surface in their ring.
3. **About 13 % (−0.10) is gas pressure.**
4. **About 24 % (−0.18) is not accounted for.** The forces flicker from step to step (standard deviation 0.155 per sample against a mean of −0.09), so integrating 39 snapshots is uncertain by about ±0.1. The viscous force is not written and is not included.
5. **Kinetic energy stays flat** because the cap loses energy to the pull of the wall (surface tension +0.32 on the cap), while the low-mass smeared metal gains momentum.

**Consequences:**
- The momentum after the tear (10.39 at 7.6 ns) is not a physical result. The physical momentum gain after about 4 ns is at most the pocket and gas-pressure part, and the pocket part assumes vapour that should have vented.
- The logged `momentumBudget` (α-weighted) must not be used after the tear.
- Run C refines only the vertical cells along the flight path. The fragments of the wall are thinner than the radial cells, so C will not remove this force [UNCERTAIN: C may reduce it on the cap].

**Possible remedies (not done):** vent the pocket rings once the film above them has opened; resolve the wall radially (10–15 nm cells for r < 3.5 µm from about 3 ns); a momentum-conserving surface-tension form, or a curvature from height functions.
