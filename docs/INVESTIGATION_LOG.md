# Mesh-study investigation log (November 2026)

This is a chronological record of the work to get mesh-converged film momentum and speed for `cases/fsLIFT_332mJ_Ti_meshStudy`. It includes the hypotheses that turned out wrong, so the reasoning can be checked. The defects themselves and their fixes are listed in [CODE_REVIEW.md](CODE_REVIEW.md) (items 8 and 9 of §1); this file explains how they were found.

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
