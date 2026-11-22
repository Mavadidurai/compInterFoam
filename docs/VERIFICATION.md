# Verification

1. [What is checked here, and what is not](#1-what-is-checked-here-and-what-is-not)
2. [The tests and their results](#2-the-tests-and-their-results)
3. [How to repeat the tests](#3-how-to-repeat-the-tests)
4. [Defects found and corrected](#4-defects-found-and-corrected)
5. [Open points](#5-open-points)

---

## 1. What is checked here, and what is not

| Question | Name | Answered here? |
|---|---|---|
| Does the code solve its equations correctly? | code verification | **largely**: the laser source, the two-temperature model (with k_e = C_e D_e), evaporation, recoil, surface tension, viscosity and Galilean invariance against exact solutions, an independent code or hand calculations (§2). Not covered by a test: the Drude k_e path, the gas–metal exchange, and the mass balance of the trapped-vapour pocket |
| How large is the numerical error of a production run? | solution verification | **partly**: a three-mesh study of the launch to 315 ps (§2.4), an energy and mass balance (§2.4), and drop benchmarks at several resolutions. The tear and the flight after about 4 ns have been computed on one mesh only (§5) |
| Are these the right equations for fs-LIFT? | validation | **no**: no comparison with measurements is included |

**Which code version the numbers belong to.** The rows marked *released code* were run with the code in this repository, including fixes 8 and 9 of [CODE_REVIEW.md §1](CODE_REVIEW.md#1-defects-that-change-results) unless stated. The other rows were run during development, before the corrections to the interface transport (the density update and the interface normal of the compression, §4) and before fixes 8 and 9. Those corrections do not touch the laser, evaporation or recoil code, and the 1-D two-temperature test repeated after fix 8 gave a peak electron temperature of 6649 K (6636 K before; finite-difference code 6672 K), so the earlier results for these parts still apply.

---

## 2. The tests and their results

### 2.1 Laser heating and two-temperature model

| Part | Test | Reference | Result |
|---|---|---|---|
| Laser energy | total energy deposited in a 1-D film | analytic value, 1473.7 J/m² | 1470.3 J/m² (−0.2 %) |
| Laser radial profile | deposited power against r at 0.5 ps (2-D) | exp(−2r²/w₀²), w₀ = 3 µm | largest deviation 0.0002 (relative) for r < 6 µm |
| Laser depth profile | decay of the deposited power with depth | exp(−a·d) | decay constant reproduced exactly inside the film; the first cell at the glass side is 1.2 % low |
| Electron and lattice temperatures | 1-D column against the independent finite-difference code `tools/ttm_fd_1d.py`, same inputs (*released code*) | finite-difference solution | peak electron temperature at the glass side 6636 K against 6672 K (−0.5 %); electrons and lattice meet (within 5 %) at 0.500 ps against 0.496 ps (+0.8 %) |
| Energy conservation of the two-temperature model | change of the electron and lattice energy against the absorbed energy, 1-D | exact balance | 99.5 % of the energy retained over 20 ps with conduction active |

### 2.2 Evaporation and recoil

| Part | Test | Reference | Result |
|---|---|---|---|
| Evaporation flux (Hertz–Knudsen) | every evaporating cell of a 2-D run (100 cells, 1986–5129 K), and a hand calculation | the formula | ratio 1.000000 (0.999999–1.000002); hand check 4249 against 4246 kg/(m² s) |
| Recoil pressure 0.54 · p_sat | every recoil cell of a 2-D run, and a hand calculation | the formula | ratio 0.9989; hand check 5.96 against 5.95 MPa. The 0.1 % difference comes from the gas constant: the recoil uses k_B/m = 173.7 J/(kg K), the evaporation the input 174 J/(kg K) |
| Latent heat | lattice energy lost against the latent plus sensible heat of the evaporated metal, 1-D | energy balance | 1.6 times the latent heat alone, against about 1.4 expected from the latent and sensible heat: 14 % apart, not explained |
| Force on the film | time integral of the logged force budget against the change of the film momentum, mesh-study case (*released code*) | Newton's second law | 0.969 (L2) and 0.977 (L3) over 2–315 ps; 0.959 on L1 over 0–305 ps. The residual shrinks with the cell size (§5). Before fix 9 the budget was 1.65 times the film's momentum change, because it also counted the push on the receiver; an earlier test on the superseded reference run (0.99, 0.96, 0.94 over 42–402 ps) used the recoil and vapour entries only |

### 2.3 Flow, interface and surface tension

Drop benchmarks: a liquid titanium drop (density 4515 kg/m³, viscosity 2.35e-3 Pa s, surface tension 1.64 N/m) of radius R = 0.2 µm in gas, with the laser, phase change, recoil, trapped vapour and gravity switched off. All other settings are those of the reference case. R/Δx is the number of cells per drop radius.

| Part | Test | Reference | Result (*released code*) |
|---|---|---|---|
| Surface tension, static | pressure jump of a static drop at 0.2 ns | Laplace law 2σ/R = 16.40 MPa | −2.75 % (R/Δx = 10), −3.80 % (20), −4.09 % (40). The error does not decrease under refinement (it grows slightly, to about −4 %): the curvature is computed from the gradient of α, a known limitation of this surface-tension method |
| Spurious velocity | largest velocity in the static drop | exactly 0 | 2.85 m/s (R/Δx = 10), 1.47 m/s (20), 1.40 m/s (40) |
| Metal volume | metal volume of the static and oscillating drops | exactly conserved | change below 1e-6 (static drop) and 3e-5 to 5e-5 (oscillating drop over 12 ns) |
| Inertia, surface tension and viscosity together | oscillation of a drop deformed into an ellipsoid (b/a = 1.2, mode l = 2) | Lamb: period 10.49 ns, decay time 15.4 ns | frequency 12.9 % below Lamb at R/Δx = 10 and 10.5 % below at R/Δx = 20 (approaching Lamb); damping rate 52 % and 70 % too low, not improving with refinement (§5) |
| Viscous term | single-phase flow along the axis in a pipe, u = U₀ J₀(kr), decaying in time | decay rate νk² = 1.0827e8 1/s | 1.0814e8 1/s over 5 ns (−0.12 %) |
| Continuity | divergence of the velocity in the film and in the drop | 0 | about 1e-3 1/s in the film; about 5 1/s in the drop, against a velocity gradient scale of 1e7 1/s |
| Galilean invariance | 1-D column (2.23 nm film cells) moving as a whole at −170 m/s, laser off | velocity −170 m/s and pressure 1e5 Pa everywhere, for all time | velocity −170.000 m/s; pressure within 0.1 % of 1e5 Pa at 10 ps and within 0.001 % at 30 and 100 ps; metal volume exactly conserved. Before the density correction (§4) the metal carried a false pressure of −1.2e8 to −1.6e8 Pa, close to −ρU² = −1.3e8 Pa |
| Interface compression on a moving interface | 1-D film with a smeared 40 nm interface layer, moving at −170 m/s; one setting changed at a time | the compression must re-sharpen the layer | after the correction of the interface normal (§4): 19 partly filled cells at 0.5 ns with `cAlpha 2` against 26 with `cAlpha 0`. Before the correction, every setting gave identical results |
| Parallel and restart | serial against 3 and 4 processes with different decompositions; restart from a written time | identical results | identical to 8 digits; a restart reproduces the uninterrupted run |

### 2.4 Solution verification of the mesh-study case

Case `cases/fsLIFT_332mJ_Ti_meshStudy`, meshes L1 / L2 / L3 with 62,010 / 124,200 / 248,252 cells (refinement ratio √2 in both directions), all fixes (*released code*). Full record: [INVESTIGATION_LOG.md](INVESTIGATION_LOG.md) §9 and §11; data: [`figures/fsLIFT_332mJ_Ti_meshStudy/data/fixes8and9/`](../figures/fsLIFT_332mJ_Ti_meshStudy/data/fixes8and9/).

| Quantity at 315 ps | L1 / L2 / L3 | Result |
|---|---|---|
| Film momentum | 9.116 / 9.104 / 9.061 × 10⁻¹³ kg m/s | spread 0.6 %; the change grows on the finer pair, so no order or GCI can be given |
| Speed within 1 µm of the axis | 396.3 / 389.6 / 389.5 m/s | converged, L2–L3 0.03 % |
| Film thickness on the axis | 68.84 / 68.88 / 68.90 nm | converged, GCI 0.02 % |
| Peak lattice temperature | 19,183 / 19,778 / 20,044 K | monotonic, GCI 1.3 % |
| Kinetic energy | 1.300 / 1.285 / 1.267 × 10⁻¹⁰ J | spread 2.6 % |

| Balance (L1) | Result |
|---|---|
| Laser energy deposited | 326.04 pJ against 326.24 pJ requested (wedge), −0.06 % |
| Evaporation heat against the latent heat of the mass lost | 92 % at 30 ps |
| Energy not accounted for | −0.8 % (0.51 ps), −2.3 % (30 ps), −4.0 % (315 ps), −7.1 % (2.16 ns) of the laser energy (§5) |
| Metal volume | −0.098 % over 0–315 ps, all by evaporation |
| Time step | L2 at three largest time steps (0.595, 1.0, 1.68 ps): momentum within 0.4 %, GCI 0.2 % (before fix 8) |

Figures: [`figures/verification/static_drop.png`](../figures/verification/static_drop.png), [`oscillating_drop.png`](../figures/verification/oscillating_drop.png), [`ttm_column_vs_fd.png`](../figures/verification/ttm_column_vs_fd.png).

---

## 3. How to repeat the tests

Load the OpenFOAM environment and build the solver first ([USER_GUIDE.md §1](USER_GUIDE.md#1-building-the-solver)). Run the commands from the repository root; `RUNS` is any empty directory.

**Drop benchmarks** (the case builder copies the reference case and switches the laser and phase change off):

```bash
RUNS=$HOME/benchmarks; mkdir -p $RUNS && cd $RUNS
B=/path/to/compInterFoam-fsLIFT/tools/benchmarks
# static drop at R/dx = 10, 20, 40: <dir> <cell size nm> <aspect b/a> <end time s> <largest time step s>
bash $B/make_drop.sh sd_dx20 20 1 4e-10 2e-12 && (cd sd_dx20 && compInterFoam > log.run)
bash $B/make_drop.sh sd_dx10 10 1 4e-10 2e-12 && (cd sd_dx10 && compInterFoam > log.run)
bash $B/make_drop.sh sd_dx5   5 1 4e-10 1e-12 && (cd sd_dx5  && compInterFoam > log.run)
# oscillating drop at R/dx = 10 and 20
bash $B/make_drop.sh os_dx20 20 1.2 1.2e-8 1e-11 && (cd os_dx20 && compInterFoam > log.run)
bash $B/make_drop.sh os_dx10 10 1.2 1.2e-8 5e-12 && (cd os_dx10 && compInterFoam > log.run)
for c in sd_dx20 sd_dx10 sd_dx5 os_dx20 os_dx10; do xvfb-run -a pvpython $B/drop_analyze.py $c; done
python3 $B/fit_osc.py os_dx20 os_dx10        # frequency and damping against Lamb
python3 $B/plot_benchmarks.py $RUNS          # figures/verification/static_drop.png, oscillating_drop.png
```

Run times on one core: 0.5 min (sd_dx20), 2 min (sd_dx10), about 10 min (sd_dx5), 3 min (os_dx20), about 25 min (os_dx10).

**Viscous term:**

```bash
bash $B/make_drop.sh vt_dx20 20 1 5e-9 1e-11
cd vt_dx20 && sed -i 's/^writeFormat .*/writeFormat     ascii;/' system/controlDict
rm -rf 0 && cp -r 0.orig 0 && postProcess -func writeCellCentres -time 0 && cd ..
python3 $B/make_vt.py vt_dx20 && (cd vt_dx20 && compInterFoam > log.run) && python3 $B/vt_analyze.py vt_dx20
```

**1-D columns:** see [USER_GUIDE.md §3](USER_GUIDE.md#3-the-one-dimensional-test-cases).

**Force on the film** in a run with `momentumBudget true`:

```bash
xvfb-run -a pvpython tools/metal_momentum.py cases/fsLIFT_332mJ_axisymmetric 4.2e-11,1.02e-10
```

The change of the printed momentum between two times is compared with the time integral of the `recoil` entry of the `momentumBudget` lines in the solver log over the same interval.

---

## 4. Defects found and corrected

These defects were found by the tests above, or by inspecting development runs. The released code contains every correction. The results of the reference case in `figures/fsLIFT_332mJ_axisymmetric/` were produced before the corrections of [CODE_REVIEW.md §1](CODE_REVIEW.md#1-defects-that-change-results) and are superseded ([RESULTS.md §3](RESULTS.md#3-the-earlier-reference-case-superseded)).

| Defect | How it was found | Correction |
|---|---|---|
| The recoil force was applied at half strength: the recoil pressure sits in the surface cell only, and the linear face average gave each of the two interface faces half of it | the force balance on the film | face value from the larger neighbour (`recoilFaceInterpolation localMax`, default) |
| A false film edge at the outer boundary: an inflow condition with α = 0 cut the film open there, and the surface tension of that edge pushed on the film | the force balance; a per-step force budget located the force at r = 7–8 µm | `zeroGradient` for α on the outer boundary |
| Surface tension applied in a band only one cell wide around the interface missed about 17 % of the force, at every resolution | static and oscillating drop | band three cells wide (`surfaceTensionBandLayers 3`, default) |
| The interface compression flux was transported like a velocity: α became unbounded (up to 2.7), and the clipping that followed deleted 9 % of an oscillating drop's metal in 12 ns | metal volume of the oscillating drop | standard OpenFOAM form of the flux (`stockCompressionFlux true`, default); the volume is now exactly conserved |
| The density was updated inside the sub-cycles of the α equation. A moving interface then carried a false pressure of −ρU² | Galilean test; a run-away of the leading film surface | density updated once, after the sub-cycles |
| The interface compression used an interface normal built once at the start of the run, so it had no effect on a moving interface. The tip of the moving dome smeared out to α ≈ 0.5 | compression test on a moving interface; a smeared dome tip | the current interface normal is used |
| The optional pressure clamp overwrote the pressure after the pressure solve, so the velocity no longer satisfied continuity, and a development run crashed | a crash at 647 ps | `pressureClamp false` |
| The electron conductivity k_e = C_e · D_e grows with T_e even when electrons and atoms are at the same temperature: 10–17 times too high in a hot film, so the film first moved the wrong way | comparison with the finite-difference code | Drude form k_e = k_e0 · T_e/T_l (option `ke0`) |
| Temperature limits of 10,000 K were reached with the 10 nm absorption depth and removed energy | energy balance | limits raised to 50,000 K |
| The first optical inputs were wrong: reflectivity 0.35 instead of 0.428, and an absorption coefficient of 1.07e7 1/m (93 nm depth) of unknown origin (it is not 2k/λ = 1.65e7, k/λ or 2πk/λ) instead of 4πk/λ (1.035e8 1/m, 9.7 nm) | a check of the inputs against the optical constants of titanium | corrected inputs in the reference case |
| The temperature equations were skipped after the pulse: a solver tolerance of 1e-6 to 1e-8 was above the residual, so 0 iterations were made and the temperature froze | 1-D energy budget; solver logs | tolerance 1e-14 |
| With the physical laser inputs the run failed within 2–6 ps | development runs | cause: the mass source of the unverified phase-explosion sub-model; switched off |
| One trapped-vapour pocket for the whole spot pushed a 7.4 µm-wide region as a rigid plug, which tore at the cold edge | development runs | one vapour column per radial ring (`mode local`) |
| The wedge received the energy of a bounding box rather than of its 5° sector | laser energy test | exact sector share (`beamCoverage geometric`, default) |
| Forced iterations of the α solver on an already converged system produced NaN | development runs | `minIter 0` |
| T_e and T_l were carried by the mixture flux. At a receding metal surface it carried gas cells (reset to ambient temperature) into the hottest metal cells, at a rate set by the cell size: the film momentum changed by about 6 % per mesh level | mesh study; a rerun without temperature transport removed the mesh dependence | carried by the metal flux of the α equation in advective form (`temperatureAdvectionFlux metal`, default); the Galilean heat loss of a moving film fell from 3.0 % to 0.06 % in 90 ps (fix 8) |
| The trapped-vapour pressure was also applied to the receiver's top surface, and the receiver area was added to the pocket rings | the logged force budget was 1.65 × the film's momentum change; the receiver gained momentum before any signal could cross the gap | only surfaces above `filmMinHeight` belong to the pocket (fix 9); the budget then closes within 2–4 % |

---

## 5. Open points

1. **Damping of the drop oscillation.** The oscillation frequency approaches Lamb's value as the mesh is refined, but the damping is several times too weak. The viscous term alone is exact (§2.3), so the cause is in the interface numerics: the spurious velocities of the surface-tension discretisation are as large as the oscillation velocity at R/Δx = 10. A more accurate curvature (height functions, or a geometric interface reconstruction) would be needed. This matters for any study of how a liquid neck breaks.
2. **Tear and flight on one mesh only.** The mesh study covers 0–315 ps. The tear of the dome wall (about 4 ns) and the flight were computed on the coarsest mesh only, where the wall is about one cell thick at the tear and the cap is thinner than one cell from about 7 ns. A run refined along the flight path is under way; the radial resolution of the tear is not yet addressed ([RESULTS.md §2](RESULTS.md#2-runs-in-progress)).
3. **Interface compression settings.** `cAlpha 2`, `compressionLimiter 0.5` and the compression band were chosen during development, before the compression worked on moving interfaces, and have not been re-tuned since. With these values the effective compression coefficient is min(2, 0.5) = 0.5 ([MODEL.md §7.1](MODEL.md#71-volume-fraction-transport)).
4. **Residual of the force budget.** The extra 60–75 % of the logged budget seen in the superseded run was the pocket pressure on the receiver (fix 9). With the fix the logged total closes to 0.96–0.98 of the film's momentum change; the residual shrinks with the cell size, which fits the α weighting of the interface cells, but this is not proven.
5. **Energy residual.** About 4 % of the laser energy is not accounted for at 315 ps and 7 % at 2.16 ns (§2.4). It is taken to be the mismatch between the temperature transport with the metal flux and the transport of α between steps, two terms of about 40 % of the laser energy each that cancel in exact arithmetic. This is not proven. Its effect on the momentum has not been measured, and the balance has not been repeated on L2 and L3.
