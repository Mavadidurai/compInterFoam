# User guide

1. [Building the solver](#1-building-the-solver)
2. [Running the reference case](#2-running-the-reference-case)
3. [The one-dimensional test cases](#3-the-one-dimensional-test-cases)
4. [The case settings explained](#4-the-case-settings-explained)
5. [Analysis and plotting tools](#5-analysis-and-plotting-tools)
6. [Setting up a new case](#6-setting-up-a-new-case)

The physics behind the settings is in [MODEL.md](MODEL.md).

---

## 1. Building the solver

The solver needs **OpenFOAM v2406** from ESI/OpenCFD (openfoam.com). The Foundation releases (openfoam.org) have different interfaces and are not supported.

```bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc   # path of your installation
cd src
wmake                                             # about 1 minute
```

The executable `compInterFoam` is installed in `$FOAM_USER_APPBIN`, which is on the `PATH` once the OpenFOAM environment is loaded. The compiler prints a few warnings about unused variables; they are harmless.

---

## 2. Running the reference case

### 2.1 With the run script

```bash
tools/run_case.sh cases/fsLIFT_332mJ_axisymmetric [number of processes, default 4] [end time in s, default 2e-8]
```

The script:
1. creates the mesh (`blockMesh`), copies the initial fields from `0.orig` to `0`, and creates the receiver zone (`topoSet`) and the initial phases (`setFields`);
2. decomposes the case for the requested number of processes (`decomposePar`);
3. runs the solver in three stages, setting `endTime`, `maxDeltaT` and `writeInterval` in `system/controlDict` before each:

| Stage | Simulated time | Largest time step | Fields written every | Why |
|---|---|---|---|---|
| 1 | 0–2 ps | 10 fs | 0.5 ps | laser pulse and electron–lattice coupling |
| 2 | 2 ps – 1 ns | 1 ps | 20 ps | launch of the film |
| 3 | 1 ns – end | 2.5 ps | 0.2 ns | flight of the film; the capillary limit is about 4.4 ps |

Within each stage the time step is usually smaller than the largest value, because the Courant limits in `system/controlDict` apply.

**Progress** is written to `run_status.txt` in the case directory, with a line every minute. Each stage has its own log (`log.stage1`, `log.stage2`, `log.stage3`).

**Droplet detection.** Every minute the script runs `tools/detect_droplet.py`. It reports a metal piece on the beam axis that has a gas gap of at least 20 nm to all other metal, and writes detections to `droplet_detections.txt`. It only reports, and never stops the run. It looks along the beam axis only: a piece that separates by tearing away from the axis is not seen. A smeared interface can also give a false detection. Confirm every detection in ParaView.

**Run time and disk space** for the reference case on 4 processes: about 20 minutes for stage 1, and about 20 hours in all to 20 ns. The written fields take about 5 GB.

### 2.2 Continuing an interrupted run

```bash
RESUME=1 tools/run_case.sh cases/fsLIFT_332mJ_axisymmetric 4 2e-8
```

The script finds the latest written time in `processor0` (or the case directory for a serial run) and continues with the right stage. A run stopped during stage 1 (before 2 ps) should be started again from the beginning instead. The log of the interrupted stage is kept as `log.<stage>.part1`, `part2`, … and the analysis tools read all parts in order. Before resuming, check that the last written time directory is complete in every `processor*` directory; delete it if a write was interrupted. The trapped-vapour state is restored from `<time>/uniform/vapourPocket`.

### 2.3 By hand

```bash
cd cases/fsLIFT_332mJ_axisymmetric
blockMesh && rm -rf 0 && cp -r 0.orig 0 && topoSet && setFields
decomposePar -force
mpirun -np 4 compInterFoam -parallel > log.run 2>&1
```

Set `endTime`, `maxDeltaT` and `writeInterval` in `system/controlDict` yourself, following the stage table above. `numberOfSubdomains` in `system/decomposeParDict` must match the number of processes.

### 2.4 Mesh study with titanium inputs from the literature

`cases/fsLIFT_332mJ_Ti_meshStudy` is the reference case with literature inputs for titanium and a mesh family for a grid-convergence study:

| Input | Reference case | Mesh-study case | Source |
|---|---|---|---|
| Density | 4515 kg/m³ | 4506 kg/m³: solid Ti at 300 K, so that with the deposited 71.4 nm the mass per area of the film is right | room-temperature value; the liquid (4193 kg/m³ at the melting point, Ozawa et al. 2017) would need a thicker film |
| Surface tension | 1.64 N/m | 1.613 N/m, the measured value at the melting point | Seimiya et al. 2022 |
| Electron heat capacity | 630 T_e J/(m³ K) | 318 T_e J/(m³ K) | Hobmaier et al. 2026 (Sommerfeld fit) |
| Electron–lattice coupling | table to 1e20 W/(m³ K) | the same table capped at 2.5e18 W/(m³ K) | Hobmaier et al. 2026 (first-principles peak) |

The two-temperature inputs are those of case C of `tools/ttm_parameter_sensitivity.py`, so its 1-D solution can be compared with the early heating. The surface tension is constant: the measured temperature dependence holds near the melting point, and the linear fit reaches zero near 9,800 K, far outside the measured range. The viscosity (2.35 mPa s) is unchanged and not yet checked against a measurement.

`tools/make_mesh_level.py <case> <1|2|3>` writes the mesh of one level. Level 2 has the axial cells of the reference case and 35 nm radial cells for r < 3.5 µm, graded by 4 to the edge; levels 1 and 3 have √2 fewer and √2 more cells in each direction (62,010, 124,200 and 248,252 cells). The largest time steps are scaled with (cell size)^1.5, the scaling of the capillary time-step limit. The whole study runs with

```bash
tools/run_mesh_study.sh <work dir> [number of processes]
```

It runs levels 1, 2 and 3 to 1 ns, then continues levels 2 and 1 to 20 ns. On 4 processes, level 3 takes about 6 hours to 1 ns and would need about 45 hours more to reach 20 ns, so it is not continued. To analyse the results:

```bash
for L in 1 2 3; do xvfb-run -a pvpython tools/mesh_study_fields.py <work dir>/L$L; done
python3 tools/mesh_study_gci.py <work dir> 1e-9
```

`mesh_study_gci.py` computes the apparent order of convergence, the extrapolated value and the grid convergence index of each quantity by the procedure of Celik et al. (2008).

The run script needs `mpirun`; in a container that runs as root, set `OMPI_ALLOW_RUN_AS_ROOT=1` and `OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1`.

---

## 3. The one-dimensional test cases

Both are single columns of cells along y through the same layer sequence as the reference case: receiver, gap, film and gas. They run on one core.

**`cases/column1D_laser_TTM`: laser heating and the two-temperature model, 0–20 ps.** Its laser and material inputs are kept identical to those of the independent finite-difference code `tools/ttm_fd_1d.py`. They are verification inputs, not the physical titanium optics: reflectivity 0.35, absorption coefficient 1.07e7 1/m, and k_e = C_e · D_e.

```bash
cd cases/column1D_laser_TTM
blockMesh && rm -rf 0 && cp -r 0.orig 0 && topoSet && setFields && compInterFoam > log.run 2>&1
cd ../.. && python3 tools/compare_ttm_1d.py cases/column1D_laser_TTM
```

The comparison prints the peak electron temperature and the time at which the electron and lattice temperatures meet, for both codes, and draws `figures/verification/ttm_column_vs_fd.png`.

**`cases/column1D_galilean_test`: a film moving as a whole, 0–100 ps.** The laser is off, and the whole column starts at −170 m/s with gas entering at the top and leaving at the bottom. The exact solution is that nothing changes: the velocity stays −170 m/s and the pressure stays 1e5 Pa everywhere. The run is set up the same way as above.

---

## 4. The case settings explained

The dictionaries of `cases/fsLIFT_332mJ_axisymmetric` are commented. This section collects the settings that matter most. The default is the value the solver uses when the entry is missing.

### 4.1 Laser: `constant/laserProperties`

| Entry | Reference value | Meaning |
|---|---|---|
| `pulseEnergy` | 4.1089e-8 J | Energy that reaches the film, after the glass (see [MODEL.md §3](MODEL.md#3-laser-energy-deposition)) |
| `pulseWidth` | 214e-15 s | Full width at half maximum of the Gaussian pulse |
| `wavelength` | 343e-9 m | |
| `spotSize` | 6e-6 m | Beam **diameter** at 1/e² intensity (radius 3 µm) |
| `focus`, `direction` | (0 20.0357e-6 0), (0 −1 0) | Beam centre on the axis at mid-depth of the film; the beam travels in −y |
| `absorptionCoeff` | 1.03469e8 1/m | 4πk/λ for titanium at 343 nm |
| `reflectivity` | 0.42832 | Glass–titanium interface |
| `beamCoverage` | `geometric` (default) | Deposits the exact share of the energy that falls into the wedge |
| `maxVolumetricSource` | 7e24 W/m³ | Safety cap, about 9 times the peak source of this case, so it does not act |

To change the fluence, change `pulseEnergy`: energy on the film = 0.87545 × F₀ · π w₀² / 2 for this glass.

### 4.2 Two-temperature model: `twoTemperatureProperties` in `system/controlDict`

| Entry | Reference value | Meaning |
|---|---|---|
| `Ce` | 630 · T_e | Electron heat capacity [J/(m³ K)] |
| `G` | table | Electron–lattice coupling against T_e |
| `Cl`, `kl` | 2.5e6, 4.5 | Lattice heat capacity [J/(m³ K)] and conductivity [W/(m K)] |
| `ke0` | 18.9 | Electron conductivity k_e = ke0 · T_e/T_l [W/(m K)]. If `ke0` is absent the older form k_e = C_e · D_e (entry `De`) is used, which is too high in a hot film |
| `keMax` | very large (default) | Optional upper bound on k_e |
| `electronSubCycles` | 4 | Sub-steps of the electron equation per time step |
| `maxTe`, `maxTl` | 50,000 K | Safety bounds. Keep them well above the peak temperatures: 10,000 K was too low for this laser and removed energy |
| `advectTemperatures` | `true` (default) | Carry T_e and T_l with the flow |

### 4.3 Evaporation, recoil and trapped vapour: `system/controlDict`

`phaseChangeCoeffs`:

| Entry | Reference value | Meaning |
|---|---|---|
| `model` | `clausius_clapeyron` | Saturation pressure model |
| `hf`, `gasConstant`, `Tvap`, `p_ref` | 9.1e6 J/kg, 174 J/(kg K), 3560 K, 101325 Pa | Inputs of the saturation pressure |
| `evaporationCoeff` | 0.18 | β in the Hertz–Knudsen flux |
| `phaseChangeVolumeSource` | `false` (default `true`) | `false`: the evaporated mass leaves without creating an artificial vapour volume in the gas |

`advancedInterfaceCapturing`:

| Entry | Reference value | Meaning |
|---|---|---|
| `recoilCoefficient` | 0.54 | Recoil pressure / p_sat on a free surface |
| `momentumAccommodationCoeff` | 1.0 | Factor on the recoil; 1 is the physical value, since 0.54 · p_sat already contains the back-flux of vapour |
| `maxRecoilPressure` | 5e9 Pa | Bound on the recoil pressure |
| `maxPhysicalTemperature` | 10,000 K | Highest temperature used in p_sat and the evaporation flux |
| `confinedVapourFactor` | 1.0 (default) | Older option that scaled the glass-side recoil; keep 1 when the trapped-vapour model is on |
| `confinedVapourPocket { }` | see below | Vapour trapped between the glass and the film |

`confinedVapourPocket` sub-dictionary:

| Entry | Reference value | Default | Meaning |
|---|---|---|---|
| `enabled` | `true` | `false` | Switch the model on |
| `mode` | `local` | `uniform` | `local`: one vapour column per radial ring (recommended). `uniform`: one pocket for the whole spot; it pushes the film as a rigid plug |
| `binWidth` | 1e-7 m | 1e-7 m | Width of a ring (2 radial cells in the reference mesh) |
| `maxRadius` | 1e-5 m | 1e-5 m | Cells further from the axis belong to no ring and feel no vapour pressure |
| `minThickness` | 1e-9 m | 1e-9 m | Initial gap h₀ between glass and film |
| `metalDensity` | 4515 kg/m³ | 4515 kg/m³ | Converts the evaporated mass into the volume it leaves behind |
| `filmMinHeight` | 9e-6 m | no limit | Only metal above this y belongs to the donor film. Set it between the receiver top and the film |

### 4.4 Interface and pressure: `system/fvSolution`

In `solvers/"alpha.*"`:

| Entry | Reference value | Meaning |
|---|---|---|
| `nAlphaSubCycles` | 2 | Sub-cycles of the α equation per time step |
| `cAlpha` | 2.0 | Interface compression coefficient. The compression speed is min(`cAlpha`, `compressionLimiter`) × the flow speed at the face, so with the values here the effective coefficient is 0.5 and `cAlpha` has no effect above 0.5 ([MODEL.md §7.1](MODEL.md#71-volume-fraction-transport)) |
| `compressionLimiter` | 0.5 | Upper bound on the compression speed as a fraction of the face flow speed; here it sets the effective compression coefficient |
| `compressionInterfaceBand` | `true` (default) | Apply compression only near the interface; compression everywhere created small dips of α inside the metal |
| `stockCompressionFlux` | `true` (default) | Standard OpenFOAM form of the compression flux; it keeps α bounded and conserves the metal volume exactly |

In `compInterFoamCoeffs`:

| Entry | Reference value | Meaning |
|---|---|---|
| `surfaceTensionInterfaceBand` | `true` (default) | Apply surface tension only near the interface |
| `surfaceTensionBandLayers` | 3 (default) | Width of that band in cells; 1 layer missed about 17 % of the force |
| `recoilFaceInterpolation` | `localMax` (default) | Face value of the recoil pressure from the larger neighbour. `linear` applied only half the recoil |
| `pressureClamp` | `false` (default `true`) | **Keep `false`.** The clamp overwrote the pressure after the solve, so the velocity no longer satisfied continuity |
| `maxPressureGradient` | 1e300 | Limiter on the pressure gradient, switched off |
| `maxVelocity` | 0 | Speed of the metal above which the run stops with an error (default 500 m/s); 0 switches the check off |
| `momentumBudget` | `true` (default `false`) | Log the forces on the metal and the change of its momentum every time step |

### 4.5 Sub-models that must stay off

In `system/controlDict`: `enablePhaseExplosion false` and `enablePlasmaModel false`. These sub-models are not verified, and the phase-explosion mass source caused a numerical runaway with the physical laser inputs.

### 4.6 Duplicated material values

The density used in the flow equations is taken from `constant/thermophysicalProperties.metal` and `.air`. `constant/transportProperties` repeats some of the values (it is printed at start-up as `rho=4515` and `rho=1.2`). Keep the files consistent when changing a material.

---

## 5. Analysis and plotting tools

All tools take the case directory as their first argument. The Python tools need NumPy and Matplotlib; the `pvpython` tools need ParaView (run them with `xvfb-run -a` on a machine without a display).

| Tool | What it does | Output |
|---|---|---|
| `tools/axis_timeseries.py <case>` | Reads the line samples along the beam axis (`postProcessing/axisLine`) and follows the largest metal piece on the axis: position, speed, lattice temperature | `tools/axis_timeseries_<case>.csv`, `figures/<case>/data/axis_timeseries.csv`, `figures/<case>/axis_position_speed_temperature.png` |
| `tools/results_package.py <case> [comparison case]` | Results table and a four-panel figure: speed on the axis, distance moved, trapped-vapour pressure, metal volume. Needs the output of `axis_timeseries.py` | `figures/<case>/results_table.md`, `results_timeseries.png`, `data/trapped_vapour_pressure.csv` |
| `tools/detect_droplet.py <case> [--all]` | Separate metal pieces along the beam axis | printed; `droplet_detections.txt` |
| `pvpython tools/count_metal_pieces.py <case> <times>` | All separate metal pieces in the whole domain at the given times: volume, position, mean velocity | printed |
| `pvpython tools/render_dome.py <case> [t1,t2,…]` | Sections through the axis coloured by vertical velocity and by temperature, with the α = 0.5 interface | `figures/<case>/dome_*.png` |
| `pvpython tools/render_gallery.py <case>` | Sections of α, speed with direction arrows, radial velocity and pressure at several times; close-ups of the laser heating; 3-D views of the revolved interface | `figures/<case>/gallery/` |
| `python3 tools/compare_ttm_1d.py <case>` | The 1-D column against the finite-difference code | `figures/verification/ttm_column_vs_fd.png` |
| `python3 tools/ttm_fd_1d.py` | Independent finite-difference solution of the 1-D two-temperature model, with its own convergence check | printed |

For a decomposed case, `count_metal_pieces.py` merges the processor pieces before counting. Without the merge, the shared points on processor boundaries make one piece look like several.

The drop benchmarks are described in [VERIFICATION.md](VERIFICATION.md#3-how-to-repeat-the-tests).

---

## 6. Setting up a new case

1. **Copy `cases/fsLIFT_332mJ_axisymmetric`.** It holds every correction made during testing.
2. **Keep these settings unless you have tested the change:**
   - `zeroGradient` for `alpha.metal` and `alpha.air` on the outer boundary `right`. The real film continues far beyond the domain. An inflow condition with α = 0 cut the film open at the boundary, and the surface tension of that false edge added a spurious force on the film;
   - `maxTe`, `maxTl` and `maxT` well above the expected peak (50,000 K here);
   - `ke0` for the electron conductivity;
   - `pressureClamp false`;
   - trapped vapour in `mode local`;
   - `enablePhaseExplosion false` and `enablePlasmaModel false`.
3. **Fluence:** change `pulseEnergy` (see §4.1). For a different glass, recompute the transmitted fraction.
4. **Mesh:** along the path of the dome tip, keep the cells at about 6 nm or smaller for as far as the dome travels. On coarser cells the tip smears out and the interface reaches a plateau of α ≈ 0.5. The reference mesh has 6 nm cells for the first 3 µm below the film and 50 nm radial cells.
5. **Time step:** raise `maxDeltaT` above 1 ps only after the laser and the fastest heating are over (about 1 ns), and keep it below the capillary limit of the mesh.
6. **Disk space:** about 5 GB for the reference case with the default write intervals.
