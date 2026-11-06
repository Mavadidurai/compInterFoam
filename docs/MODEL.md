# The model

This document describes what the solver computes: the physical sequence it represents, the equations, the input values of the reference case and the order of operations in one time step. The values quoted are those of `cases/fsLIFT_332mJ_axisymmetric`.

Symbols used throughout:

| Symbol | Meaning | Unit |
|---|---|---|
| α | metal volume fraction of a cell: 1 = metal, 0 = gas | – |
| T_e | electron temperature | K |
| T_l | lattice temperature: the temperature of the atoms | K |
| p_sat | saturation (vapour) pressure of the metal at the surface temperature | Pa |
| j | evaporation mass flux | kg/(m² s) |
| r, y | distance from the beam axis, and position along it | m |

---

## 1. The physical sequence and its time scales

The laser enters through the glass and travels in −y. It is absorbed in the titanium film just below the glass. The film surface facing the glass is called the *glass-side surface* below, and the surface facing the gap the *gap-side surface*.

| Stage | Time after the pulse | What happens | Model part |
|---|---|---|---|
| 1. Absorption | 0–0.5 ps | The 214 fs pulse is absorbed by the electrons within ~10 nm of the glass-side surface | laser source (§2) |
| 2. Electron–lattice heating | ~0.5–5 ps | The hot electrons heat the atoms and carry heat deeper into the film | two-temperature model (§3) |
| 3. Melting and evaporation | ps to ns | Hot metal evaporates, first at the glass-side surface | evaporation (§4) |
| 4. Pressure build-up | ps to ns | On the gap side the vapour escapes freely and pushes back on the surface (recoil). On the glass side it is trapped and builds a pressure | recoil (§4) and trapped vapour (§5) |
| 5. Film motion | 0.1 to tens of ns | The pressure difference pushes the molten centre of the film towards the receiver and forms a dome | flow and interface (§6) |
| 6. Ejection, or not | ns to tens of ns | The dome thins. It either separates into a droplet or is pulled back by surface tension | flow and surface tension (§6) |

The problem is numerically demanding because of these ranges:
- **length:** 10 nm (absorption depth) to 71.4 nm (film) to 12 µm (gap);
- **time:** femtoseconds (laser) to picoseconds (electrons) to tens of nanoseconds (ejection).

---

## 2. Laser heat source (`femtosecondLaserModel.C`)

Absorbed power per unit volume:

  Q = (1 − R) · I(t) · exp(−2r²/w₀²) · a · exp(−a·d)

where:
- I(t) is a Gaussian pulse in time, with the full width at half maximum given by `pulseWidth` (214 fs);
- w₀ = 3 µm is the beam radius at which the intensity falls to 1/e²;
- d is the depth below the glass-side surface;
- a is the absorption coefficient and R the reflectivity of the glass–titanium interface.

**Optical inputs (titanium at 343 nm, refractive index n = 1.7283, k = 2.8242):**
- a = 4πk/λ = 1.0347e8 1/m, i.e. an absorption depth of 9.66 nm;
- R = 0.428 at the glass–titanium interface (glass index 1.5463).

**Energy reaching the film.** The nominal peak fluence 332 mJ/cm² is defined before the 1.05 mm glass plate. The corresponding pulse energy is E = F₀ · πw₀²/2 = 3320 J/m² · π · (3e-6 m)² / 2 = 4.694e-8 J. The air–glass surface reflects 4.60 % and 91.77 % of the rest crosses the glass, so the film receives (1 − 0.0460) · 0.9177 = 0.87545 of it: 4.109e-8 J. The derivation is repeated in the comments of `constant/laserProperties`.

**Wedge geometry.** The computational domain is a 5° slice of the cylinder. `beamCoverage geometric` deposits exactly the share of the pulse energy that falls into the slice.

---

## 3. Two-temperature model (`twoTemperatureModel.C`, `TEqn.H`)

A femtosecond pulse heats the electrons far faster than they can pass the energy to the atoms. The model therefore carries two temperatures:

  C_e(T_e) ∂T_e/∂t = ∇·(k_e ∇T_e) − G(T_e) (T_e − T_l) + Q

  C_l ∂T_l/∂t = ∇·(k_l ∇T_l) + G(T_e) (T_e − T_l) − (latent heat removed by evaporation at the surface)

| Input | Value in the reference case |
|---|---|
| Electron heat capacity C_e | 630 · T_e J/(m³ K) |
| Electron–lattice coupling G(T_e) | table from 1e18 W/(m³ K) at 300 K to 1e20 W/(m³ K) at 20,000 K |
| Lattice heat capacity C_l | 2.5e6 J/(m³ K) |
| Lattice conductivity k_l | 4.5 W/(m K) |
| Electron conductivity k_e | k_e0 · T_e / T_l with k_e0 = 18.9 W/(m K) (option `ke0`) |

**Electron conductivity.** The form k_e = k_e0 · T_e/T_l is the usual Drude-type form for T_e well below the Fermi temperature. The older form k_e = C_e · D_e (still available) grows with T_e even when the electrons and atoms are at the same temperature. In a hot film this made it 10–17 times too high, and the film moved the wrong way in early tests. Use `ke0`.

**Time stepping.** The electron equation is stiff and is sub-stepped within each flow time step (`electronSubCycles 4`).

The temperatures are also carried with the flow (`advectTemperatures true`), so ejected metal keeps its heat.

---

## 4. Evaporation and recoil (`twoPhaseMixtureThermo.C`, `advancedInterfaceCapturing.C`)

**Saturation pressure** (Clausius–Clapeyron):

  p_sat = p_ref · exp[ (h_v / R_v) · (1/T_v − 1/T_l) ]

with p_ref = 101,325 Pa at T_v = 3560 K, latent heat of vaporisation h_v = 9.1e6 J/kg and gas constant of titanium vapour R_v = 174 J/(kg K). The temperature used in p_sat and j is capped at `maxPhysicalTemperature` = 10,000 K.

**Evaporation flux** (Hertz–Knudsen):

  j = β · p_sat / √(2π R_v T_l),  with β = 0.18 (`evaporationCoeff`).

The interface area used with j is the area of the cell faces projected on the interface normal. The evaporated mass leaves the metal; no artificial vapour volume is added to the gas (`phaseChangeVolumeSource false`).

**Recoil pressure on a free surface** (Anisimov's result for the vapour leaving through a thin non-equilibrium layer): p_r = 0.54 · p_sat. It is applied to the fluid as a force on the cell faces at the interface. The face value is taken from the larger neighbour (`recoilFaceInterpolation localMax`). A plain average applied only half the recoil and was corrected during testing. The recoil is bounded at 5 GPa (`maxRecoilPressure`), a value the reference run does not reach.

---

## 5. Vapour trapped between the glass and the film

**Why it is needed.** In the experiment the glass-side surface of the film is bonded to the glass. Vapour formed there cannot escape: it is trapped and pushes the film away from the glass. Both phases in the flow solver are incompressible (see [Limitations](RESULTS.md#limitations)), so the field equations cannot hold a compressible vapour volume. The vapour is therefore held in a small model of its own, solved alongside the flow.

**One vapour column per radial ring** (`mode local`). The glass-side surface is divided into rings 100 nm wide (`binWidth`). In each ring b, at every time step:

```
m_b += Δt · Σ j·dA                  evaporated mass from the glass-side interface cells of the ring
d_b += Δt · (mean film velocity away from the glass)
V_b  = A_b · (h₀ + max(d_b, 0)) + m_b / ρ_metal     vapour volume: initial gap h₀ = 1 nm, the opened gap,
                                                     and the volume of the metal that evaporated
p_b  = min( p_sat,b ,  m_b · R_v · T_v,b / V_b )    ideal gas, at most the saturation pressure
```

- If the ideal-gas pressure would exceed the saturation pressure, the excess vapour condenses (m_b is reduced).
- p_sat,b and the vapour temperature T_v,b are averages over the ring, weighted by the evaporation flux, so the hottest surface dominates.
- The pressure applied to the glass-side surface is the larger of p_b and the free-surface recoil, through the same face force as in §4.
- Only metal above `filmMinHeight` = 9 µm counts as donor film; the receiver top is at 8 µm.

**Why one column per ring, and not one pocket for the whole spot.** A single pocket for the whole spot (`mode uniform`) pushed the whole 7.4 µm-wide heated area as a rigid plug and tore the film at its cold edge. Within a few nanoseconds the real vapour cannot flow far sideways, so independent rings are the better approximation.

**Restarts.** The vapour state is written to `<time>/uniform/vapourPocket`. A restart with a different `mode` or `binWidth` stops with an error rather than continue with an inconsistent state.

**Simplifications:**
- no sideways vapour flow between rings;
- no heat loss of the vapour to the glass;
- the vapour temperature equals the flux-weighted surface temperature;
- the initial gap h₀ is not known; 1 nm is assumed.

---

## 6. Flow and interface (`alphaEqn.H`, `UEqn.H`, `pEqn.H`)

**Interface capturing.** The volume-of-fluid method with OpenFOAM's MULES scheme transports α. An interface compression term keeps the interface sharp:
- `cAlpha 2`, applied only near the interface (`compressionInterfaceBand true`);
- `compressionLimiter 0.5`;
- 2 sub-cycles of the α equation per time step.

The compression uses the interface normal of the current time step. An earlier version used a normal frozen at the start of the run, so the compression had no effect on a moving interface; this was found and corrected during testing (see [VERIFICATION.md](VERIFICATION.md#4-defects-found-and-corrected)).

**Phases.** Both phases are incompressible and Newtonian:

| Phase | Density | Viscosity |
|---|---|---|
| Metal (titanium, liquid everywhere) | 4515 kg/m³ | 2.35e-3 Pa s |
| Gas (in the gap, and standing in for the glass) | 1.16 kg/m³ | 1.8e-5 Pa s |

**Surface tension.** Continuum surface force with σ = 1.64 N/m. It is applied in a band three cells wide around the interface (`surfaceTensionBandLayers 3`); a band one cell wide missed about 17 % of the force.

**Pressure–velocity coupling.** PIMPLE, with 3 outer and 3 pressure correctors. The momentum equation contains the recoil and trapped-vapour face forces and the surface tension.

---

## 7. What happens in one time step

Order of operations in `compInterFoam.C`:

1. **Time step.** Δt is the smallest of the limits from the flow Courant number (`maxCo 0.2`), the interface Courant number (`maxAlphaCo 0.5`), the thermal Courant number (`maxThermalCourant 0.5`) and `maxDeltaT`.
   - `tools/run_case.sh` sets `maxDeltaT` to 10 fs for the first 2 ps (pulse and electron–lattice coupling), 1 ps up to 1 ns, and 2.5 ps after that.
   - The capillary time-step limit on this mesh is about 4.4 ps.
2. **Laser:** the source Q is updated.
3. **PIMPLE outer loop:**
   1. mixture properties and phase-change sources are updated;
   2. **α transport** (2 sub-cycles). In each sub-cycle the recoil and trapped-vapour state is refreshed, then α is advected with interface compression and bounded to [0, 1]. The density is updated **once, after** the sub-cycles. Updating it inside them gave a moving interface a false pressure of −ρU²; this was found by the Galilean test and corrected;
   3. **temperatures:** T_e and T_l (two-temperature model with electron sub-steps);
   4. **recoil and trapped-vapour pressure** (§4–5);
   5. **momentum predictor**, with the recoil force and surface tension;
   6. **pressure correctors**, re-evaluating the recoil in each one. The optional pressure clamp must be off (`pressureClamp false`): it overwrote the pressure after the solve, so the velocity no longer satisfied continuity.
4. **Diagnostics and output:**
   - energy audit, and the force budget on the metal when `momentumBudget` is on;
   - a droplet-breakup indicator, for information only (it changes α only if `applyAlphaRedistribution` is set, which no case does);
   - writing of the fields.

---

## 8. Sub-models that are present but switched off

The source also contains a phase-explosion model, a plasma model and a breakup model (`enhancedLIFTPhysics.C`). They have not been verified. With the physical laser inputs, the mass source of the phase-explosion model caused a numerical runaway within a few picoseconds. They are switched off in every provided case (`enablePhaseExplosion false`, `enablePlasmaModel false`); the trapped-vapour model represents the vaporisation instead.
