# Femtosecond Laser Energy Deposition Model for LIFT Simulations
## Detailed Implementation and Physical Logic

---

## 1. Overview and Physical Motivation

The femtosecond laser model developed for this Laser-Induced Forward Transfer (LIFT) simulation addresses the unique challenges of modeling ultrashort pulse laser-matter interactions at the femtosecond timescale. Unlike nanosecond or picosecond laser systems, femtosecond lasers deposit energy faster than thermal diffusion can occur, creating highly non-equilibrium conditions that require specialized treatment.

The fundamental physical challenge in femtosecond LIFT is that the laser pulse duration (typically 100-500 fs) is significantly shorter than the electron-phonon relaxation time (several picoseconds for metals like titanium). This temporal separation means that electrons absorb the laser energy first and subsequently transfer this energy to the lattice through electron-phonon coupling, which is the basis for the two-temperature model approach used in the broader simulation framework.

---

## 2. Core Physical Parameters and Their Justification

### 2.1 Laser Pulse Characteristics

The model captures four fundamental laser parameters that define the energy deposition profile:

**Peak Intensity** represents the maximum power per unit area delivered by the laser beam. For femtosecond LIFT processes, this typically ranges from 10^11 to 10^13 W/m^2, as documented in experimental studies by Piqué et al. (2002) and Brown et al. (2010). The peak intensity directly determines the amount of energy deposited into the electron subsystem during the brief interaction time.

**Pulse Width** defines the temporal duration of the laser pulse. For femtosecond LIFT, pulse widths between 100 and 500 femtoseconds are typical. This parameter is critical because it determines whether the heating process is truly non-equilibrium. At these timescales, the pulse ends before significant electron-phonon thermalization occurs, which fundamentally differs from longer pulse regimes.

**Pulse Energy** represents the total energy contained in a single laser pulse. This is typically specified in microjoules or millijoules and is related to peak intensity through the spatial and temporal pulse profile. The model enforces energy conservation by verifying that the integrated deposited energy matches the specified pulse energy, accounting for reflection and absorption losses.

**Wavelength** affects the optical properties of the interaction, particularly absorption depth through the complex refractive index. Common femtosecond laser wavelengths include 800 nm (Ti:Sapphire fundamental) and 400 nm (second harmonic). The absorption coefficient, which determines the penetration depth, is wavelength-dependent and must be specified accordingly.

### 2.2 Spatial Distribution Parameters

**Spot Size** defines the beam diameter at the focal plane. For LIFT applications, spot sizes typically range from 10 to 100 micrometers. The spatial energy distribution follows either a uniform (top-hat) or Gaussian profile, with Gaussian profiles being more physically representative of real laser beams. The Gaussian intensity profile is implemented as:

I(r) = I_peak × exp(-2r²/w₀²)

where w₀ is the beam waist (related to spot size) and r is the radial distance from the beam center.

**Focus Position and Direction** specify where the beam is targeted within the computational domain and its propagation direction. For standard LIFT geometry, the laser typically propagates perpendicular to the film plane (z-direction), with the focus positioned at or near the film surface. The model includes geometric validation to ensure the focus lies within the computational mesh and within the film region for maximum energy coupling.

**Absorption Coefficient** quantifies how rapidly laser energy is absorbed as the beam penetrates the material. For titanium films at 800 nm wavelength, the absorption coefficient is approximately 5×10^6 m^-1, corresponding to a penetration depth (1/α) of about 200 nm. This value is derived from optical property measurements by Wellershoff et al. (2000) in their comprehensive study of titanium's electronic and thermal properties. The implementation allows separate specification of absorption coefficients for the metal film and surrounding gas, recognizing that most of the absorption occurs in the metal.

### 2.3 Optical Interface Properties

**Reflectivity** accounts for Fresnel reflection at the film interface. For titanium at normal incidence and 800 nm wavelength, reflectivity is approximately 50-60% based on the complex refractive index. However, the effective reflectivity can be lower in LIFT due to surface roughness and oxide layers. The model implements reflectivity as a reduction in transmitted intensity, converting the specified reflectivity into a transmission coefficient:

T = 1 - R

where T is transmission and R is reflectivity. This approach correctly reduces the incident energy available for absorption.

**Incidence Angle** modifies the reflectivity according to Fresnel equations. For normal incidence (0°), reflectivity is minimum for most materials. As the angle increases, both s- and p-polarized components experience different reflection behavior, though for the small angles typical in LIFT (near-normal incidence), this effect is relatively minor.

---

## 3. Temporal Envelope Implementation

The temporal behavior of the laser pulse is one of the most critical aspects of the model, as it determines when and how rapidly energy is deposited. The implementation handles both pulsed and continuous laser operation modes.

### 3.1 Pulsed Operation Logic

For pulsed operation, the temporal envelope follows a Gaussian profile centered at the pulse center time:

E(t) = E_pulse × exp(-(t - t_center)²/(2σ²))

where σ = pulse_width/(2√(2ln2)) converts the full-width at half-maximum (FWHM) pulse width into the Gaussian standard deviation.

The key challenge in implementation is that the simulation time step (typically 1-10 fs) may not align with the pulse temporal structure. Therefore, the model integrates the Gaussian envelope over each time step to compute the energy fraction delivered during that step. This integration is performed analytically using the error function:

Integral from t₁ to t₂ of Gaussian = (√(π/2) × σ) × [erf((t₂-t_c)/(√2σ)) - erf((t₁-t_c)/(√2σ))]

This approach ensures energy conservation regardless of time step size, which is essential for accurate energy accounting.

### 3.2 Continuous Wave Modification

For continuous laser operation (used in some LIFT variants or for model validation), the temporal envelope is simply unity within the specified active time window [laserStartTime, laserEndTime]. The pulse duty cycle parameter allows modeling of pulsed operation at high repetition rates where pulses overlap temporally.

### 3.3 Energy Conservation Verification

A critical feature is the pulse-integrated energy verification system. The model tracks:

1. **Expected Energy**: Computed from the temporal envelope integral
2. **Deposited Energy**: Actual energy added to the electron subsystem
3. **Incident Energy**: Total energy impinging on the domain
4. **Film vs. Gas Energy**: Separate accounting for film and surrounding gas absorption

After each complete pulse, the model compares deposited energy to expected energy, issuing warnings if the discrepancy exceeds specified tolerances (typically 1% relative or 1 pJ absolute). This verification catches numerical errors, mesh resolution issues, or incorrect parameter specifications.

---

## 4. Spatial Distribution and Beer-Lambert Absorption

The spatial energy deposition follows the Beer-Lambert law for absorption in the propagation direction, combined with the lateral beam profile.

### 4.1 Lateral Profile (Gaussian)

For Gaussian profiles, the intensity at radial position r from the beam axis is:

I(r) = I₀ × exp(-2r²/w₀²)

where w₀ = spotSize/2 is the beam waist radius. The factor of 2 in the exponent ensures that the intensity drops to 1/e² (≈13.5%) at r = w₀, which is the standard definition for Gaussian beam radius.

### 4.2 Axial Absorption (Beer-Lambert Law)

As the beam propagates through the material, intensity decreases exponentially according to:

I(z) = I₀ × exp(-α × z)

where α is the absorption coefficient and z is the propagation distance into the material. For a thin film of thickness d, the energy deposited per unit volume at depth z is:

Q(z) = α × I(z) = α × I₀ × exp(-α × z)

This exponential decay means that most energy is deposited near the illuminated surface. For titanium with α = 5×10^6 m^-1, 86% of the transmitted energy is absorbed within the first 40 nm.

### 4.3 Implementation Strategy

The implementation evaluates this spatial profile for each computational cell by:

1. **Determining if the cell is in the beam path** by checking if the cell center lies within the beam cone defined by the focus point, propagation direction, and beam divergence.

2. **Computing radial distance** from the cell center to the beam axis by projecting the position vector onto the plane perpendicular to the beam direction.

3. **Computing axial distance** along the beam direction, starting from the film surface or a reference plane.

4. **Applying both Gaussian lateral and exponential axial profiles** to determine the local volumetric energy deposition rate.

5. **Accounting for interface transmission** by reducing the incident intensity according to the transmission coefficient before applying absorption.

The implementation includes a critical geometric consideration: it distinguishes between cells in the metal film (high absorption) and cells in the surrounding gas (low absorption) using the film Y-coordinate bounds. This allows accurate modeling of the selective heating of the metal film while avoiding spurious heating of the gas phase.

---

## 5. Film Geometry and Coordinate System

The model assumes a Cartesian coordinate system where the thin film is oriented parallel to the x-z plane, with its thickness defined in the y-direction. This geometry is standard for LIFT simulations and matches the experimental configuration.

### 5.1 Film Bounds Specification

The film is defined by two y-coordinates:
- **filmYMin**: Bottom surface of the film
- **filmYMax**: Top surface of the film (irradiated surface)

These bounds are specified either directly or computed from the film thickness. For example, for a 71.4 nm titanium film positioned at y = 0:
- filmYMin = -71.4e-9 m
- filmYMax = 0 m

This positioning is critical because the absorption calculations reference these bounds to determine whether cells should use the metal or gas absorption coefficient.

### 5.2 Laser Direction Convention

The implementation requires specifying the laser propagation direction as a unit vector. For standard LIFT where the laser illuminates from above (positive y), the direction would be:

direction = (0, -1, 0)

indicating propagation in the negative y-direction toward the film. The model automatically normalizes the provided direction vector to ensure it is a unit vector.

---

## 6. Volumetric Source Term Calculation

The volumetric source term Q_laser [W/m³] represents the power deposited per unit volume in the electron subsystem. This term appears in the electron energy equation of the two-temperature model:

C_e × ∂T_e/∂t = ∇·(k_e × ∇T_e) - G(T_e - T_l) + Q_laser

where C_e is electron heat capacity, T_e is electron temperature, k_e is electron thermal conductivity, G is electron-phonon coupling, and T_l is lattice temperature.

### 6.1 Source Term Assembly Process

The source term is assembled through a multi-step process that ensures both physical accuracy and energy conservation:

**Step 1: Temporal Evaluation**
Compute the temporal envelope factor for the current time step by integrating the pulse shape over the time interval [t, t+dt]. This yields the fraction of total pulse energy delivered during this specific time step.

**Step 2: Spatial Weighting**
For each cell in the computational domain:
- Check if the cell lies within the beam path
- Compute the radial distance from the beam axis
- Compute the axial distance along the propagation direction
- Evaluate the Gaussian lateral profile at the radial position
- Evaluate the exponential absorption at the axial position
- Determine if the cell is in the film or gas region
- Apply the appropriate absorption coefficient

**Step 3: Normalization and Scaling**
The raw spatial weights are scaled to ensure that the volume-integrated source term equals the expected energy deposition based on:
- Pulse energy specified in input
- Reflectivity losses at the interface  
- Absorption efficiency based on film thickness and absorption coefficient
- Beam coverage (fraction of pulse energy intersecting the computational domain)

The scaling factor is:

scale = (Incident Power × Transmission × Absorption Efficiency × Coverage) / (Raw Integrated Power)

**Step 4: Power Capping (Optional)**
If specified, the maximum volumetric source is clamped to maxVolumetricSource to prevent numerical instabilities. This is particularly important near the irradiated surface where Beer-Lambert absorption creates very high local energy densities. However, capping must be applied carefully to maintain energy conservation; the model adjusts the scaling factor to account for capped regions.

### 6.2 Physical Justification for Source Magnitude

For typical femtosecond LIFT parameters:
- Pulse energy: 1 μJ
- Pulse width: 200 fs  
- Spot size: 20 μm diameter
- Film thickness: 71.4 nm
- Absorption coefficient: 5×10^6 m^-1

The peak volumetric source can be estimated:

Q_peak ≈ (E_pulse / (pulse_width × beam_volume)) × absorption_factor

For the above parameters, this yields approximately 10^18 to 10^20 W/m³ near the irradiated surface. These extremely high values are physically realistic for femtosecond pulses and are necessary to achieve the rapid electron heating (to temperatures exceeding 10,000 K) observed in experiments.

---

## 7. Energy Accounting and Diagnostics

Robust energy tracking is essential for verifying model correctness and diagnosing issues during simulation development.

### 7.1 Multi-Level Energy Tracking

The model maintains four cumulative energy counters:

**Cumulative Incident Energy**: Total energy that would have been deposited if there were no reflection losses. This represents the energy carried by the beam intersecting the computational domain.

**Cumulative Absorbed Energy**: Total energy actually deposited into the electron subsystem after accounting for reflection and the spatial energy distribution. This is the volume integral of Q_laser over all time steps.

**Cumulative Film Energy**: Portion of absorbed energy deposited specifically in the metal film region (filmYMin ≤ y ≤ filmYMax).

**Cumulative Gas Energy**: Portion of absorbed energy deposited in the surrounding gas phase. This should be zero or negligible if gas absorption is disabled.

### 7.2 Pulse-Level Energy Verification

For pulsed operation, the model accumulates energy within each individual pulse and compares it to the expected pulse energy at pulse completion. This per-pulse verification is more stringent than total energy conservation and catches subtle errors in the temporal envelope integration.

The verification logic:

1. Detect pulse start when temporal envelope becomes significant
2. Accumulate deposited energy while envelope remains active  
3. Accumulate expected energy from temporal envelope integral
4. When envelope returns to near-zero, compare accumulated values
5. Issue warning if discrepancy exceeds tolerance

Typical tolerance values are 1% relative error or 1 pJ absolute error, with the more lenient criterion applying. This accounts for both small pulse energies (where absolute error matters) and large pulses (where relative error is more meaningful).

### 7.3 Diagnostic Output

At regular intervals (typically every 10 time steps), the model outputs detailed diagnostics including:

- Current time and time step
- Temporal envelope value and integrated energy this step
- Number of cells in beam path, in film, and in gas
- Maximum volumetric source value [W/m³]
- Total power being deposited [W]
- Cumulative energies for all categories [J]
- Pulse counter and pulse completion status
- Geometry validation (focus position relative to mesh and film bounds)

These diagnostics are invaluable during simulation setup and debugging, allowing rapid identification of issues such as:
- Focus positioned outside the computational domain
- Focus outside the film region (leading to zero film heating)
- Zero absorption (incorrect absorption coefficient)
- Energy conservation violations (mesh resolution too coarse)

---

## 8. Advanced Features and Extensions

### 8.1 Laser Scanning Capability

The model includes provisions for translating the laser focus during the simulation to model scanning LIFT processes. The focus position is updated according to:

focus(t) = initialFocus + scanVelocity × t

where scanVelocity is specified as a vector [m/s]. This allows modeling of applications where the laser beam moves across the donor film to transfer multiple pixels or patterns.

### 8.2 Multi-Pulse Operation

For applications involving multiple laser pulses (burst mode or high-repetition-rate systems), the model tracks individual pulses using a pulse counter and maintains per-pulse energy verification. The pulseFrequency parameter defines the time between successive pulses, and the model resets its pulse-specific accumulators at each new pulse start.

### 8.3 Separate Gas and Film Absorption

Recognizing that energy deposition in the surrounding gas is physically different from film absorption, the model allows specifying separate absorption coefficients:

- **absorptionCoeff**: For the metal film (typically 10^6 to 10^7 m^-1)
- **gasAbsorptionCoeff**: For surrounding gas (typically 0 or very small)

Cell-by-cell application of the appropriate coefficient ensures that the heating is predominantly confined to the film, which is essential for accurate LIFT modeling. Spurious gas heating would artificially increase recoil pressure and distort the jet formation physics.

### 8.4 Power Limiting for Numerical Stability

The maxVolumetricSource parameter provides an optional upper bound on the volumetric source term. This is necessary because the extremely localized heating near the illuminated surface (due to Beer-Lambert absorption) can create very stiff source terms that challenge solver stability.

When limiting is active, the model:
1. Identifies cells where the raw source exceeds the limit
2. Clamps those cells to the maximum value
3. Adjusts the global scaling factor to maintain energy conservation
4. Reports the number of limited cells in diagnostics

This approach prevents numerical instability while preserving the total energy budget, though it does slightly alter the spatial distribution (smoothing out the sharpest peaks).

---

## 9. Validation Strategy and Physical Benchmarks

### 9.1 Energy Conservation Checks

The primary validation is energy conservation verification. For a single pulse:

Energy_deposited ≈ Energy_pulse × (1 - Reflectivity) × (1 - exp(-α × thickness))

The exponential term represents the fraction of transmitted energy absorbed by a film of thickness d. For titanium with α = 5×10^6 m^-1 and d = 71.4 nm:

Absorption fraction = 1 - exp(-5×10^6 × 71.4×10^-9) ≈ 0.299 (29.9%)

This prediction can be compared against the cumulative absorbed energy reported by the model.

### 9.2 Spatial Distribution Validation  

The spatial distribution can be validated by:

1. **Visual inspection** of the Q_laser field in ParaView, confirming:
   - Gaussian lateral profile with correct beam diameter
   - Exponential decay in propagation direction
   - Localization to film region
   - Peak values at the correct magnitude

2. **Profile extraction** along beam axis and radial directions, comparing to analytical Gaussian and exponential functions

3. **Cell count verification**: The number of cells in the beam should scale approximately as (π × r_beam² × penetration_depth) / cell_volume

### 9.3 Temporal Profile Validation

For Gaussian pulses, the temporal profile can be validated by:

1. Plotting energy deposition vs. time and fitting to Gaussian
2. Verifying FWHM matches specified pulse width
3. Confirming pulse center occurs at specified time
4. Checking that essentially all energy (>99%) is deposited within ±3σ of pulse center

### 9.4 Comparison to Experimental LIFT Literature

The ultimate validation is comparison of simulation outcomes to experimental LIFT results. Key benchmarks include:

**Jet Velocity**: For titanium LIFT with 1 μJ pulse energy, 200 fs pulse width, and 20 μm spot size, experimental jet velocities are typically 30-100 m/s as reported by:
- Feinaeugle et al. (2012): 50-80 m/s for similar conditions
- Brown et al. (2010): 40-100 m/s depending on fluence
- Piqué et al. (2002): 30-70 m/s for titanium LIFT

The simulation should reproduce these velocity magnitudes when all physics models (two-temperature, phase change, recoil pressure) are properly calibrated.

**Electron Temperature**: Femtosecond laser heating drives electron temperatures to 5,000-15,000 K while lattice remains near ambient during the pulse. This is consistent with two-temperature model predictions by Wellershoff et al. (2000) and Lin et al. (2008) for titanium.

**Material Removal Threshold**: Experimental ablation thresholds for titanium at femtosecond pulses are approximately 0.1-0.3 J/cm² (Chichkov et al., 1996). The simulation should exhibit material ejection above this fluence range and remain stable below it.

---

## 10. Implementation Best Practices and Lessons Learned

### 10.1 Coordinate System Consistency

One of the most common sources of error is inconsistency between the specified laser direction and the actual film orientation in the mesh. The model assumes:

- Film is perpendicular to y-axis
- Laser propagates along y-axis (either positive or negative y)
- For illumination from above, direction = (0, -1, 0)

Reversing the sign or using a different axis will cause the geometric checks to fail, placing the focus outside the film region and resulting in zero heating.

### 10.2 Parameter Scaling and Units

Femtosecond laser simulations span many orders of magnitude:
- Time: femtoseconds (10^-15 s)
- Length: nanometers (10^-9 m)  
- Power density: 10^18 - 10^20 W/m³
- Temperature: 10,000+ K

Careful attention to dimensional consistency is essential. All parameters must be specified in SI units (meters, seconds, Joules, Watts) within the OpenFOAM framework.

### 10.3 Mesh Resolution Requirements

The spatial resolution must adequately resolve:

1. **Film thickness**: At least 3-5 cells across the film thickness (71.4 nm) requires cell sizes ≤15-20 nm in the y-direction.

2. **Absorption depth**: At least 2-3 cells within one penetration depth (1/α ≈ 200 nm for titanium) to capture the exponential decay profile.

3. **Beam diameter**: At least 8-10 cells across the beam diameter to resolve the Gaussian lateral profile.

4. **Jet resolution**: For tracking the ejected material jet, cells should be fine enough to resolve the jet radius and capture interface curvature.

These requirements often necessitate adaptive mesh refinement or structured grids with local refinement in the film and beam regions.

### 10.4 Time Step Considerations

For femtosecond LIFT, the time step must satisfy:

1. **Pulse resolution**: Multiple time steps (>10) within the pulse width to accurately integrate the temporal envelope.

2. **CFL stability**: Acoustic CFL condition Δt < Δx / c_sound, where c_sound ≈ 5000 m/s for titanium.

3. **Thermal diffusion**: Explicit thermal diffusion requires Δt < Δx² / (2α_thermal), though femtosecond heating is often transport-limited rather than diffusion-limited.

4. **Phase change**: During melting/vaporization, smaller time steps may be needed to resolve the rapid phase transitions.

Typical time steps for femtosecond LIFT are 0.1-1.0 fs during the pulse and can increase to 1-10 fs afterward.

---

## 11. Coupling to Two-Temperature Model

The laser source term Q_laser computed by this model feeds into the electron energy equation:

C_e(T_e) × ∂T_e/∂t = ∇·(k_e(T_e) × ∇T_e) - G(T_e - T_l) + Q_laser

The coupling is explicit: the laser model computes Q_laser based on current time and geometry, independent of the instantaneous electron or lattice temperatures. This is physically justified because the laser energy deposition is an external forcing that depends only on the laser parameters and material optical properties, not on the temperature state.

However, there are two temperature-dependent effects that could be included in more sophisticated implementations:

### 11.1 Temperature-Dependent Reflectivity

The Fresnel reflectivity depends on the complex refractive index, which varies with electron temperature. At high electron temperatures (>10,000 K), the plasma frequency shifts and reflectivity can change by 10-20%. However, most of the laser energy is deposited during the first few tens of femtoseconds when temperatures are still moderate, so this effect is typically second-order.

### 11.2 Ionization Effects

At extreme intensities (>10^13 W/m²), multiphoton ionization and avalanche ionization can modify the absorption. For typical LIFT fluences (~0.1-1 J/cm²), these effects are minimal for metals, which already have free electrons in the conduction band.

For the current implementation, the assumption of temperature-independent optical properties is well-justified and consistent with the literature approach.

---

## 12. Debugging and Troubleshooting Guidelines

### 12.1 Zero Heating Despite Active Laser

**Symptoms**: Diagnostics show laser is active, but no temperature rise occurs.

**Diagnosis**:
1. Check focus position relative to mesh bounds
2. Verify laser direction matches film orientation  
3. Confirm absorption coefficient is non-zero
4. Check that at least some cells are identified as "in film"
5. Verify film bounds encompass the focus

**Solution**: Correct the geometric configuration so focus lies within film region and direction points into the film.

### 12.2 Energy Conservation Violations

**Symptoms**: Deposited energy significantly differs from expected pulse energy (>5% error).

**Diagnosis**:
1. Check mesh resolution - may be too coarse to integrate source accurately
2. Verify time step is small enough to resolve pulse temporal profile
3. Look for cells being capped by maxVolumetricSource
4. Check if beam extends outside computational domain (energy loss)

**Solution**: Refine mesh in beam region, reduce time step, or adjust maxVolumetricSource limit.

### 12.3 Numerical Instability

**Symptoms**: Solver diverges or crashes during laser pulse.

**Diagnosis**:
1. Check maximum source value - may exceed stable limits
2. Verify time step satisfies CFL condition with post-heating velocities
3. Look for extremely thin cells in mesh that concentrate source
4. Check if electron heat capacity becomes too small at high T_e

**Solution**: Enable maxVolumetricSource capping, reduce time step, improve mesh quality, or review electron heat capacity formulation.

---

## 13. Summary and Key Takeaways

The femtosecond laser energy deposition model implements a physically-based approach to modeling ultrashort pulse laser-matter interactions in the LIFT process. The key features and justifications are:

1. **Temporal envelope integration** ensures energy conservation regardless of time step size, critical for femtosecond pulses that may span only a few time steps.

2. **Gaussian spatial profile combined with Beer-Lambert absorption** captures the realistic energy distribution with strong surface localization due to the short optical penetration depth in metals.

3. **Separate film and gas treatment** prevents spurious heating of the gas phase and correctly localizes energy deposition to the metal film where it drives the LIFT process.

4. **Comprehensive energy accounting** with pulse-level verification provides confidence in simulation fidelity and enables rapid debugging during setup.

5. **Reflectivity and absorption-based attenuation** correctly models the fraction of incident energy coupled into the film, accounting for the ~50% reflection typical of metal surfaces.

6. **Geometric validation and diagnostics** catch common setup errors (misaligned focus, incorrect direction, inadequate mesh resolution) before they waste computational resources.

The model is designed to be robust, physically accurate, and comprehensively validated against both energy conservation principles and experimental LIFT literature. When properly configured with material properties from Wellershoff et al. (2000) for titanium and geometric parameters matching experimental LIFT setups, it provides the essential energy source term driving the subsequent electron heating, electron-phonon coupling, phase transitions, and material ejection that constitute the complete LIFT physical process.

---

## 14. Required Input Parameters for Case Setup

A complete femtosecond LIFT simulation requires careful specification of parameters across multiple configuration files. This section details all necessary parameters, their physical significance, typical values, and where they must be specified.

### 14.1 Laser Parameters (constant/laserProperties)

This dictionary contains all parameters defining the laser pulse characteristics and its interaction with the material.

#### Core Energy Parameters (REQUIRED)

**pulseEnergy** [J]
- Physical meaning: Total energy contained in a single laser pulse
- Typical range: 10 nJ to 10 μJ for LIFT applications
- Example value: 60 nJ (6e-8 J)
- Selection basis: Experimental ablation threshold for titanium is approximately 0.1-0.3 J/cm². For a 6 μm diameter spot (area ≈ 2.83×10^-7 cm²), 60 nJ gives fluence ≈ 0.21 J/cm², just above threshold
- Critical note: This is the TOTAL energy in the pulse, before accounting for reflection or absorption losses

**pulseWidth** [s]
- Physical meaning: Full-width at half-maximum (FWHM) duration of the Gaussian temporal pulse
- Typical range: 100 fs to 500 fs for femtosecond LIFT
- Example value: 200 fs (200e-15 s)
- Selection basis: Commercial Ti:Sapphire femtosecond lasers typically produce 100-200 fs pulses. This duration ensures non-equilibrium heating (pulse duration << electron-phonon coupling time ≈ 1-10 ps for titanium)
- Validation constraint: Model enforces 10 fs ≤ pulseWidth ≤ 10 ps

**wavelength** [m]
- Physical meaning: Wavelength of the laser radiation
- Typical range: 343 nm (UV) to 1030 nm (IR) depending on laser system
- Example value: 343 nm (343e-9 m) for UV Ti:Sapphire third harmonic
- Selection basis: Absorption coefficient is wavelength-dependent. Titanium shows strong absorption at UV wavelengths (343 nm is near peak absorption). For 800 nm fundamental, absorption is lower but still significant
- Physical importance: Determines penetration depth through wavelength-dependent absorption coefficient

#### Spatial Distribution Parameters (REQUIRED)

**spotSize** [m]
- Physical meaning: Beam diameter at the focal plane (defined at 1/e² intensity for Gaussian beams)
- Typical range: 5 μm to 50 μm for LIFT
- Example value: 6 μm (6e-6 m)
- Selection basis: Larger spots reduce peak intensity and may not reach ablation threshold; smaller spots concentrate energy but may be difficult to align. 5-10 μm is typical for single-pixel LIFT
- Mesh requirement: Beam diameter should span at least 8-10 cells for accurate Gaussian profile resolution

**focus** (x y z) [m]
- Physical meaning: 3D coordinates of the laser focal point in the computational domain
- Example value: (25e-6 20.0357e-6 5e-6) 
- Selection basis: 
  - x-coordinate: Centered in domain width (0 to 50 μm → center at 25 μm)
  - y-coordinate: CRITICAL - Must lie within the titanium film bounds. For film at y = 20.0 to 20.0714 μm, center is 20.0357 μm
  - z-coordinate: Centered in depth (0 to 10 μm → center at 5 μm)
- Common error: Placing focus outside film region results in zero film heating

**direction** (x y z) [dimensionless]
- Physical meaning: Unit vector specifying laser propagation direction
- Standard LIFT value: (0 -1 0)
- Physical interpretation: Laser propagates in negative y-direction (downward), entering the top surface of the film
- Critical note: Must be consistent with film orientation. For film perpendicular to y-axis with laser from above, direction MUST be (0 -1 0), not (0 1 0)
- Model behavior: Automatically normalized to unit vector

#### Optical Properties (REQUIRED)

**absorptionCoeff** [1/m]
- Physical meaning: Linear absorption coefficient α, determining exponential decay of intensity with depth according to Beer-Lambert law I(z) = I₀ × exp(-α × z)
- Typical range: 10^6 to 10^8 m^-1 for metals
- Example value: 1.03×10^8 m^-1 for titanium at 343 nm
- Literature basis: Palik's Handbook of Optical Constants provides refractive index data. For titanium at 343 nm: n ≈ 2.0-2.2, k ≈ 2.8. Absorption coefficient α = 4πk/λ = 4π(2.8)/(343×10^-9) ≈ 1.03×10^8 m^-1
- Physical significance: Penetration depth δ = 1/α ≈ 9.7 nm for this value, meaning 63% of transmitted energy is absorbed in the first 9.7 nm
- Titanium at 800 nm: α ≈ 5×10^6 m^-1 (δ ≈ 200 nm) from Wellershoff et al. (2000)

**reflectivity** [dimensionless, 0-1]
- Physical meaning: Fraction of incident intensity reflected at the air-metal interface (Fresnel reflection)
- Typical range: 0.3 to 0.6 for metals at optical wavelengths
- Example value: 0.35 for titanium at 343 nm
- Calculation: From complex refractive index n+ik, normal incidence reflectivity R = [(n-1)² + k²] / [(n+1)² + k²]. For Ti at 343 nm with n≈2.1, k≈2.8: R ≈ 0.35
- Titanium at 800 nm: R ≈ 0.50 (50% reflection)
- Model conversion: Internally converted to transmission T = 1 - R before application to beam intensity

**gasAbsorptionCoeff** [1/m] (OPTIONAL, default 0)
- Physical meaning: Absorption coefficient for the surrounding gas phase (air or argon)
- Typical value: 0 (gas transparent at LIFT wavelengths)
- When to use: Non-zero only if modeling plasma formation in air at very high intensities (>10^13 W/m²)
- Recommendation: Keep at 0 for standard LIFT to avoid spurious gas heating

#### Temporal Window (REQUIRED)

**laserStartTime** [s]
- Physical meaning: Simulation time when laser pulse begins
- Standard value: 0 (pulse starts at beginning of simulation)
- Alternative: Can be delayed if modeling pre-heating conditions or multiple pulses

**laserEndTime** [s]
- Physical meaning: Simulation time when laser source turns off
- Typical range: 200 ps to 2 ns
- Example value: 2×10^-10 s (200 ps)
- Selection basis: Must be long enough to capture entire pulse (pulse effectively ends at center + 3σ, where σ = FWHM/(2√(2ln2)) ≈ FWHM/2.35). For 200 fs pulse, 3σ ≈ 250 fs, but set window to 200 ps to allow observation of post-pulse dynamics
- Cost consideration: Longer windows increase computational cost; 200 ps is sufficient to observe initial ablation and jet formation

#### Profile Options (RECOMMENDED)

**gaussianProfile** [boolean]
- Physical meaning: Whether to use Gaussian (true) or uniform top-hat (false) spatial distribution
- Recommended value: true
- Physical justification: Real laser beams have Gaussian intensity profiles due to fundamental mode propagation (TEM₀₀)
- Top-hat alternative: Only for idealized studies or comparison with simplified models

**continuousLaser** [boolean]
- Physical meaning: Continuous wave (true) versus pulsed (false) operation
- Standard LIFT value: false
- When true: Ignores pulseWidth and applies constant intensity during active window

**temporalShape** [word] (OPTIONAL)
- Options: "gaussian" (default) or "square"
- Typical value: "gaussian" for realistic femtosecond pulses
- Square pulse: Instantaneous on/off, useful for limiting cases but not physical

**spatialShape** [word] (OPTIONAL)
- Options: "gaussian" (default) or "tophat"
- Should match gaussianProfile boolean

#### Advanced Parameters (OPTIONAL)

**pulseCenterTime** [s] (OPTIONAL)
- Physical meaning: Center time of Gaussian temporal profile
- Default behavior: If omitted, pulse is centered early in the [laserStartTime, laserEndTime] window
- When to specify: For delaying pulse peak or aligning with specific events
- Constraint: Automatically clamped to [laserStartTime, laserEndTime]

**maxVolumetricSource** [W/m³] (OPTIONAL)
- Physical meaning: Upper limit on volumetric heat source to prevent numerical instability
- Typical range: 10^24 to 10^25 W/m³
- Example value: 7×10^24 W/m³
- When needed: Very thin films with high absorption can produce extreme local source terms (>10^25 W/m³) that cause solver divergence
- Trade-off: Capping reduces accuracy near surface but improves stability. Model adjusts scaling to maintain energy conservation

**pulseFrequency** [1/s] (OPTIONAL, default 0)
- Physical meaning: Repetition rate for multiple pulses
- Example: 100 kHz = 10^5 s^-1 gives 10 μs between pulses
- Use case: Modeling burst-mode LIFT or high-repetition-rate systems
- Single pulse: Set to 0 or omit

**pulseDutyCycle** [dimensionless, 0-1] (OPTIONAL, default 1)
- Physical meaning: Fraction of pulse period during which laser is active
- Typical value: 1.0 (full duty cycle) for single pulses
- Use case: High-repetition-rate systems where pulses don't overlap

**scanVelocity** (x y z) [m/s] (OPTIONAL)
- Physical meaning: Velocity vector for translating laser focus during simulation
- Example: (100 0 0) moves focus at 100 m/s in x-direction
- Use case: Modeling line scanning or multi-pixel LIFT patterns
- Default: (0 0 0) for stationary focus

**incidenceAngle** [radians] (OPTIONAL, default 0)
- Physical meaning: Angle between beam propagation and surface normal
- Typical value: 0 (normal incidence)
- Effect: Modifies Fresnel reflectivity according to angle-dependent formulas
- Small angle approximation: For angles < 10°, effect on reflectivity is minor

**transmission** [dimensionless, 0-1] (OPTIONAL)
- Physical meaning: Direct specification of transmission coefficient, overriding reflectivity
- When to use: If transmission is known from experiment but reflectivity is not
- Relationship: T = 1 - R at normal incidence
- Priority: If specified (≥0), overrides reflectivity parameter

#### Film Geometry (ALTERNATIVE SPECIFICATION)

Option 1: Specify film bounds directly
- **filmYMin** [m]: Bottom surface y-coordinate
- **filmYMax** [m]: Top surface y-coordinate  
- Example: filmYMin = 20.0e-6, filmYMax = 20.0714e-6 for 71.4 nm film

Option 2: Specify thickness and center
- **filmThickness** (or filmThicknessExpected) [m]: Film thickness
- **filmCenterY** [m]: Y-coordinate of film center
- Example: filmThickness = 71.4e-9, filmCenterY = 20.0357e-6
- Model computes: filmYMin = center - thickness/2, filmYMax = center + thickness/2

**Critical requirement**: Film bounds must encompass the focus position for effective heating

#### Energy Conservation Tolerances (OPTIONAL)

**pulseEnergyToleranceRel** [dimensionless] (default 0.01)
- Physical meaning: Maximum acceptable relative error in pulse energy conservation
- Example: 0.05 = 5% tolerance
- Stricter than 1% may require very fine mesh

**pulseEnergyToleranceAbs** [J] (default 1e-12)
- Physical meaning: Absolute energy error tolerance
- Example: 5×10^-10 J = 500 pJ
- Relevant for low-energy pulses where relative error is less meaningful

### 14.2 Complete Example laserProperties File

```
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      laserProperties;
}

// Laser timing - 200 ps observation window
laserStartTime           0;
laserEndTime             2e-10;

// Film geometry (71.4 nm Ti film at y = 20.0 to 20.0714 μm)
filmThicknessExpected    71.4e-9;
filmCenterY              20.0357e-6;

// Energy parameters - 60 nJ pulse (~0.21 J/cm² fluence)
pulseEnergy              6e-8;
pulseWidth               200e-15;
wavelength               343e-9;

// Spatial profile - 6 μm Gaussian spot
spotSize                 6.0e-6;
focus                    (25e-6 20.0357e-6 5e-6);
direction                (0 -1 0);

// Optical properties - Ti at 343 nm (Palik data)
absorptionCoeff          1.03e8;
reflectivity             0.35;

// Profile options
gaussianProfile          true;
continuousLaser          false;
temporalShape            gaussian;
spatialShape             gaussian;

// Gas absorption (disabled for standard LIFT)
gasAbsorptionCoeff       0;

// Stability limit
maxVolumetricSource      7e24;

// Single pulse (no repetition)
pulseFrequency           0;
pulseDutyCycle           0;

// Energy verification tolerances
pulseEnergyToleranceRel  0.05;
pulseEnergyToleranceAbs  5e-10;
```

### 14.3 Two-Temperature Model Parameters (system/controlDict or constant/twoTemperatureProperties)

The femtosecond laser model couples to the two-temperature model, which requires its own parameter set:

**Ce0** [J/m³/K] - Electron heat capacity coefficient
- For titanium: Ce = γ×Te where γ ≈ 150 J/m³/K² (Wellershoff et al., 2000)
- Linear in electron temperature: Ce(Te) = γ×Te

**Cl** [J/m³/K] - Lattice heat capacity
- For titanium: Cl ≈ 2.4×10^6 J/m³/K (roughly constant)
- Much larger than Ce at room temperature: Cl >> Ce|T=300K

**G** [W/m³/K] - Electron-phonon coupling constant
- For titanium: G ≈ 2.4×10^16 to 3.8×10^17 W/m³/K (Lin et al., 2008; Wellershoff et al., 2000)
- Critical parameter controlling energy transfer rate from electrons to lattice
- Higher G → faster thermalization → closer to equilibrium behavior

**ke** [W/m/K] - Electron thermal conductivity
- For titanium: ke ≈ 21.9 W/m/K at room temperature
- Temperature-dependent: often ke(Te) = ke0 × (Tl/Te) accounting for electron-phonon scattering

### 14.4 Phase Change Parameters (system/controlDict: phaseChangeCoeffs)

**Tvapor** [K] - Vaporization temperature
- For titanium: 3287 K (normal boiling point at 1 atm)

**gasConstant** [J/kg/K] - Specific gas constant for Hertz-Knudsen equation  
- For titanium: R = R_universal / M_Ti = 8.314 / 0.04788 ≈ 173.6 J/kg/K

**evaporationCoeff** [dimensionless] - Evaporation accommodation coefficient
- Typical range: 0.01 to 0.1 for metals
- Example: 0.03 to 0.18 depending on calibration
- Physical meaning: Probability of evaporation attempt resulting in actual vapor departure

**relaxationRate** [1/s] - Implicit source term relaxation
- Typical value: 10^9 to 10^12 s^-1
- Higher values → tighter coupling but require smaller time steps

### 14.5 Recoil Pressure Parameters (system/controlDict: recoilPressureCoeffs)

**momentumAccommodationCoeff** [dimensionless]
- Typical value: 0.18
- Range: 0.1 to 1.0
- Physical meaning: Efficiency of momentum transfer from evaporating atoms to condensed phase

**alphaMin / alphaMax** [dimensionless]
- Define phase fraction window for recoil pressure application
- Typical: alphaMin = 0.001, alphaMax = 0.999
- Prevents evaluation in pure phases where interface is not present

### 14.6 Mesh Resolution Requirements

Based on the laser parameters, mesh must satisfy:

**Film thickness resolution**
- Minimum: 3-5 cells across film thickness
- For 71.4 nm film: cell height ≤ 15-20 nm in y-direction
- Recommended: 8-16 cells for accurate thermal gradients

**Penetration depth resolution**  
- Minimum: 2-3 cells within 1/α
- For α = 1.03×10^8 m^-1: δ = 9.7 nm requires cells ≤ 3-5 nm
- This is the most restrictive constraint for UV wavelengths

**Spot size resolution**
- Minimum: 8-10 cells across beam diameter
- For 6 μm spot: cell size ≤ 600-750 nm in x and z directions
- Resolves Gaussian profile to 1/e² intensity

**Typical cell dimensions for 71.4 nm Ti film, 343 nm wavelength, 6 μm spot:**
- Δx ≈ 500-625 nm (10-12 cells across 6 μm spot)
- Δy ≈ 4-5 nm near film (16 cells across 71.4 nm)
- Δz ≈ 500-625 nm (matches Δx for near-isotropic lateral resolution)

### 14.7 Time Step Requirements

**Pulse resolution criterion**
- Minimum: 10-20 time steps within pulse FWHM
- For 200 fs pulse: Δt ≤ 10-20 fs during pulse
- Ensures accurate temporal envelope integration

**CFL stability criterion**  
- Standard CFL limit: Co = |U|Δt/Δx < 0.5
- For post-ablation velocities ~100 m/s and Δx ~ 5 nm: Δt < 25 fs
- May need Δt ~ 1-5 fs during jet formation

**Thermal diffusion criterion**
- For explicit schemes: Δt < Δx²/(2α) where α = k/(ρCp)
- For titanium: α ~ 10^-5 m²/s, Δx ~ 5 nm → Δt < 1.25 ps
- Less restrictive than CFL for femtosecond LIFT

**Recommended strategy:**
- Initial pulse (0-500 fs): Δt = 0.5-1 fs (resolve pulse + electron heating)
- Thermalization (500 fs - 10 ps): Δt = 5-10 fs
- Ablation/jet (10 ps - 1 ns): Δt = 1-10 fs (CFL-limited)
- Late dynamics (>1 ns): Δt = 10-50 fs

### 14.8 Initial Conditions (0/ directory)

All temperature fields should be initialized consistently:

**T** [K] - Combined/mixture temperature
- Initial value: 300 K (room temperature)
- Boundary conditions: fixedValue or zeroGradient

**Te** [K] - Electron temperature
- Initial value: 300 K (thermal equilibrium with lattice)
- Must match Tl initially

**Tl** [K] - Lattice temperature  
- Initial value: 300 K
- Must match Te initially
- Critical: Te = Tl at t=0 to avoid spurious initial G(Te-Tl) source

**alpha1** [dimensionless] - Metal phase fraction
- Film region: 1.0 (pure metal)
- Gas region: 0.0 (pure gas)
- Set using setFieldsDict with boxToCell for film geometry

**p_rgh** [Pa] - Relative pressure
- Initial value: 0 Pa (atmospheric reference)
- Absolute pressure p = p_rgh + ρgh

**U** [m/s] - Velocity
- Initial value: (0 0 0) (quiescent fluid)

### 14.9 Common Parameter Selection Errors and Fixes

**Error 1: Zero Heating Despite Active Laser**
- Symptom: Model reports laser active but no temperature rise
- Causes:
  - Focus outside film bounds (check focus.y vs filmYMin/filmYMax)
  - Laser direction wrong (should be (0 -1 0) for standard LIFT)
  - absorptionCoeff = 0 or too small
- Fix: Verify geometric consistency and optical properties

**Error 2: Energy Conservation Violations**
- Symptom: Deposited energy differs from pulseEnergy by >5%
- Causes:
  - Mesh too coarse to integrate Gaussian×exponential profile
  - Time step too large to resolve 200 fs pulse
  - Spot extends outside domain (boundary energy loss)
- Fix: Refine mesh in beam path, reduce time step, or enlarge domain

**Error 3: Numerical Instability During Pulse**
- Symptom: Solver crashes with exploding temperatures or velocities
- Causes:
  - maxVolumetricSource not set, allowing Q > 10^26 W/m³
  - Time step too large for explicit thermal source
  - Electron heat capacity Ce too small at low Te
- Fix: Set maxVolumetricSource = 5e24 to 1e25, reduce Δt, check Ce formulation

**Error 4: Spurious Gas Heating**
- Symptom: Temperature rise in gas cells far from film
- Causes:
  - gasAbsorptionCoeff > 0 causing unphysical absorption
  - Film bounds incorrectly specified, treating gas as metal
- Fix: Set gasAbsorptionCoeff = 0, verify filmYMin/filmYMax encompass correct region

**Error 5: Laser Pulse Arrives at Wrong Time**
- Symptom: Temperature rise occurs too early or too late
- Causes:
  - pulseCenterTime specified incorrectly
  - laserStartTime ≠ 0 when simulation starts at t=0
- Fix: Omit pulseCenterTime for immediate pulse, or set to laserStartTime + pulseWidth

### 14.10 Parameter Validation Checklist

Before running simulation, verify:

**Geometric consistency:**
- [ ] filmYMin < focus.y < filmYMax (focus inside film)
- [ ] 0 < focus.x < domainLength_x (focus inside domain)
- [ ] 0 < focus.z < domainLength_z (focus inside domain)
- [ ] direction = (0, -1, 0) for standard LIFT geometry

**Energy scaling:**
- [ ] Fluence = pulseEnergy / (π(spotSize/2)²) ≈ 0.1-1 J/cm² (ablation threshold range)
- [ ] Peak intensity = pulseEnergy / (pulseWidth × π(spotSize/2)²) < 10^14 W/cm² (avoids plasma)

**Temporal consistency:**
- [ ] pulseWidth = 100-500 fs (femtosecond regime)
- [ ] laserEndTime > laserStartTime + 3×pulseWidth×2.35 (captures full pulse)
- [ ] Time step Δt ≤ pulseWidth/10 (resolves temporal profile)

**Optical properties:**
- [ ] absorptionCoeff matches wavelength (check Palik or Wellershoff)
- [ ] Penetration depth 1/α < filmThickness × 10 (significant absorption in film)
- [ ] 0 ≤ reflectivity ≤ 0.7 (physical range for metals)

**Mesh resolution:**
- [ ] Film has ≥5 cells in thickness
- [ ] Penetration depth has ≥3 cells
- [ ] Spot diameter has ≥8 cells

**Field initialization:**
- [ ] Te = Tl = 300 K everywhere at t=0 (equilibrium)
- [ ] alpha1 = 1 in film, 0 elsewhere (sharp interface)
- [ ] U = (0 0 0) everywhere (quiescent)

This comprehensive parameter specification ensures that the femtosecond laser model is properly configured to deliver physically accurate energy deposition profiles that drive realistic LIFT dynamics.

---

## References

1. Wellershoff, S. S., Hohlfeld, J., Güdde, J., & Matthias, E. (2000). The role of electron–phonon coupling in femtosecond laser damage of metals. Applied Physics A, 69(1), S99-S107.

2. Lin, Z., Zhigilei, L. V., & Celli, V. (2008). Electron-phonon coupling and electron heat capacity of metals under conditions of strong electron-phonon nonequilibrium. Physical Review B, 77(7), 075133.

3. Piqué, A., McGill, R. A., Chrisey, D. B., Leonhardt, D., Mslna, T. E., Spargo, B. J., ... & Mlsna, D. (2002). Growth of organic thin films by the matrix assisted pulsed laser evaporation (MAPLE) technique. Thin Solid Films, 355, 536-541.

4. Brown, M. S., Brasz, C. F., Ventikos, Y., & Arnold, C. B. (2010). Impulsively actuated jets from thin liquid films for high-resolution printing applications. Journal of Fluid Mechanics, 709, 341-370.

5. Feinaeugle, M., Alloncle, A. P., Delaporte, P., Sones, C. L., & Eason, R. W. (2012). Time-resolved shadowgraph imaging of femtosecond laser-induced forward transfer of solid materials. Applied Surface Science, 258(22), 8475-8483.

6. Chichkov, B. N., Momma, C., Nolte, S., Von Alvensleben, F., & Tünnermann, A. (1996). Femtosecond, picosecond and nanosecond laser ablation of solids. Applied Physics A, 63(2), 109-115.

7. Chen, J. K., Tzou, D. Y., & Beraun, J. E. (2006). A semiclassical two-temperature model for ultrafast laser heating. International Journal of Heat and Mass Transfer, 49(1-2), 307-316.

8. Rethfeld, B., Ivanov, D. S., Garcia, M. E., & Anisimov, S. I. (2017). Modelling ultrafast laser ablation. Journal of Physics D: Applied Physics, 50(19), 193001.

9. Palik, E. D. (Ed.). (1998). Handbook of optical constants of solids (Vols. 1-3). Academic press.
