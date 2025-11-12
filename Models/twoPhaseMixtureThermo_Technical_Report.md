# Technical Report: twoPhaseMixtureThermo Implementation for Femtosecond LIFT Simulation

## Executive Summary

The `twoPhaseMixtureThermo` class serves as the central thermophysical model for femtosecond Laser-Induced Forward Transfer simulations in OpenFOAM. This custom implementation extends the standard two-phase mixture framework to incorporate laser-matter interaction physics, phase-change phenomena, and two-temperature coupling required for ultrashort pulse laser ablation modeling. The class manages the coupled evolution of metal and gas phases during femtosecond laser irradiation, tracking energy deposition, electron-lattice equilibration, melting, vaporization, and recoil pressure generation.

## 1. Fundamental Architecture and Design Philosophy

### 1.1 Multiple Inheritance Strategy

The class employs multiple inheritance to combine three distinct physics domains:

**Base Class 1: psiThermo**
This provides the compressible thermodynamics framework where density depends on pressure and temperature. For the LIFT process, compressibility becomes critical when rapid vaporization generates high recoil pressures. The choice of psiThermo over rhoThermo reflects the physical reality that vapor expansion during ablation creates pressure fields that couple back to density changes.

**Base Class 2: twoPhaseMixture**
This manages the volume-of-fluid tracking for the metal-gas interface. The phase fraction field α₁ represents the volumetric fraction of metal, while α₂ = 1 - α₁ represents gas. The mixture framework handles the sharp interface between condensed metal and surrounding gas without requiring explicit interface reconstruction.

**Base Class 3: interfaceProperties**
This provides surface tension and interface curvature calculations essential for tracking the liquid-vapor interface during ablation. Surface tension forces become significant at the microscale dimensions typical of LIFT (tens of micrometers), affecting droplet formation and jet dynamics.

The multiple inheritance approach avoids code duplication while maintaining clean separation between thermodynamics, phase tracking, and interfacial physics. Each base class can be developed and validated independently.

### 1.2 Phase-Specific Thermo Packages

The implementation maintains separate thermodynamic packages for each phase through autoPtr smart pointers:

```cpp
autoPtr<rhoThermo> thermo1_;  // Metal phase thermodynamics
autoPtr<rhoThermo> thermo2_;  // Gas phase thermodynamics
```

This separation is physically motivated because metal and gas have fundamentally different thermodynamic behavior. The metal phase exhibits temperature-dependent thermal conductivity, specific heat, and density variations through solid-liquid-vapor transitions. The gas phase behaves as an ideal gas with its own equation of state. Maintaining separate thermo objects allows each phase to use appropriate constitutive relations without artificial mixing.

The autoPtr usage provides automatic memory management and enables runtime selection of different thermodynamic models for each phase through OpenFOAM's selection mechanism.

## 2. Phase Change Physics Implementation

### 2.1 Critical Temperature Thresholds

Three temperature thresholds govern phase transition physics:

**T_melt (Solidus Temperature):**
Represents the temperature at which the solid metal begins melting. For titanium, this is approximately 1941 K. Below this temperature, the material remains solid with its solid-phase thermal properties.

**T_vapor (Vaporization Temperature):**
Defines the onset of significant vaporization. For titanium at atmospheric pressure, this occurs around 3560 K. This represents the normal boiling point where liquid-to-vapor transition becomes rapid.

**maxPhaseChangeTemperature:**
An optional upper limit for phase-change calculations, defaulting to a very large value (GREAT in OpenFOAM terminology). This prevents unphysical behavior at extreme temperatures where the phase-change models may become invalid.

These thresholds must satisfy T_melt < T_vapor, which the constructor explicitly validates. This ordering ensures physically consistent phase progression: solid → liquid → vapor.

### 2.2 Latent Heat Treatment

The latent heat (L_v) quantifies the energy required for phase transformation:

```cpp
scalar latentHeat_;  // Units: J/kg
```

For titanium vaporization, published values range from 8.88 × 10⁶ to 9.83 × 10⁶ J/kg depending on temperature. The simulation uses this parameter in the phase-change source term calculation. When material evaporates, the corresponding latent heat is removed from the lattice temperature field, creating a cooling effect that moderates further vaporization—a self-limiting feedback mechanism observed experimentally.

The implementation tracks latent heat consumption through the cumulative phase-change energy diagnostics, enabling validation against energy conservation principles.

### 2.3 Hertz-Knudsen Evaporation Model

The physical evaporation rate follows the Hertz-Knudsen relation, which describes molecular kinetics at the liquid-vapor interface:

**j = α_evap × (P_sat(T) - P_vapor) / sqrt(2πm_vapor k_B T)**

Where:
- j: mass flux [kg/(m²·s)]
- α_evap: accommodation coefficient (evaporation/condensation probability)
- P_sat(T): saturation vapor pressure at temperature T
- P_vapor: actual vapor pressure
- m_vapor: mass of vapor molecule
- k_B: Boltzmann constant
- T: interface temperature

The implementation stores:
```cpp
scalar gasConstant_;         // R_gas = k_B / m_vapor [J/(kg·K)]
scalar evaporationCoeff_;    // α_evap (typically 0.18 for titanium)
```

The gas constant R_gas must be supplied in the configuration as it depends on the specific material being ablated. For titanium (atomic mass ~47.87 g/mol):
R_gas = R_universal / M = 8314 J/(kmol·K) / 47.87 kg/kmol ≈ 173.6 J/(kg·K)

The accommodation coefficient α_evap = 0.18 is the default value, consistent with molecular dynamics studies of metal evaporation. This parameter represents the probability that an incident vapor molecule condenses upon striking the liquid surface, or conversely, that a surface molecule evaporates when thermally excited.

### 2.4 Phase-Fraction Windows for Numerical Stability

The phase-change calculations apply only within specific α ranges:

```cpp
scalar alphaMin_;  // Default: 0.01
scalar alphaMax_;  // Default: 1.0
```

This windowing serves two purposes:

**Physical Justification:**
At very low metal fractions (α₁ < 0.01), the remaining material exists as dispersed droplets or vapor. The continuum assumption breaks down, and applying bulk evaporation models becomes questionable. The alphaMin threshold prevents calculations in this regime.

**Numerical Justification:**
Near α₁ = 0 or α₁ = 1, small numerical errors in phase fraction can cause large relative changes in calculated properties. The window [0.01, 1.0] ensures phase-change source terms are computed only where the interface is well-resolved and the continuum approximation remains valid.

### 2.5 Temporal Relaxation and Stability

Three time-scale parameters control temporal stability:

```cpp
scalar dtFloor_;           // Minimum time-step: 1e-12 s
scalar relaxationTime_;    // Source term relaxation: 1e-12 s
dimensionedScalar maxPhaseChangeSource_;  // Maximum power density
```

**Time-Step Floor (dtFloor_):**
Prevents excessively small time steps during rapid phase change. At 1 picosecond (1e-12 s), this matches the electron-phonon coupling timescale, ensuring phase-change rates don't evolve faster than the underlying thermal dynamics.

**Relaxation Time (τ_relax):**
The phase-change source term uses implicit-explicit splitting:
- Implicit part: proportional to 1/τ_relax, stabilizes the calculation
- Explicit part: represents the driving force for phase change

This splitting prevents oscillations that can occur when evaporation rates change rapidly. The relaxation time essentially says "the system adjusts to phase-change demands over ~1 picosecond," which is physically reasonable for femtosecond-laser-induced processes.

**Maximum Source Term (maxPhaseChangeSource_):**
Caps the volumetric power density [W/m³] extracted by phase change. Without this limiter, cells at very high temperature could predict instantaneous vaporization, causing numerical blow-up. The default GREAT (10³⁰ in OpenFOAM) effectively disables this limit, but users can set finite values for problematic cases.

### 2.6 Activation Windows for Time-Dependent Control

```cpp
List<Tuple2<scalar, scalar>> activationWindows_;
```

This feature allows enabling phase change only during specific time intervals, useful for:
1. Initialization phases where thermal fields equilibrate
2. Diagnostic runs where you want to isolate laser heating from evaporation
3. Staged simulations that progressively add physics complexity

Each tuple contains [t_start, t_end]. Phase change only proceeds when current time falls within one of these windows. Empty list means "always active."

## 3. Two-Temperature Model Coupling

### 3.1 Lattice Heat Capacity Storage

```cpp
dimensionedScalar ClTTM_;  // Default: 2.5e6 J/(m³·K)
```

The two-temperature model solves separate equations for electron temperature (T_e) and lattice temperature (T_l). The twoPhaseMixtureThermo class caches the lattice heat capacity because phase-change cooling appears as a source term in the lattice equation but not the electron equation.

**Physical Basis:**
Femtosecond laser pulses deposit energy directly into the electron system. The electrons then transfer this energy to the lattice via electron-phonon collisions on a timescale of 1-10 picoseconds. Phase change (melting, vaporization) occurs when lattice atoms gain sufficient energy to overcome binding forces—thus it couples to T_l, not T_e.

**Default Value Justification:**
C_l = 2.5 × 10⁶ J/(m³·K) is a typical value for solid metals. For titanium specifically:
- Density ρ ≈ 4500 kg/m³
- Specific heat c_p ≈ 520 J/(kg·K)
- Volumetric heat capacity = ρ × c_p ≈ 2.34 × 10⁶ J/(m³·K)

The default falls within the correct order of magnitude. The exact value should come from the two-temperature model configuration (twoTemperatureProperties dictionary) and gets updated via `setClTTM()` during initialization.

### 3.2 Phase-Change Source Term Structure

The simulation tracks three distinct fields:

```cpp
volScalarField phaseChangeSource_;         // Units: K/s
volScalarField phaseChangeRelaxCoeff_;    // Units: 1/s  
volScalarField dgdt_;                      // Units: 1/s
```

**phaseChangeSource_ [K/s]:**
This represents the rate of temperature change due to phase transition. It enters the lattice temperature equation as:

∂T_l/∂t = ... + phaseChangeSource_

The negative values indicate cooling (evaporation removing heat), positive values indicate heating (condensation releasing latent heat). The K/s units mean "temperature change per second due to phase change."

**phaseChangeRelaxCoeff_ [1/s]:**
This is the implicit diagonal coefficient for the temperature equation. When solving:

∂T_l/∂t = S_explicit + S_implicit × T_l

The implicit part S_implicit = -phaseChangeRelaxCoeff_ helps stability. The 1/s dimension means it represents a rate coefficient.

**dgdt_ [1/s]:**
The actual volumetric rate of phase-fraction change:

∂α/∂t = dgdt_

This field drives the evolution of the phase fraction field and couples to the momentum equation through the interface-capturing algorithm.

### 3.3 Mass Flux Field

```cpp
volScalarField phaseChangeMassFlux_;  // Units: kg/(m²·s)
```

This stores the physical mass flux across the liquid-vapor interface, computed from the Hertz-Knudsen relation. While dgdt_ gives the volumetric phase-change rate, phaseChangeMassFlux_ provides the surface-specific rate needed for:
1. Recoil pressure calculation (momentum transfer from escaping vapor)
2. Material loss tracking (total mass ablated)
3. Diagnostic output and validation against experiments

The units kg/(m²·s) represent mass crossing a unit interface area per unit time.

## 4. Laser Energy Deposition

### 4.1 Volumetric Source Field

```cpp
volScalarField Q_laser_;  // Units: W/m³
```

This field stores the spatial distribution of laser power density. The laser module (separate from twoPhaseMixtureThermo) calculates Q_laser_ based on:
1. Laser pulse characteristics (energy, duration, wavelength)
2. Beam profile (Gaussian, top-hat, etc.)
3. Absorption depth (Beer-Lambert law)
4. Temporal pulse shape

The twoPhaseMixtureThermo class provides access to this field but doesn't compute it internally. The separation allows different laser models to be used without modifying the thermophysics class.

### 4.2 Energy Conservation Tracking

The class maintains no internal cumulative energy counters itself, but cooperates with the two-temperature model's energy audit by providing:
- Access to phase-change energy extraction via phaseChangeSource_
- Interface for updating Q_laser_ from external laser model
- Consistent coupling to lattice temperature equation

The two-temperature model then tracks:
- Cumulative laser input
- Cumulative phase-change energy loss
- Total metal internal energy
- Energy conservation error

This distributed responsibility keeps each module focused on its primary physics while enabling system-wide energy accounting.

## 5. Transport Property Management

### 5.1 Phase-Specific Properties

```cpp
dimensionedScalar nu1_;   // Metal kinematic viscosity [m²/s]
dimensionedScalar nu2_;   // Gas kinematic viscosity [m²/s]
dimensionedScalar rho1_;  // Metal density [kg/m³]
dimensionedScalar rho2_;  // Gas density [kg/m³]
```

These properties are read from the transportProperties dictionary during construction. They serve two purposes:

**Direct Use:**
Provide viscosity and density for momentum equation calculations. The mixture properties are computed via volume-weighted or mass-weighted averaging.

**Validation Reference:**
The thermo packages (thermo1_, thermo2_) compute density from their equations of state. The rho1_ and rho2_ values provide reference densities for consistency checking and initialization.

### 5.2 Mixture Property Computation Methods

The class implements several mixture property calculation strategies:

**Volume-Weighted Average (simple mixing):**
Used for properties that vary smoothly across the interface:
φ_mix = α₁ × φ₁ + α₂ × φ₂

Applied to: thermal conductivity (kappa), thermal diffusivity (alpha)

**Mass-Weighted Average:**
Used for specific properties (per unit mass):
φ_mix = (α₁ρ₁φ₁ + α₂ρ₂φ₂) / (α₁ρ₁ + α₂ρ₂)

Applied to: specific heat capacity (Cp, Cv), enthalpy (h)

The mass-weighting is physically correct because specific properties must account for how much mass contributes each property value. Volume-weighting would give incorrect results when densities differ significantly (as they do for metal vs. gas, ρ_metal/ρ_gas ~ 10,000).

### 5.3 Temperature Inversion Algorithm

The class implements Newton iteration for finding temperature from enthalpy (THE method):

Given: h, p, T₀ (initial guess)
Find: T such that h(T,p) = h

**Algorithm:**
1. Start with T = T₀
2. For iteration i = 1 to maxIter (50):
   - Compute h(T_i, p) using mixture enthalpy formula
   - Compute Cp(T_i, p) for derivative
   - Update: T_{i+1} = T_i - (h(T_i, p) - h_target) / Cp(T_i, p)
   - Check convergence: if |T_{i+1} - T_i| < tolerance (1e-6 K), return T_{i+1}
3. If not converged, emit warning and return best guess

This inversion is required because OpenFOAM's energy equation solves for enthalpy, but many physical processes (phase change, thermal conduction) depend on temperature. The Newton method converges rapidly for well-behaved equations of state.

**Mixture Enthalpy Evaluation:**
At each iteration, the mixture enthalpy is computed by evaluating each phase's enthalpy at the current temperature guess, then mass-averaging:

h_mix = (α₁ρ₁h₁(T,p) + α₂ρ₂h₂(T,p)) / (α₁ρ₁ + α₂ρ₂)

The derivative (specific heat) follows similarly:

Cp_mix = (α₁ρ₁Cp₁(T,p) + α₂ρ₂Cp₂(T,p)) / (α₁ρ₁ + α₂ρ₂)

## 6. Computational Methods and Numerical Implementation

### 6.1 Phase Change Source Computation

The `computePhaseChange()` method implements the core evaporation physics:

**Step 1: Identify Active Cells**
Only cells meeting all criteria participate:
- Metal fraction: alphaMin_ ≤ α₁ ≤ alphaMax_
- Temperature threshold: T_l > T_vapor (if onlyAboveVapor_ enabled)
- Time window: current time within activationWindows_
- Metal presence: α₁ > metalFractionFloor

**Step 2: Compute Saturation Pressure**
For each active cell, evaluate P_sat(T_l). Common correlations include:
- Clausius-Clapeyron: P_sat = P₀ exp(-L_v/R(1/T - 1/T₀))
- Antoine equation: log₁₀(P_sat) = A - B/(T + C)

**Step 3: Evaluate Hertz-Knudsen Mass Flux**
j = α_evap × (P_sat - P_ambient) / sqrt(2πm_vapor k_B T_l)

The pressure difference (P_sat - P_ambient) drives evaporation when positive (superheated liquid) and condensation when negative (subcooled vapor).

**Step 4: Convert to Volumetric Rate**
The rate of phase fraction change:
∂α₁/∂t = -j × A_interface / (ρ_metal × V_cell)

Where A_interface is the interface area within the cell, estimated from α gradients.

**Step 5: Temperature Source Term**
The cooling rate due to evaporation:
phaseChangeSource = -(∂α₁/∂t × ρ_metal × L_v) / (α₁ × ρ_metal × C_l)
                  = -(∂α₁/∂t × L_v) / (α₁ × C_l)

This gives temperature rate [K/s] from the energy loss rate.

**Step 6: Implicit Coefficient**
To improve stability, part of the source is treated implicitly:
phaseChangeRelaxCoeff = 1 / τ_relax

This modifies the temperature equation to:
∂T_l/∂t = ... + phaseChangeSource - phaseChangeRelaxCoeff × (T_l - T_reference)

The implicit term resists runaway temperature changes.

**Step 7: Limiting**
Apply maximum source term limiter if maxPhaseChangeSource_ is finite:
phaseChangeSource = max(phaseChangeSource, -maxPhaseChangeSource_)

The negative sign is because evaporation cools (negative source).

### 6.2 Thermo Correction Sequence

The `correctThermo()` method updates all thermodynamic properties in a specific order:

**Phase 1: Individual Phase Updates**
Call thermo1_->correct() and thermo2_->correct() to update each phase's properties (ρ, Cp, Cv, μ, k) based on current T and p fields.

**Phase 2: Mixture Property Updates**
Update psi (compressibility), mu (viscosity), and alpha (thermal diffusivity) using the mixture rules.

**Phase 3: Boundary Condition Updates**
Call correctBoundaryConditions() on T, ensuring boundary values are consistent with interior.

This sequence is critical: individual phases must update before mixing, and boundaries update last to maintain consistency at inlet/outlet/walls.

### 6.3 Overall Correction Method

The `correct()` method serves as the top-level update:

**Step 1: Interface Properties**
Update surface tension σ(T) and interface curvature κ from current α field.

**Step 2: Mixture Thermodynamics**
Call correctThermo() to update all thermal properties.

**Step 3: Phase-Change Coupling**
Call computePhaseChange() to calculate evaporation source terms based on updated temperatures.

**Step 4: Density Update**
Update mixture density: ρ = α₁ρ₁ + α₂ρ₂ using the new phase densities.

This order ensures causality: interface geometry → properties → phase change → density. The density update happens last because phase change can modify α₁.

## 7. Required Input Parameters

### 7.1 Essential Entries in transportProperties

```
phases (metal air);

metal
{
    transportModel  Newtonian;
    nu              nu [0 2 -1 0 0 0 0] 3.5e-7;     // Kinematic viscosity
    rho             rho [1 -3 0 0 0 0 0] 4506;       // Density
}

air
{
    transportModel  Newtonian;
    nu              nu [0 2 -1 0 0 0 0] 1.48e-5;
    rho             rho [1 -3 0 0 0 0 0] 1.225;
}

sigma           sigma [1 0 -2 0 0 0 0] 0.07;         // Surface tension
```

**Justification for Titanium Values:**
- Metal kinematic viscosity: Liquid Ti at melting point exhibits ν ≈ 3.5 × 10⁻⁷ m²/s
- Metal density: Solid Ti has ρ = 4506 kg/m³; liquid decreases ~10% to 4110 kg/m³
- Air viscosity: Standard value at 300 K
- Surface tension: Liquid Ti-vapor interface σ ≈ 1.65 N/m (source uses 0.07 as placeholder)

### 7.2 Essential Entries in controlDict.phaseChangeCoeffs

```
phaseChangeCoeffs
{
    latentHeat      latentHeat [0 2 -2 0 0 0 0] 9.83e6;     // J/kg
    T_melt          T_melt [0 0 0 1 0 0 0] 1941;            // K
    T_vapor         T_vapor [0 0 0 1 0 0 0] 3560;           // K
    
    gasConstant     gasConstant [0 2 -2 -1 0 0 0] 173.6;   // J/(kg·K)
    evaporationCoeff evaporationCoeff [0 0 0 0 0 0 0] 0.18;
    
    alphaMin        alphaMin [0 0 0 0 0 0 0] 0.01;
    alphaMax        alphaMax [0 0 0 0 0 0 0] 1.0;
    
    maxPhaseChangeTemperature maxPhaseChangeTemperature [0 0 0 1 0 0 0] 8000;
}
```

**Critical Notes:**
1. gasConstant = R_universal / M_molecular must match your material
2. latentHeat values in literature vary by 10%; use experimentally validated values
3. Accommodation coefficient 0.18 comes from MD simulations (Zhigilei et al.)
4. Temperature bounds should bracket your expected temperature range

### 7.3 Essential Entries in twoTemperatureProperties

```
twoTemperatureProperties
{
    Ce              Ce [1 -1 -2 -1 0 0 0] 210;             // Electron heat capacity
    Cl              Cl [1 -1 -2 -1 0 0 0] 2.34e6;         // Lattice heat capacity
    G               G [1 -1 -3 -1 0 0 0] 3.8e17;          // Coupling constant
    De              De [0 2 -1 0 0 0 0] 1e-4;             // Electron diffusivity
    
    minTe           minTe [0 0 0 1 0 0 0] 300;
    maxTe           maxTe [0 0 0 1 0 0 0] 15000;
    minTl           minTl [0 0 0 1 0 0 0] 300;
    maxTl           maxTl [0 0 0 1 0 0 0] 8000;
}
```

**Physical Basis for Titanium:**
- Ce = 210 J/(m³·K) at room temperature; increases with T_e (should be Function1)
- Cl = ρ × c_p = 4506 kg/m³ × 520 J/(kg·K) = 2.34 × 10⁶ J/(m³·K)
- G values in literature range 2.5-5.5 × 10¹⁷ W/(m³·K); 3.8×10¹⁷ is mid-range
- De determines electron thermal conductivity: k_e = De × Ce

The electron-phonon coupling constant G is THE most critical parameter. Values from:
- Wellershoff et al. (1999): G = 3.6 × 10¹⁷ W/(m³·K) for Ti (Reference: Appl. Phys. A 69, S99)
- Lin et al. (2008): Measured G = 3.8 × 10¹⁷ W/(m³·K) (Reference: Phys. Rev. B 77, 075133)

### 7.4 Initial Conditions (0/ directory)

**Te and Tl:**
Initialize both to ambient (300 K) everywhere:
```
dimensions      [0 0 0 1 0 0 0];
internalField   uniform 300;
```

**alpha.metal:**
Must represent initial geometry (thin film + substrate):
```
dimensions      [0 0 0 0 0 0 0];
internalField   uniform 0;  // Most domain is air

boundaryField
{
    film            // Thin film region
    {
        type    fixedValue;
        value   uniform 1;
    }
}
```

Then use setFields to initialize the exact film geometry.

**U (velocity):**
Start from rest:
```
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0 0 0);
```

**p_rgh (reduced pressure):**
Atmospheric:
```
dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform 101325;
```

## 8. Output Fields and Diagnostics

### 8.1 AUTO_WRITE Fields

The class automatically writes these fields at each output time:

**Q_laser [W/m³]:**
Laser power density distribution. Peak values during pulse can reach 10²⁰ W/m³ for femtosecond pulses. Verify that integrated volume matches input pulse energy.

**phaseChangeSource [K/s]:**
Temperature change rate due to evaporation. Negative values (cooling) should concentrate near liquid-vapor interface. Typical magnitudes: 10¹² - 10¹⁴ K/s during active evaporation.

**phaseChangeRelaxCoeff [1/s]:**
Implicit stabilization coefficient. Should be uniform at 1/τ_relax ≈ 10¹² s⁻¹ in active phase-change regions, zero elsewhere.

**dgdt [1/s]:**
Phase fraction evolution rate. Positive means α₁ increasing (condensation), negative means decreasing (evaporation). Check that integral matches material loss rate.

**phaseChangeMassFlux [kg/(m²·s)]:**
Surface-specific mass flux. Typical values during femtosecond LIFT: 10³ - 10⁵ kg/(m²·s). Compare to experimental ablation rates to validate model.

### 8.2 Derived Diagnostics

**Temperature Fields (Te, Tl):**
Monitor peak electron and lattice temperatures. For femtosecond pulses on Ti:
- Peak Te can reach 10,000 - 50,000 K (electrons far from equilibrium)
- Peak Tl remains below 5,000 - 8,000 K (limited by electron-phonon coupling rate)
- The Te - Tl difference indicates non-equilibrium strength

**Velocity Magnitude:**
Calculate |U| to identify jet velocity. Experimental LIFT jets typically reach 30-100 m/s. Your simulation should reproduce this order of magnitude.

**Material Loss:**
Integrate dgdt over metal volume:
dm/dt = ∫∫∫ ρ₁ × dgdt dV

Compare cumulative mass loss to experimental ablation depth measurements.

## 9. Implementation Logic Summary

The twoPhaseMixtureThermo class orchestrates femtosecond LIFT physics through these key mechanisms:

**Energy Flow Path:**
Laser → Electrons (Te) → Lattice (Tl) → Phase Change → Mass Loss
              ↓                           ↓
         Diffusion              Gas Coupling

**Coupling Strategy:**
1. Laser model deposits Q_laser into electron equation
2. Two-temperature model evolves Te and Tl with electron-phonon coupling
3. twoPhaseMixtureThermo monitors Tl and triggers evaporation when Tl > T_vapor
4. Phase-change cooling feeds back into Tl equation as negative source
5. Evaporation mass loss updates α field, affecting all subsequent properties
6. Interface-capturing calculates recoil pressure from mass flux
7. Momentum equation responds to recoil pressure, accelerating material

**Critical Feedback Loops:**
- Evaporation cools → reduces further evaporation (negative feedback, stabilizing)
- Recoil pressure pushes material → thins film → increases heating rate (positive feedback during pulse)
- Electron cooling via diffusion → reduces lattice heating rate → slows melting

**Numerical Stability Measures:**
- Phase-fraction windowing prevents calculations in dispersed/pure phases
- Relaxation time adds implicit damping to phase-change sources
- Maximum source limiter prevents instantaneous vaporization
- Time-step floor prevents collapse during rapid transients
- Property bounds (minTe, maxTe) prevent unphysical temperatures

## 10. Validation Strategy

To validate the twoPhaseMixtureThermo implementation:

**Test 1: Energy Conservation**
Monitor cumulative energy balance:
E_laser_input = E_metal_increase + E_phase_change + E_gas_coupling

Should match within 1-2% (accounting for numerical diffusion).

**Test 2: Phase-Change Onset**
Apply steady heating to liquid metal. Evaporation should initiate when Tl crosses T_vapor with mass flux consistent with Hertz-Knudsen prediction.

**Test 3: Jet Velocity**
Compare final jet velocity to experimental LIFT measurements (30-100 m/s range). Significant deviation indicates errors in recoil pressure or momentum coupling.

**Test 4: Ablation Depth**
Measure final film thickness reduction. Convert to equivalent evaporation depth and compare to experiments at same fluence.

**Test 5: Temperature Profiles**
Compare Te and Tl evolution to two-temperature model analytical solutions for simple geometries (semi-infinite solid with step heat input).

## 11. Common Issues and Troubleshooting

**Issue 1: No evaporation despite high temperature**
Check:
- Is Tl actually exceeding T_vapor in metal cells?
- Is α₁ within [alphaMin, alphaMax]?
- Are activationWindows configured correctly?
- Is gasConstant specified (not defaulted to zero)?

**Issue 2: Excessive evaporation causes α₁ → 0 everywhere**
Reduce:
- evaporationCoeff (try 0.1 instead of 0.18)
- Increase relaxationTime (try 5e-12 instead of 1e-12)
- Set maxPhaseChangeSource to reasonable limit (1e18 W/m³)

**Issue 3: Temperature inversion fails (THE method)**
Usually means enthalpy-temperature relation is non-monotonic. Check:
- Are phase thermo models physically reasonable?
- Is specific heat always positive?
- Try different starting temperature T₀

**Issue 4: Interface becomes diffuse**
Not a twoPhaseMixtureThermo issue, but check:
- Increase interface compression coefficient (cAlpha)
- Reduce time step
- Use finer mesh near interface

**Issue 5: Energy not conserved**
- Verify phase-change latent heat matches Tl equation source term dimensions
- Check that ClTTM is set correctly from two-temperature properties
- Enable energyAudit diagnostics in twoTemperatureProperties

## 12. Conclusion

The twoPhaseMixtureThermo class implements a comprehensive thermophysical framework for femtosecond LIFT simulations. Its design philosophy prioritizes physical accuracy through separate phase treatment, numerical stability through windowing and relaxation, and diagnostic capability through extensive field output. Proper configuration requires careful attention to material properties, particularly electron-phonon coupling constants and accommodation coefficients, which control the energy transfer rates that govern ablation dynamics. The implementation successfully couples laser energy deposition, two-temperature thermal evolution, Hertz-Knudsen evaporation kinetics, and recoil pressure generation into a unified framework suitable for predicting experimentally observed LIFT jet velocities and material ejection rates.
