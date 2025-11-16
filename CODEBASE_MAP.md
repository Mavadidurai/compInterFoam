# compInterFoam Codebase Structure - Comprehensive Map

## Overview
**compInterFoam** is an advanced OpenFOAM-based solver for Laser-Induced Forward Transfer (LIFT) simulations. It extends the standard `compressibleInterFoam` solver with femtosecond laser physics, two-temperature electron-lattice coupling, and advanced evaporation/recoil modeling for metal-gas phase systems.

**Root Directory:** `/home/user/compInterFoam`

---

## 1. MAIN SOLVER FILE (ENTRY POINT)

### `/home/user/compInterFoam/compInterFoam.C` (1,243 lines)
**Primary solver executable and main loop**

- **Purpose:** Main application entry point; orchestrates time-stepping loop and couples all physics models
- **Key Functions:**
  - Initializes OpenFOAM framework (fvCFD.H, dynamicFvMesh.H)
  - Sets up PIMPLE algorithm for pressure-velocity coupling
  - Manages adaptive time-stepping with Courant number limits
  - Includes LIFT process tracking and diagnostics output
  - Implements `liftProcessTracker()` function for real-time state snapshots
  - Reads configuration from `system/controlDict`

- **Key Includes:**
  - `twoPhaseMixtureThermo.H` - Thermophysical properties & phase change
  - `femtosecondLaserModel.H` - Laser energy deposition
  - `twoTemperatureModel.H` - Electron-lattice temperature coupling
  - `advancedInterfaceCapturing.H` - Recoil pressure and interface treatment
  - `pimpleControl.H` - PIMPLE pressure-velocity algorithm control

- **Output:** Executable compiled to `$(FOAM_USER_APPBIN)/compInterFoam`

---

## 2. EQUATION FILES (COUPLED PHYSICS)

### `/home/user/compInterFoam/UEqn.H` (654 lines)
**Momentum equation assembly**

- **Purpose:** Assembles the momentum equation using the PIMPLE/SIMPLE algorithm
- **Physics Included:**
  - Unsteady term: `∂(ρU)/∂t`
  - Convection: `∇·(ρU⊗U)`
  - Reynolds stress divergence from transport model
  - Recoil pressure force coupling (from `advancedInterfaceCapturing`)
  - MRF (Multiple Reference Frame) support for rotating regions
- **Solver Configuration:** 
  - Reads from `system/fvSolution` → `solvers["U.*"]`
  - Uses smoothSolver with symGaussSeidel smoother
- **Recoil Force Integration:** Checks for `recoilForceReady` and adds recoil pressure gradient to `HbyA`

### `/home/user/compInterFoam/alphaEqn.H` (374 lines)
**Phase-fraction (volume of fluid) equation**

- **Purpose:** Advects the metal phase fraction (alpha.metal) across domains
- **Physics Included:**
  - VOF interface tracking using MULES (Multidimensional Universal Limiter for Explicit Solution)
  - Interface compression with `interfaceCompression` scheme
  - Sub-cycling support for enhanced stability on small time-steps
  - Phase-change mass source coupling (Sp/Su terms from `twoPhaseMixtureThermo`)
  - Thin-film correction (optional)
- **Solver Configuration:**
  - Reads from `system/fvSolution` → `solvers["alpha.*"]`
  - Supports `nAlphaSubCycles` for sub-cycling stability
  - Compression limiter to prevent overshooting
- **Mass Conservation:**
  - Updates `rhoPhi` (mass flux) based on alpha transport
  - Bounds alpha between 0 and 1 after each iteration

### `/home/user/compInterFoam/TEqn.H` (1,431 lines)
**Temperature equation (complex multi-temperature coupling)**

- **Purpose:** Evolves the gas-phase temperature field, integrating laser heating and two-temperature coupling
- **Physics Included:**
  - Unsteady thermal diffusion with anisotropic conductivity
  - Laser heat source from `femtosecondLaserModel`
  - Gas-metal interfacial heat exchange via `twoTemperatureModel`
  - Phase-change energy sink from evaporation/condensation
  - Temperature bounding (optional clamping to physical limits)
  - Fallback to lattice temperature (Tl) in metal-dominated regions
- **Solver Configuration:**
  - Reads from `system/fvSolution` → `solvers["T"]` / `solvers["TFinal"]`
  - Uses smoothSolver with symGaussSeidel smoother
  - Thermal relaxation factor for stability
- **Coupling:**
  - Receives `gasMetalHeatFlux` from `twoTemperatureModel`
  - Feeds temperature field back to phase-change model

### `/home/user/compInterFoam/pEqn.H` (517 lines)
**Pressure equation (PIMPLE pressure-correction)**

- **Purpose:** Pressure-velocity coupling using PIMPLE/SIMPLE iterative scheme
- **Physics Included:**
  - Pressure-Poisson equation assembly via implicit momentum coupling
  - Recoil pressure effects on RHS through `advancedInterfaceCapturing`
  - Gravity/hydrostatic reference pressure (p_rgh, not absolute p)
  - Density-weighted velocity correction
  - Interface normal stress (surface tension) term
- **Solver Configuration:**
  - Reads from `system/fvSolution` → `solvers["p_rgh"]` / `solvers["p_rghFinal"]`
  - Uses GAMG (Geometric Algebraic Multi-Grid) with DIC smoothing
  - Reference cell: cell 0, reference value: 101325 Pa
- **PIMPLE Coupling:**
  - Iteration loop: `nOuterCorrectors` (typically 10) × `nCorrectors` (typically 3)
  - Pressure relaxation factor from `relaxationFactors`

---

## 3. PHYSICS MODEL FILES

### `/home/user/compInterFoam/twoPhaseMixtureThermo.H` (451 lines)
### `/home/user/compInterFoam/twoPhaseMixtureThermo.C` (2,037 lines)
**Thermophysical properties and phase-change model**

- **Purpose:** Manages two-phase mixture thermodynamics, phase-change (evaporation), and recoil pressure source
- **Key Data Members:**
  - `thermo1_`, `thermo2_` - Individual phase thermodynamic packages (rhoThermo)
  - `latentHeat_` - Latent heat of vaporization (J/kg)
  - `T_melt_`, `T_vapor_` - Phase-change temperature thresholds (K)
  - `Q_laser_` - Volumetric laser heating source (W/m³)
  - `phaseChangeSource_` - Explicit phase-change mass source (1/s)
  - `phaseChangeRelaxCoeff_` - Implicit phase-change relaxation coefficient
  - `gasConstant_` - Hertz-Knudsen relation constant (J/(kg·K))
  - `evaporationCoeff_` - Evaporation accommodation coefficient (0–1)
- **Key Methods:**
  - `correct()` - Updates mixture properties based on current fields
  - `latentHeat()` - Returns latent heat
  - `T_melt()`, `T_vapor()` - Return thresholds
  - `phaseChangeAlphaMin()` - Alpha range for phase-change activation
  - `Q_laser()`, `gasMetalExchangeCoeffField()` - Access source fields
  - `setClTTM(Cl)` - Cache lattice heat capacity from two-temperature model
- **Configuration:**
  - Primary: `system/controlDict` → `phaseChangeCoeffs` block
  - Secondary: `constant/thermophysicalProperties*` files (per-phase thermo data)
  - Activation windows: Time ranges when phase-change is active

**Phase-Change Model (Hertz-Knudsen):**
- Evaporation rate calculated from Clausius-Clapeyron or kinetic theory
- Source term limited by `maxPhaseChangeSource` (W/m³)
- Applied only within `alphaMin` to `alphaMax` range (typically 0.01–0.99)
- Relaxation time controls explicit vs. implicit balance

### `/home/user/compInterFoam/twoTemperatureModel.H` (289 lines)
### `/home/user/compInterFoam/twoTemperatureModel.C` (2,206 lines)
**Electron-lattice temperature decoupling and coupling**

- **Purpose:** Solves coupled electron (Te) and lattice (Tl) temperature equations with laser absorption and gas-metal interfacial heat exchange
- **Key Data Members:**
  - `Te_`, `Tl_` - Electron and lattice temperature fields (K)
  - `Ce_` - Electron heat capacity (J/m³/K), optionally temperature-dependent
  - `Cl_` - Lattice heat capacity (J/m³/K), typically ~2.3e6 for titanium
  - `G_` - Electron-phonon coupling factor (W/m³/K), optionally temperature-dependent
  - `De_` - Electron thermal diffusivity (m²/s), ~10⁻⁴ for metals
  - `gasMetalExchangeCoeff_` - Gas-metal interfacial exchange (W/m³/K), ~5e17 for strong coupling
  - `ambientTemperature_` - Fallback temperature in non-metal cells (K)
  - `metalFraction_` - Reference to alpha1 (metal phase fraction)
  - Temperature bounds: `minTe`, `maxTe`, `minTl`, `maxTl` (K)
- **Key Methods:**
  - `correct()` - Main solve routine; updates Te and Tl via coupled equations
  - `solveElectronEquation()` - Assembles and solves electron energy equation
  - `solveLatticeEquation()` - Assembles and solves lattice energy equation
  - `checkEnergyConservation()` - Audits total energy balance
  - `gasMetalExchangeCoeffField()` - Returns spatially varying exchange coefficient
- **Configuration:**
  - Primary: `system/controlDict` → `twoTemperatureProperties` block
  - Alternate: `constant/twoTemperatureProperties` (optional dedicated file)
  - Key parameters: Ce, Cl, G, De, gasMetalExchangeCoeff, minTe, maxTe, minTl, maxTl
  - Sub-cycling: `electronSubCycles` for stiff coupling on small time-steps
- **Physics:**
  - Electron equation: `Ce·∂Te/∂t + ... = G(Tl - Te) + laser source + gas-metal exchange`
  - Lattice equation: `Cl·∂Tl/∂t + ... = G(Te - Tl) + phase-change sink + gas-metal exchange`
  - Gas-metal coupling: `h_gm·(Tl - T_gas)` where h_gm can use Kapitza resistance model

### `/home/user/compInterFoam/femtosecondLaserModel.H` (236 lines)
### `/home/user/compInterFoam/femtosecondLaserModel.C` (2,341 lines)
**Femtosecond laser energy deposition model**

- **Purpose:** Computes spatially and temporally resolved laser volumetric heat source for LIFT processes
- **Key Parameters:**
  - `peakIntensity_` - Peak laser intensity (W/m²), ~10¹⁶ for femtosecond LIFT
  - `pulseWidth_` - FWHM of pulse (s), typically 100–300 fs
  - `wavelength_` - Laser wavelength (m), e.g., 343 nm for UV on titanium
  - `spotSize_` - Gaussian beam diameter (m), typically 1–5 µm for LIFT
  - `absorptionCoeff_` - Material absorption coefficient (1/m), ~10⁸ for metal
  - `gasAbsorptionCoeff_` - Optional gas absorption (1/m)
  - `pulseEnergy_` - Total pulse energy (J), derived from intensity and spot
  - `direction_` - Unit vector of laser propagation direction
  - `focus_` - Beam focus point in domain (m)
  - `pulseFrequency_` - Pulse repetition rate (1/s), 0 for single pulse
  - `reflectivity_` - Surface reflectivity (0–1), ~0.35 for titanium at 343 nm
  - `incidenceAngle_` - Angle of incidence (rad) for Fresnel calculations
  - `continuousLaser_` - If true, treats as continuous wave rather than pulsed
  - `gaussianProfile_` - If true, uses Gaussian spatial profile; else flat-top
- **Key Methods:**
  - `correct()` - Updates laser source term based on current simulation time
  - `Q_laser()` - Returns volumetric heat source (W/m³)
  - `setFocusPosition()` - Moves focus for scanning/raster operations
  - `getAbsorbedPower()` - Returns total power absorbed by domain (W)
  - `validatePulseEnergy()` - Checks energy conservation within tolerance
- **Configuration:**
  - File: `constant/laserProperties`
  - Key tolerances: `pulseEnergyToleranceRel` (%), `pulseEnergyToleranceAbs` (J)
  - Temporal shape: Gaussian envelope with `pulseWidth` FWHM
  - Spatial shape: Gaussian or flat-top, optionally with transmission through metal film
  - Reflectivity → transmission: Fresnel corrections applied automatically

**Physics:**
- Absorptivity = 1 - reflectivity (Fresnel at interface)
- Penetration depth: δ = 1/absorptionCoeff
- Energy deposition: Q(r,t) = I(r,t) × absorptivity × exp(-z/δ) within spotSize
- Pulse envelope: Gaussian in time; transmitted or reflected at interfaces
- Continuous mode: Uniform time evolution

### `/home/user/compInterFoam/advancedInterfaceCapturing.H` (115 lines)
### `/home/user/compInterFoam/advancedInterfaceCapturing.C` (1,075 lines)
**Recoil pressure and interface stabilization**

- **Purpose:** Calculates evaporation-driven recoil pressure and applies safeguards to prevent numerical instabilities during interface ejection
- **Key Data Members:**
  - `recoilPressure_` - Volumetric recoil pressure field (Pa)
  - `previousRecoilPressure_` - Cached recoil pressure for temporal smoothing
  - `meltingTemp_`, `vaporTemp_` - Temperature activation thresholds
  - `recoilMax_` - Clamping ceiling on recoil pressure (typically 100 MPa)
  - `pressureScale_` - Scaling factor for evaporation-rate-to-pressure conversion (Pa·s)
  - `alphaMin_`, `alphaMax_` - Alpha range for recoil activation
  - `clampRecoil_` - Enable/disable pressure clamping
  - `recoilRelax_` - Temporal smoothing factor (0–1)
  - `boltzmannConstant_`, `vaporParticleMass_` - Kinetic theory parameters
  - `momentumAccommodationCoeff_` - Accommodation coefficient for momentum transfer
  - `rampProgress_`, `rampIncrement_` - Gradual ramp-up during start-up
- **Key Methods:**
  - `correct()` - Compute and apply recoil pressure
  - `calculateRecoilPressure()` - Core kinetic-theory recoil calculation
  - `recoilPressure()` - Access current recoil field
  - `write()` - Output diagnostics
- **Configuration:**
  - Primary: `system/controlDict` → `advancedInterfaceCapturing` block
  - Recoil model: Kinetic theory (Hertz-Knudsen + momentum accommodation)
  - Diagnostics: Detailed logging available with `verbose = true`

**Recoil Physics:**
- Evaporation rate drives recoil: P_recoil ≈ √(2πnM kT) × Γ_evap × accommodation
- Temperature-dependent activation: linear ramp from T_melt to T_vap
- Clamping: Prevents unphysical pressure spikes
- Ramping: Gradual ramp-in over first ~10–50 time-steps to avoid shock

---

## 4. SUPPORTING MODEL FILES

### `/home/user/compInterFoam/compressibleInterPhaseTransportModel.H` (106 lines)
### `/home/user/compInterFoam/compressibleInterPhaseTransportModel.C` (224 lines)
**Mixture turbulence and transport properties**

- **Purpose:** Manages turbulent transport closures for compressible two-phase flows
- **Key Functions:**
  - `divDevRhoReff(U)` - Divergence of stress tensor for momentum equation
  - Property blending: Mixture or per-phase turbulence model selection
  - Support for RAS (RANS) and LES turbulence models
- **Configuration:**
  - `constant/turbulenceProperties` - Turbulence model selection

### `/home/user/compInterFoam/createFields.H` (794 lines)
**Field initialization and management**

- **Purpose:** Initializes all solution fields at simulation start
- **Fields Created:**
  - Pressure: `p_rgh`, optionally `p`
  - Velocity: `U`
  - Phase fractions: `alpha1` (metal), `alpha2` (gas)
  - Temperatures: `T`, `Te`, `Tl` (from two-temperature model)
  - Density: `rho`, `rho1`, `rho2` (phase-specific)
  - Thermodynamic fields: `h` (enthalpy), `psi` (compressibility)
  - Turbulent fields: `k`, `epsilon`, etc. (if applicable)
  - Diagnostic fields: `Q_laser`, `phaseChangeSource`, `recoilPressure`
  - Time-stepping: `rDeltaT`, `rSubDeltaT` (for Local Time-Stepping, LTS)
- **Configuration:** Reads from:
  - `0.orig/` or `0/` directories for initial conditions
  - `constant/` for property dictionaries
  - `system/controlDict` for run parameters

### `/home/user/compInterFoam/compressibleAlphaEqnSubCycle.H` (202 lines)
**Sub-cycling helper for phase-fraction equation**

- **Purpose:** Enables multiple sub-iterations of alpha equation per PIMPLE step
- **Use Case:** Stabilizes interface on very small time-steps (femtosecond scale)

### `/home/user/compInterFoam/postProcess.H` (165 lines)
**Post-processing and diagnostic output**

- **Purpose:** Writes additional diagnostic fields to results
- **Outputs:**
  - LIFT process state snapshot (via `liftProcessTracker`)
  - Residual information
  - Phase-change and recoil diagnostics

### Supporting Header Files:
- **`alphaScheme.H`** (73 lines) - Interface compression scheme selection
- **`alphaSuSp.H`** (106 lines) - Phase-change source term assembly (Sp, Su)
- **`rhofs.H`** (100 lines) - Face density field for flux calculations

---

## 5. BUILD SYSTEM

### `/home/user/compInterFoam/Make/files`
```
twoPhaseMixtureThermo.C
femtosecondLaserModel.C
twoTemperatureModel.C
advancedInterfaceCapturing.C
compressibleInterPhaseTransportModel.C
compInterFoam.C

EXE = $(FOAM_USER_APPBIN)/compInterFoam
```

### `/home/user/compInterFoam/Make/options`
**Compiler flags and linked libraries**

- **Include paths (EXE_INC):**
  - Compressible transport models
  - Thermophysical models (basic)
  - Two-phase mixture properties
  - Interface properties
  - Turbulence models (both incompressible and compressible)
  - Finite volume discretization
  - Dynamic mesh support
  - fvOptions for source term handling

- **Link libraries (EXE_LIBS):**
  - `-lcompressibleTransportModels` - Compressible CFD closures
  - `-lfluidThermophysicalModels` - Thermophysical package
  - `-lspecie` - Species/mixture data
  - `-lturbulenceModels` - Turbulence model library
  - `-lcompressibleTurbulenceModels` - Compressible turbulence closures
  - `-ltwoPhaseProperties` - Two-phase mixture properties
  - `-linterfaceProperties` - Surface tension, VOF discretization
  - `-lfiniteVolume` - FV discretization
  - `-ldynamicFvMesh` - Mesh motion support
  - `-ldynamicMesh` - Dynamic mesh infrastructure
  - `-ltopoChangerFvMesh` - Mesh topology changes
  - `-lmeshTools` - Mesh utilities
  - `-lfvOptions` - fvOptions framework

**Build Command:**
```bash
source /opt/openfoam/etc/bashrc  # or similar OpenFOAM installation
cd /home/user/compInterFoam
wmake
```

---

## 6. CASE FILE STRUCTURE (EXAMPLE: TestCase)

### Case Directory Layout:
```
/home/user/compInterFoam/TestCase/
├── system/           # Solver configuration
├── constant/         # Static properties
├── 0.orig/          # Initial conditions (template)
└── 0/               # Initial conditions (working copy)
```

### `/home/user/compInterFoam/TestCase/system/`

#### **`controlDict`** (Primary Simulation Control)
- **Application & Timing:**
  - `application compInterFoam`
  - `startTime 0`, `endTime 2e-9` (2 nanoseconds typical)
  - `deltaT 1e-16` (0.1 fs initial), `maxDeltaT 5e-15` (5 fs maximum)
  - `adjustTimeStep yes` - Adaptive stepping enabled
  
- **Courant Number Limits:**
  - `maxCo 0.15` - Flow Courant number limit
  - `maxAlphaCo 0.25` - Interface Courant limit (stricter)
  - `maxThermalCourant 0.25` - Thermal Courant limit
  - `maxDi 8` - Interface dilation limit

- **Output Control:**
  - `writeControl adjustableRunTime`
  - `writeInterval 2e-13` (200 fs between outputs)
  - `purgeWrite 200` - Keep most recent 200 time-steps
  - `writeFormat binary`, `writePrecision 10`
  - `fileHandler uncollated` - Do not collect parallel files

- **Model Configuration Blocks:**
  
  **`phaseChangeCoeffs`:**
  - `model clausius_clapeyron` - Evaporation model
  - `hf [0 2 -2 0 0 0 0] 9.1e6` - Latent heat (J/kg) for Ti
  - `gasConstant [0 2 -2 -1 0 0 0] 174` - Hertz-Knudsen constant (J/(kg·K))
  - `Tsol 1941 K`, `Tvap 2200 K` - Phase-change thresholds
  - `evaporationCoeff 0.3` - Evaporation coefficient
  - `alphaMin 0.001`, `alphaMax 0.999` - Phase-fraction window
  - `maxSource [1 -1 -3 0 0 0 0] 1e22` - Max evaporation source (W/m³)
  - `relaxationTime 1e-12` - Explicit source relaxation
  - `activationTime ((0 2e-10))` - Time windows when phase-change is active

  **`twoTemperatureProperties`:**
  - `Ce` - Electron heat capacity (J/m³/K), e.g., 210 or Function1
  - `Cl` - Lattice heat capacity (J/m³/K), e.g., 2.3e6
  - `G` - Electron-phonon coupling (W/m³/K), e.g., 5e17
  - `De` - Electron thermal diffusivity (m²/s), e.g., 1e-4
  - `gasMetalExchangeCoeff` - Gas-metal coupling (W/m³/K), e.g., 5e17
  - `minTe`, `maxTe`, `minTl`, `maxTl` - Temperature bounds (K)
  - `metalFractionFloor` - Minimum alpha for activation
  - `electronSubCycles` - Number of electron sub-iterations per step
  
  **`advancedInterfaceCapturing`:**
  - `meltingTemperature`, `vaporTemperature` - Thresholds
  - `recoilMax 1e8` (100 MPa) - Pressure ceiling
  - `clampRecoil true` - Enable clamping
  - `recoilRelax 1.0` - Temporal relaxation
  - `alphaMin 0.001`, `alphaMax 0.999` - Active interface range

#### **`fvSolution`** (Linear Solver & Iteration Settings)
```
solvers:
  alpha.*      - PBiCGStab + DILU preconditioner
               - tolerance 1e-12, relTol 0
               - nAlphaCorr 5, nAlphaSubCycles 5
  p_rgh        - GAMG with DIC smoother
               - tolerance 1e-6, relTol 0.05
  U.*          - smoothSolver + symGaussSeidel
               - tolerance 1e-6, relTol 0.5
  T / Te / Tl  - smoothSolver + symGaussSeidel
               - tolerance 1e-8, relTol 0.01
  
PIMPLE:
  momentumPredictor true
  nOuterCorrectors 10
  nCorrectors 3
  nNonOrthogonalCorrectors 2
  p_rghRefCell 0, p_rghRefValue 101325 Pa
  
relaxationFactors:
  fields: p_rgh 0.3, alpha.* 0.5
  equations: U 0.5, T 0.5, Te 0.5, Tl 0.5
```

#### **`fvSchemes`** (Discretization Schemes)
- **Time Derivative (`ddtSchemes`):** Euler (first-order, stable)
- **Gradients:** Gauss linear with limiting
- **Divergence:**
  - Alpha: `Gauss vanLeer` + `interfaceCompression` for compression flux
  - Momentum: `Gauss upwind`
  - Energy: `Gauss upwind`
  - Pressure: `Gauss upwind`
- **Laplacian:** `Gauss linear corrected` / `orthogonal` for thermal diffusion
- **Interpolation:** Linear
- **Surface Gradient:** Corrected
- **Wall Distance:** meshWave method

#### **`blockMeshDict`**
Mesh generation for domain:
- Typical LIFT domain: x ∈ [0, 50 µm], y ∈ [0, 28.07 µm], z ∈ [0, 10 µm]
- 3D structured hex mesh
- Refined near film interface and laser focus region

#### **`decomposeParDict`**, **`refineMeshDict`**, **`setFieldsDict`**, **`topoSetDict`**
Domain decomposition, mesh refinement, field initialization, and topology tools

### `/home/user/compInterFoam/TestCase/constant/`

#### **`thermophysicalProperties`** (Global Mixture Data)
```
phases (metal air)
sigma 1.64 [N/m]  # Surface tension at metal-air interface
thermoType: hePsiThermo + pureMixture + const transport + hConst thermo + perfectGas
mixture: default air properties
```

#### **`thermophysicalProperties.metal`** (Metal Phase Properties)
```
thermoType: heRhoThermo + Boussinesq
mixture:
  specie: molWeight 47.867 kg/kmol (titanium)
  thermodynamics: Cp 650 J/(kg·K), Href 0
  transport: mu 2.35e-3 Pa·s (molten), Pr 0.032
  equationOfState: rho0 4515 kg/m³, beta 7.6e-5 1/K, T0 300 K
LIFT Properties:
  Tsol 1941 K (melting), Tliq 1941 K, Tvap 3560 K
  hf 9.1e6 J/kg (latent heat)
  kappa 17.2 W/(m·K) (thermal conductivity)
```

#### **`thermophysicalProperties.air`** (Gas Phase Properties)
```
thermoType: hePsiThermo + pureMixture + const + hConst + perfectGas
mixture:
  specie: molWeight 28.97 kg/kmol
  thermodynamics: Cp 1000 J/(kg·K)
  transport: mu 1.8e-5 Pa·s, Pr 0.7
  equationOfState: R 287 J/(kg·K)
```

#### **`laserProperties`** (Femtosecond Laser Model)
```
Laser Parameters:
  laserStartTime 0, laserEndTime 2e-10 (200 ps)
  filmThicknessExpected 71.4e-9 m (71.4 nm Ti film)
  filmCenterY 28.0357e-6 m (focus location)
  
Energy & Pulse:
  pulseEnergy 2.0e-8 J (20 nJ)
  pulseWidth 200e-15 s (200 fs FWHM)
  wavelength 343e-9 m (343 nm, UV)
  spotSize 3.2e-6 m (3.2 µm diameter)
  
Geometry:
  focus (25e-6, 28.0357e-6, 5e-6) m - beam center
  direction (0, -1, 0) - downward through substrate
  
Material Interaction:
  absorptionCoeff 1.03e8 1/m (Ti at 343 nm)
  reflectivity 0.35 (45% absorption)
  gasAbsorptionCoeff 0 (disabled)
  
Profile:
  gaussianProfile true
  continuousLaser false (pulsed)
  maxVolumetricSource 7e24 W/m³
  pulseFrequency 0 (single pulse)
  
Energy Validation:
  pulseEnergyToleranceRel 0.05 (5%)
  pulseEnergyToleranceAbs 5e-10 J (0.5 nJ)
```

#### **`transportProperties`** (Mixture Transport)
- Kinematic viscosities: `nu1` (metal), `nu2` (gas)
- Densities: `rho1`, `rho2`

#### **`turbulenceProperties`**
- Turbulence model: Typically `laminar` for LIFT or `RASProperties`/`LESProperties`

#### **`fvOptions`** (Optional source terms, e.g., gravity)
```
g (0, -9.81, 0) m/s²
```

### `/home/user/compInterFoam/TestCase/0.orig/` (Initial Conditions Template)

#### **Field Files (all with uniform initial values + boundary conditions):**
- **`p`** (101325 Pa) - Pressure
- **`p_rgh`** (0 Pa) - Hydrostatic reference pressure
- **`U`** (0 m/s) - Velocity
- **`T`** (300 K) - Gas temperature (if used)
- **`Te`** (300 K) - Electron temperature
- **`Tl`** (300 K) - Lattice temperature
- **`alpha.metal`** (1 in donor film, 0 in air) - Metal phase fraction
- **`alpha.air`** (0 in donor film, 1 in air) - Air phase fraction

Each file specifies:
- Internal field values (uniform or spatially varying)
- Boundary conditions (typically zero-gradient, fixed value, or cyclic)

---

## 7. TEST CASES / EXAMPLE SIMULATIONS

### Directory Structure:
```
TEST1/     - Baseline LIFT simulation (full 2 ns run)
TEST2/     - Alternative configuration (different laser/material parameters)
TestCase/  - Primary documented example (described in detail above)
```

### Common Features:
- Identical system/ and constant/ structure as TestCase
- Same mesh generation procedure (blockMeshDict-based)
- Variations in laser properties or phase-change coefficients for sensitivity studies

---

## 8. DOCUMENTATION (Models/ Directory)

### `/home/user/compInterFoam/Models/`

Comprehensive technical reports (markdown format):
- **`compInterFoam_Implementation_Report.md`** - Solver architecture overview
- **`FemtosecondLaserModel_Implementation_Report.md`** - Laser physics detail
- **`TwoTemperatureModel_Detailed_Report.md`** - Electron-lattice coupling
- **`AdvancedInterfaceCapturing_Thesis_Report.md`** - Recoil pressure & interface stabilization
- **`twoPhaseMixtureThermo_Technical_Report.md`** - Phase-change & thermodynamics
- **`compressibleInterPhaseTransportModel_Thesis_Report.md`** - Transport closure details

---

## 9. KEY PHYSICS SUMMARY

### Coupled Multiphysics Stack:

```
compInterFoam (Main Solver)
  ├── Momentum (UEqn.H)
  │   └─ Recoil Pressure (advancedInterfaceCapturing)
  │
  ├── Pressure-Velocity Coupling (pEqn.H) via PIMPLE
  │   └─ Reference cell at domain inlet
  │
  ├── Phase Fraction / Interface (alphaEqn.H)
  │   ├─ VOF interface tracking with MULES
  │   └─ Phase-change mass source (twoPhaseMixtureThermo)
  │
  ├── Gas Temperature (TEqn.H)
  │   ├─ Laser heating (femtosecondLaserModel)
  │   └─ Gas-metal exchange (twoTemperatureModel)
  │
  ├── Two-Temperature Physics (twoTemperatureModel)
  │   ├─ Electron Temperature (Te equation)
  │   ├─ Lattice Temperature (Tl equation)
  │   └─ Electron-phonon coupling (G term)
  │
  └── Phase-Change Modeling (twoPhaseMixtureThermo)
      ├─ Evaporation (Hertz-Knudsen)
      └─ Recoil Pressure (kinetic theory)
```

### Time-Stepping Strategy:
- **Outer loop:** PIMPLE iterations (pressure-velocity correction)
- **Inner loop (per PIMPLE):**
  - Solve momentum equation (UEqn.H)
  - Pressure-velocity coupling (pEqn.H, nCorrectors iterations)
  - Solve phase-fraction (alphaEqn.H, with sub-cycling)
  - Update thermodynamic properties (mixture.correct())
  - Solve energy equations (TEqn.H, with two-temperature coupling)
- **Adaptive control:** Courant number limits scale dt to maintain stability

---

## 10. DIRECTORY TREE SUMMARY

```
/home/user/compInterFoam/
├── README.md                                  # User guide & configuration reference
├── Thesis_MinhTriNguyen_10754416.pdf         # Original research thesis
├── EXPERIMENTAL_REPLICATE_ANALYSIS.md        # Comparison with literature
│
├── Make/
│   ├── files        # Source file listing
│   └── options      # Compiler flags & libraries
│
├── Models/          # Technical documentation (markdown reports)
│   ├── compInterFoam_Implementation_Report.md
│   ├── FemtosecondLaserModel_Implementation_Report.md
│   ├── TwoTemperatureModel_Detailed_Report.md
│   ├── AdvancedInterfaceCapturing_Thesis_Report.md
│   ├── twoPhaseMixtureThermo_Technical_Report.md
│   └── compressibleInterPhaseTransportModel_Thesis_Report.md
│
├── SOURCE FILES (Main Solver):
│   ├── compInterFoam.C                        # Main solver executable
│   ├── createFields.H                         # Field initialization
│   ├── postProcess.H                          # Output & diagnostics
│   │
│   ├── EQUATION FILES:
│   │   ├── UEqn.H                             # Momentum equation
│   │   ├── pEqn.H                             # Pressure-correction
│   │   ├── alphaEqn.H                         # Phase-fraction (VOF)
│   │   ├── TEqn.H                             # Temperature equation
│   │   ├── alphaScheme.H                      # Interface scheme
│   │   ├── alphaSuSp.H                        # Phase-change source
│   │   └── rhofs.H                            # Face density
│   │
│   └── PHYSICS MODELS:
│       ├── twoPhaseMixtureThermo.H/.C         # Mixture thermo & phase-change
│       ├── twoTemperatureModel.H/.C           # Electron-lattice coupling
│       ├── femtosecondLaserModel.H/.C         # Laser heat source
│       ├── advancedInterfaceCapturing.H/.C    # Recoil pressure & interface
│       ├── compressibleInterPhaseTransportModel.H/.C  # Turbulence/transport
│       └── compressibleAlphaEqnSubCycle.H     # Sub-cycling helper
│
├── TEST CASES:
│   ├── TestCase/    # Primary documented example
│   │   ├── system/
│   │   │   ├── controlDict              # Simulation control & time-stepping
│   │   │   ├── fvSolution               # Linear solver & PIMPLE settings
│   │   │   ├── fvSchemes                # Discretization schemes
│   │   │   ├── blockMeshDict            # Mesh generation
│   │   │   ├── decomposeParDict         # Parallelization
│   │   │   ├── setFieldsDict            # Initial field setup
│   │   │   ├── topoSetDict              # Topology definition
│   │   │   └── refineMeshDict           # Mesh refinement
│   │   ├── constant/
│   │   │   ├── thermophysicalProperties           # Global thermo data
│   │   │   ├── thermophysicalProperties.metal     # Metal phase (Ti)
│   │   │   ├── thermophysicalProperties.air       # Gas phase (air)
│   │   │   ├── laserProperties                    # Femtosecond laser config
│   │   │   ├── transportProperties                # Transport data
│   │   │   ├── turbulenceProperties               # Turbulence model
│   │   │   ├── fvOptions                          # Optional sources
│   │   │   └── dynamicMeshDict                    # Dynamic mesh control
│   │   └── 0.orig/
│   │       ├── p, p_rgh, U               # Pressure & velocity IC
│   │       ├── T, Te, Tl                 # Temperature fields IC
│   │       └── alpha.metal, alpha.air    # Phase-fraction IC
│   │
│   ├── TEST1/       # Baseline simulation
│   │   └── [same structure as TestCase]
│   │
│   └── TEST2/       # Alternative configuration
│       └── [same structure as TestCase]
│
├── meshConvergenceStudy/  # Mesh sensitivity study
│
└── .git/            # Git repository metadata
```

---

## 11. KEY FILE PATHS (QUICK REFERENCE)

| Component | File Path | Purpose |
|-----------|-----------|---------|
| **Main Solver** | `/home/user/compInterFoam/compInterFoam.C` | Entry point & time loop |
| **Momentum** | `/home/user/compInterFoam/UEqn.H` | Momentum equation |
| **Pressure** | `/home/user/compInterFoam/pEqn.H` | Pressure-velocity coupling |
| **Interface** | `/home/user/compInterFoam/alphaEqn.H` | Phase fraction evolution |
| **Temperature** | `/home/user/compInterFoam/TEqn.H` | Thermal energy equation |
| **Thermodynamics** | `/home/user/compInterFoam/twoPhaseMixtureThermo.C` | Phase-change & latent heat |
| **Two-Temperature** | `/home/user/compInterFoam/twoTemperatureModel.C` | Te-Tl coupling |
| **Laser** | `/home/user/compInterFoam/femtosecondLaserModel.C` | Laser heat source |
| **Recoil** | `/home/user/compInterFoam/advancedInterfaceCapturing.C` | Recoil pressure |
| **Build Files** | `/home/user/compInterFoam/Make/` | Compilation config |
| **Example Case** | `/home/user/compInterFoam/TestCase/` | Template simulation |
| **Simulation Control** | `/home/user/compInterFoam/TestCase/system/controlDict` | Time-stepping & models |
| **Solvers Config** | `/home/user/compInterFoam/TestCase/system/fvSolution` | PIMPLE & linear solvers |
| **Schemes** | `/home/user/compInterFoam/TestCase/system/fvSchemes` | Discretization schemes |
| **Laser Config** | `/home/user/compInterFoam/TestCase/constant/laserProperties` | Laser parameters |
| **Material Props** | `/home/user/compInterFoam/TestCase/constant/thermophysicalProperties.metal` | Metal thermophysics |
| **Initial Conditions** | `/home/user/compInterFoam/TestCase/0.orig/` | Starting fields |

---

## 12. TYPICAL SIMULATION WORKFLOW

1. **Setup:**
   - Copy TestCase → MyCase
   - Modify `constant/laserProperties` (energy, spotSize, focus)
   - Adjust `system/controlDict` phase-change & two-temp parameters
   - Tune mesh in `system/blockMeshDict` if needed

2. **Mesh Generation:**
   ```bash
   cd MyCase
   blockMesh
   ```

3. **Field Initialization:**
   ```bash
   setFields
   ```

4. **Solver Compilation (once):**
   ```bash
   cd /home/user/compInterFoam
   source /opt/openfoam/etc/bashrc
   wmake
   ```

5. **Run Simulation:**
   ```bash
   cd MyCase
   compInterFoam > log.out 2>&1 &
   # Monitor with: tail -f log.out
   ```

6. **Post-Processing:**
   - Output in `MyCase/postProcessing/` and time directories (0.1e-10, etc.)
   - ParaView/VisIt visualization of fields: p, U, T, Te, Tl, alpha.metal, recoilPressure
   - Python scripts to extract LIFT diagnostics (max Te/Tl, recoil vs time, material loss)

---

## 13. CRITICAL CONFIGURATION NOTES

### For Femtosecond LIFT on Titanium:
- **Laser:** 343 nm, 200 fs, ~20 nJ, 3.2 µm spot
- **Material:** Titanium (MW 47.9, Cp 650 J/(kg·K), ρ 4515 kg/m³, Tvap 3560 K)
- **Time-step:** dt ~ 0.1–5 fs (1e-16 to 5e-15 s)
- **Mesh:** ~100 nm cells in laser region, 1 µm far-field
- **Phase-change:** Clasius-Clapeyron, αmin 0.001, αmax 0.999
- **Two-temperature:** Ce~210, Cl~2.3e6, G~5e17, De~1e-4 (J/m³/K and W/m³/K)
- **Recoil:** 50–100 MPa max pressure, temperature-ramped from Tmelt to Tvap

---

**Document Generated:** 2025-11-16
**Codebase Analyzed:** compInterFoam OpenFOAM Solver with LIFT Physics Extensions
