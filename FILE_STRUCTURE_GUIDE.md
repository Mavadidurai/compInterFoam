# compInterFoam File Structure & Key Configuration Files

## Project Overview

```
/home/user/compInterFoam/
├── Make/                              # Build system
│   ├── files
│   └── options
│
├── TEST1/                             # PRIMARY LIFT SIMULATION CASE
│   ├── constant/                      # MATERIAL PROPERTIES & CONSTANTS
│   │   ├── thermophysicalProperties     (main thermo type: hePsiThermo + perfectGas)
│   │   ├── thermophysicalProperties.metal  ⚠️ ISSUE: hf mislabeled
│   │   ├── thermophysicalProperties.air    (perfect gas for air)
│   │   ├── transportProperties            (viscosity, surface tension)
│   │   ├── laserProperties                 (30 nJ, 343 nm, 3.2 µm spot)
│   │   ├── fvOptions                      (source terms)
│   │   ├── dynamicMeshDict                (mesh motion if enabled)
│   │   ├── turbulenceProperties           (laminar)
│   │   └── g                              (gravity vector)
│   │
│   ├── system/                        # NUMERICAL SOLVER SETTINGS
│   │   ├── controlDict                   (time stepping: 1fs base, 2ns simulation)
│   │   ├── fvSchemes                     (Euler ddt, vanLeer alpha, etc.)
│   │   ├── fvSolution                    (PIMPLE, GAMG, relaxation factors)
│   │   ├── blockMeshDict                 (mesh: 50µm × 28.07µm × 10µm)
│   │   ├── topoSetDict                   (region selection)
│   │   ├── setFieldsDict                 (initial phase fractions)
│   │   └── refineMeshDict                (optional mesh refinement)
│   │
│   └── 0.orig/                        # INITIAL & BOUNDARY CONDITIONS
│       ├── U                            (velocity: no-slip walls)
│       ├── p                            (pressure)
│       ├── p_rgh                        (pressure - hydrostatic)
│       ├── T                            (temperature: 300 K ambient)
│       ├── Te                           (electron temperature: 300 K)
│       ├── Tl                           (lattice temperature: 300 K)
│       ├── alpha.metal                  (metal volume fraction: 0 initially)
│       ├── alpha.air                    (air volume fraction: 1 initially)
│       └── laserProperties              (time-varying laser source)
│
├── RealisticLIFT/                     # ALTERNATIVE GEOMETRY CASE
│   ├── constant/                      (same structure, different laserProperties)
│   ├── system/                        (same settings)
│   └── 0.orig/                        (same IC/BC)
│
├── Source Code Files (C++):
│   ├── compInterFoam.C                (main solver)
│   ├── createFields.H                 (field initialization)
│   ├── UEqn.H                         (momentum equation)
│   ├── pEqn.H                         (pressure equation)
│   ├── alphaEqn.H                     (phase fraction transport)
│   ├── TEqn.H                         (temperature equation)
│   ├── rhofs.H                        (density field management)
│   ├── postProcess.H                  (output diagnostics)
│   ├── femtosecondLaserModel.H/C      (laser source term)
│   ├── twoTemperatureModel.H/C        (electron-phonon coupling)
│   ├── twoPhaseMixtureThermo.H/C      (thermodynamic mixture)
│   └── advancedInterfaceCapturing.H/C (recoil pressure)
│
├── README.md                          (comprehensive documentation)
├── Thesis_MinhTriNguyen_10754416.pdf  (research thesis)
└── meshConvergenceStudy/              (parametric mesh study scripts)
```

---

## Critical Configuration Files for LIFT Physics

### 1. Laser Parameters: `constant/laserProperties`

**Location**: 
- TEST1: `/home/user/compInterFoam/TEST1/constant/laserProperties`
- RealisticLIFT: `/home/user/compInterFoam/RealisticLIFT/constant/laserProperties`

**Key Parameters**:
```
pulseEnergy              [J]        30 nJ (TEST1) vs 60 nJ (RealisticLIFT)
pulseWidth              [s]        200 fs
wavelength              [m]        343 nm (UV-A)
spotSize                [m]        3.2 µm (TEST1) vs 6.0 µm (RealisticLIFT)
absorptionCoeff         [1/m]      6e7 (constant, no T-dependence)
reflectivity            [-]        0.35
direction               (vector)   (0 1 0) TEST1 vs (0 -1 0) RealisticLIFT
focus                   (point)    Center of Ti film at each case
gasAbsorptionCoeff      [1/m]      0 (air doesn't absorb)
maxVolumetricSource     [W/m³]     7e24 (cap on power density)
laserStartTime          [s]        0
laserEndTime            [s]        2e-10 (200 ps)
```

### 2. Material Properties: `constant/thermophysicalProperties.metal`

**Location**:
- TEST1: `/home/user/compInterFoam/TEST1/constant/thermophysicalProperties.metal`
- RealisticLIFT: `/home/user/compInterFoam/RealisticLIFT/constant/thermophysicalProperties.metal`

**Titanium Properties**:
```
TEST1:
  ρ = 4515 kg/m³ (constant)           ⚠️ Should be T-dependent
  Cp = 560 J/kg·K                     ⚠️ Should vary with phase
  κ = 17.2 W/m·K                      (thermal conductivity)
  μ = 2.25e-3 Pa·s                    (dynamic viscosity)
  Pr = 0.03                           (Prandtl number)
  hf = 3.65e5 J/kg                    🚨 MISLABELED (is fusion, not vapor)
  Tsol/Tliq = 1941 K
  Tvap = 3560 K

RealisticLIFT:
  ρ(T) via Boussinesq: ρ = ρ0(1 - β(T-T0))    ✅ BETTER
    ρ0 = 4515 kg/m³
    β = 7.6e-5 K⁻¹ (thermal expansion)
    T0 = 300 K
  Cp = 650 J/kg·K                     ✅ Higher (liquid phase)
  hf = 9.1e6 J/kg                     ✅ CORRECT vaporization heat
```

### 3. Time Stepping Control: `system/controlDict`

**Location**: `/home/user/compInterFoam/TEST1/system/controlDict` (same for both cases)

**Critical Time Parameters**:
```
startTime               0
endTime                 2e-9          (2 ns total simulation)
deltaT                  1e-12         (1 fs base time step)
adjustTimeStep          yes           (adaptive stepping enabled)
maxCo                   0.5           (Courant number limit)
maxDeltaT               1e-11         (10 fs max growth)
minDeltaT               1e-14         (0.01 fs minimum)
writeInterval           1e-10         (1 ps between outputs)
```

### 4. Phase Change Model: `system/controlDict` → `phaseChangeCoeffs`

**Location**: `/home/user/compInterFoam/TEST1/system/controlDict` (lines 51-73)

**Key Settings**:
```
model                   clausius_clapeyron
hf                      9.1e6 J/kg            (vaporization enthalpy)
gasConstant             174 J/kg/K            (Ti vapor R)
Tsol/Tvap               1941 K / 3560 K
evaporationCoeff        0.03                  (accommodation)
relaxationTime          1e-11 s               (10 ps response time)
maxSource               1e22 W/m³
alphaMin                0.001                 (allow phase change in thin cells)
onlyAboveVapor          false (TEST1) / true (RealisticLIFT)
activationTime          ((0 2e-10))           (active 0-200 ps)
```

### 5. Two-Temperature Model: `system/controlDict` → `twoTemperatureProperties`

**Location**: `/home/user/compInterFoam/TEST1/system/controlDict` (lines 119-195)

**Key Settings**:
```
Ce (electron heat capacity)
  Ce = γ·Te = 630·Te [J/m³/K²]

G (electron-phonon coupling) - Temperature lookup table:
  300 K:    1.0e18 W/m³/K
  1000 K:   3.0e18 W/m³/K
  3000 K:   1.0e19 W/m³/K
  5000 K:   2.0e19 W/m³/K
  10000 K:  5.0e19 W/m³/K
  20000 K:  1.0e20 W/m³/K

De (electron diffusivity)    1e-4 m²/s
Cl (lattice heat capacity)   2.5e6 J/m³/K
maxTe                        20000 K
maxTl                        8000 K (TEST1) / 10000 K (RealisticLIFT)

gasMetalExchangeCoeff
  type        kapitza         (acoustic mismatch model)
  Z_metal     2.3e7 Pa·s/m    (Ti acoustic impedance)
  Z_gas       383 Pa·s/m      (Ar acoustic impedance)
```

### 6. Recoil Pressure Model: `system/controlDict` → `advancedInterfaceCapturing`

**Location**: `/home/user/compInterFoam/TEST1/system/controlDict` (lines 84-111)

**Key Settings**:
```
model                   kinetic_theory        (recoil from evaporation)
stickingCoeff           0.18
momentumAccom           0.18
recoilMax               3.0e9 Pa              (3 GPa pressure cap)
clampRecoil             false                 (allow natural saturation)
scaleRecoilMax          false
recoilRelax             0.5                   (temporal relaxation)
alphaMin                0.001                 (minimum metal fraction)
maxPhysicalTemperature  10000 K
maxRecoilPressure       3.0e9 Pa
```

### 7. Numerical Solver Settings: `system/fvSolution`

**Location**: `/home/user/compInterFoam/TEST1/system/fvSolution`

**PIMPLE Coupling** (lines 118-155):
```
nOuterCorrectors        3               (outer pressure-velocity loops)
nCorrectors             3               (inner velocity corrections)
nNonOrthogonalCorrectors 2              (non-orthogonal mesh corrections)
momentumPredictor       true

residualControl:
  alpha.*     → 1e-9 (strict for interface)
  p_rgh       → 1e-4 (pressure convergence)
  U           → 1e-5 (momentum convergence)
  Te, Tl      → 5e-9 (thermal convergence)
```

**Solver Details** (lines 18-115):
```
alpha.metal   → PBiCGStab, tolerance 1e-12, nAlphaSubCycles 5
p_rgh         → GAMG, tolerance 1e-6
U             → smoothSolver, Pr=0.5, maxIter=200
T, Te, Tl     → smoothSolver, tolerance 1e-8, maxIter=1000
```

### 8. Discretization Schemes: `system/fvSchemes`

**Location**: `/home/user/compInterFoam/TEST1/system/fvSchemes`

**Key Schemes**:
```
Time (ddt):
  default           Euler               (1st-order implicit)

Gradient (grad):
  default           Gauss linear
  grad(U)           cellLimited Gauss linear 0.5    (slope limiter)
  grad(alpha.metal) cellLimited Gauss linear 1.0

Divergence (div):
  div(phi,alpha)         Gauss vanLeer        (TVD for interface)
  div(phir,alpha)        Gauss interfaceCompression  (compression flux)
  div(rhoPhi,U)          Gauss linearUpwind grad(U)
  div(phi,T)             Gauss upwind         (conservative for thermal)
  div(phi,Te)            Gauss upwind
  div(phi,Tl)            Gauss upwind

Laplacian (laplacian):
  default                Gauss linear corrected
  laplacian(alphaEff,T)  Gauss linear orthogonal  (thermal diffusion)
  laplacian(ke,Te)       Gauss linear orthogonal
  laplacian(kl,Tl)       Gauss linear orthogonal
```

### 9. Mesh Definition: `system/blockMeshDict`

**Location**: `/home/user/compInterFoam/TEST1/system/blockMeshDict`

**Domain Structure**:
```
TEST1 Geometry (bottom to top):
  y0 = 0.0000 µm      Donor substrate bottom
  y1 = 8.0000 µm      Donor substrate top / Ti film bottom
  y2 = 8.0714 µm      Ti film top / air gap bottom  [71.4 nm thick]
  y3 = 20.0714 µm     Air gap top / receiver bottom [12 µm gap]
  y4 = 28.0714 µm     Receiver substrate top [8 µm thick]

Mesh Resolution:
  X: 40 cells × 250 nm/cell = 10 µm
  Y: Variable:
    - Donor: 200 cells / 8 µm = 40 nm/cell
    - Ti:    200 cells / 71.4 nm ≈ 0.357 nm/cell ✅ Excellent
    - Gap:   200 cells / 12 µm = 60 nm/cell
    - Receiver: 200 cells / 8 µm = 40 nm/cell
  Z: 40 cells × 250 nm/cell = 10 µm

Total cells: ~12.8 million (large for fs simulations)
```

### 10. Initial & Boundary Conditions: `0.orig/` Files

**Velocity** (`0.orig/U`):
```
internalField: uniform (0 0 0)
donorSubstrate: fixedValue (0 0 0)    [no-slip wall]
receiver:       zeroGradient          [free surface]
Sides:          symmetryPlane         [reduce domain]
```

**Temperature** (`0.orig/T`, `0.orig/Te`, `0.orig/Tl`):
```
internalField:  uniform 300 K         [ambient]
donorSubstrate: fixedValue 300 K      [heat sink]
receiver:       zeroGradient          [insulated]
Sides:          symmetryPlane
```

**Phase Fractions** (`0.orig/alpha.metal`, `0.orig/alpha.air`):
```
internalField:  depends on layer:
                [0] if air gap
                [1] if Ti film in Ti film region
                [0] elsewhere (substrate initially)
Walls:          zeroGradient
Sides:          symmetryPlane
```

---

## How to Use This Information

### For Running Simulations:
1. **Copy TEST1 directory** and modify `constant/laserProperties` for your parameters
2. **Use RealisticLIFT geometry** if laser enters from top
3. **Check times in controlDict** - `laserEndTime` should cover full pulse + relaxation

### For Post-Processing:
1. Monitor **recoil pressure time history** in output logs
2. Extract **ejection velocity** from U field at later times
3. Check **energy balance** - sum of absorbed energy should match temperature rise

### For Validation:
1. Compare recoil pressure peak to **~80 MPa** (literature value)
2. Verify fluence: 
   - TEST1: 0.38 J/cm² (marginal)
   - RealisticLIFT: 0.2 J/cm² (optimal)
3. Check that Tvap penetration depth matches ablation depth observed

---

## Critical File Absolute Paths

For reference in automation scripts:

```
/home/user/compInterFoam/TEST1/constant/laserProperties
/home/user/compInterFoam/TEST1/constant/thermophysicalProperties.metal
/home/user/compInterFoam/TEST1/constant/transportProperties
/home/user/compInterFoam/TEST1/system/controlDict
/home/user/compInterFoam/TEST1/system/fvSolution
/home/user/compInterFoam/TEST1/system/fvSchemes
/home/user/compInterFoam/TEST1/system/blockMeshDict
/home/user/compInterFoam/TEST1/0.orig/U
/home/user/compInterFoam/TEST1/0.orig/T
/home/user/compInterFoam/TEST1/0.orig/Te
/home/user/compInterFoam/TEST1/0.orig/Tl

/home/user/compInterFoam/RealisticLIFT/constant/laserProperties
/home/user/compInterFoam/RealisticLIFT/constant/thermophysicalProperties.metal
/home/user/compInterFoam/RealisticLIFT/system/controlDict
/home/user/compInterFoam/RealisticLIFT/system/blockMeshDict
```

