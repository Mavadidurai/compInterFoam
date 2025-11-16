# compInterFoam Codebase - Documentation Index

This directory now contains comprehensive documentation to help you understand the **compInterFoam** Laser-Induced Forward Transfer (LIFT) solver implementation.

## Documentation Files

### 1. **CODEBASE_MAP.md** (792 lines, 34 KB)
**Comprehensive reference guide covering:**
- Main solver file and physics models overview
- Detailed description of each equation file (UEqn, pEqn, alphaEqn, TEqn)
- Full physics model implementation details:
  - Two-phase mixture thermodynamics & phase-change
  - Two-temperature electron-lattice coupling
  - Femtosecond laser energy deposition
  - Advanced interface capturing & recoil pressure
  - Transport closures and supporting files
- Build system configuration
- Complete case file structure with all configuration parameters
- Test case examples
- Key physics summary and time-stepping hierarchy
- Critical configuration notes for LIFT on titanium

**Best for:** In-depth understanding of every component, detailed reference lookups

---

### 2. **QUICK_REFERENCE.txt** (230 lines, 11 KB)
**Quick lookup guide with:**
- Main solver entry point location
- All coupled equation files with one-line descriptions
- Core physics models with file locations and key parameters
- Supporting files and build system overview
- Case configuration structure summary
- Key algorithms and formulations (concise)
- Time-stepping hierarchy diagram
- Critical file paths for easy navigation
- Code statistics (line counts, documentation)

**Best for:** Fast navigation, finding specific files, understanding overall structure

---

### 3. **ARCHITECTURE_DIAGRAM.txt** (254 lines, 13 KB)
**Visual representation of:**
- Complete main solver loop flow diagram
- PIMPLE iteration structure with all coupled equations
- Data flow dependencies between models
- Configuration hierarchy and parameter organization

**Best for:** Understanding execution flow, code dependencies, parameter propagation

---

## Source Code Files (14,594 lines total)

### Main Solver Entry Point
- **`compInterFoam.C`** (1,243 lines) - Main time loop & orchestration

### Equation Files
- **`UEqn.H`** (654 lines) - Momentum equation
- **`pEqn.H`** (517 lines) - Pressure-velocity coupling (PIMPLE)
- **`alphaEqn.H`** (374 lines) - Phase-fraction transport (VOF + MULES)
- **`TEqn.H`** (1,431 lines) - Temperature equation with multi-physics coupling

### Core Physics Models
- **`twoPhaseMixtureThermo.H/C`** (451+2,037 lines) - Phase-change & evaporation
- **`twoTemperatureModel.H/C`** (289+2,206 lines) - Electron-lattice coupling
- **`femtosecondLaserModel.H/C`** (236+2,341 lines) - Laser heat source
- **`advancedInterfaceCapturing.H/C`** (115+1,075 lines) - Recoil pressure

### Supporting Files
- **`createFields.H`** (794 lines) - Field initialization
- **`compressibleInterPhaseTransportModel.H/C`** (106+224 lines) - Transport closures
- **`postProcess.H`** (165 lines) - Diagnostics output
- **`alphaScheme.H`**, **`alphaSuSp.H`**, **`rhofs.H`** - Helper files

---

## Configuration Files (Example: TestCase)

### system/ Directory
- **`controlDict`** - Primary simulation control, time-stepping, model parameters
- **`fvSolution`** - Linear solver settings, PIMPLE configuration
- **`fvSchemes`** - Discretization schemes (ddt, div, laplacian)
- **`blockMeshDict`** - Mesh generation
- **`decomposeParDict`**, **`setFieldsDict`**, **`topoSetDict`** - Domain setup

### constant/ Directory
- **`thermophysicalProperties`** - Global mixture data
- **`thermophysicalProperties.metal`** - Metal phase (titanium)
- **`thermophysicalProperties.air`** - Gas phase (air)
- **`laserProperties`** - Femtosecond laser parameters
- **`transportProperties`** - Mixture transport data
- **`turbulenceProperties`** - Turbulence model selection
- **`fvOptions`** - Optional source terms (gravity, etc.)

### 0.orig/ Directory
- **`p`**, **`p_rgh`** - Pressure fields
- **`U`** - Velocity field
- **`T`**, **`Te`**, **`Tl`** - Temperature fields (gas, electron, lattice)
- **`alpha.metal`**, **`alpha.air`** - Phase fractions

---

## Additional Documentation

### Models/ Directory
Comprehensive technical reports (markdown):
- `compInterFoam_Implementation_Report.md` - Solver architecture
- `FemtosecondLaserModel_Implementation_Report.md` - Laser physics
- `TwoTemperatureModel_Detailed_Report.md` - Electron-lattice coupling
- `AdvancedInterfaceCapturing_Thesis_Report.md` - Recoil pressure
- `twoPhaseMixtureThermo_Technical_Report.md` - Phase-change
- `compressibleInterPhaseTransportModel_Thesis_Report.md` - Transport closures

### Other Documents
- **`README.md`** - User guide and configuration parameter reference
- **`Thesis_MinhTriNguyen_10754416.pdf`** - Original research thesis (8.7 MB)
- **`EXPERIMENTAL_REPLICATE_ANALYSIS.md`** - Comparison with literature

---

## How to Use This Documentation

### I want to understand the overall architecture
1. Start with **ARCHITECTURE_DIAGRAM.txt** - See the execution flow
2. Read **Section 9** of **CODEBASE_MAP.md** - Key physics summary

### I need to find a specific file or function
1. Use **QUICK_REFERENCE.txt** - Lists all main files with locations
2. Cross-reference **CODEBASE_MAP.md** for detailed descriptions

### I want to modify simulation parameters
1. Check **QUICK_REFERENCE.txt section 6** - Configuration structure
2. Consult **CODEBASE_MAP.md section 6** - Complete parameter descriptions
3. See **ARCHITECTURE_DIAGRAM.txt** - Configuration hierarchy

### I want to understand a specific physics model (e.g., laser, phase-change, two-temperature)
1. Find the model in **QUICK_REFERENCE.txt section 3**
2. Go to **CODEBASE_MAP.md section 3** - Detailed physics description
3. Read the corresponding report in **Models/** directory for deep dives

### I need to add/modify equations
1. Look at the relevant equation file:
   - Momentum: **UEqn.H**
   - Pressure: **pEqn.H**
   - Phase-fraction: **alphaEqn.H**
   - Temperature: **TEqn.H**
2. Check **ARCHITECTURE_DIAGRAM.txt** - Understand data dependencies
3. Consult **CODEBASE_MAP.md section 7** - Coupling details

### I want to understand time-stepping and PIMPLE coupling
1. See **ARCHITECTURE_DIAGRAM.txt** - Complete time-stepping hierarchy
2. Check **QUICK_REFERENCE.txt section 7** - Key algorithms
3. Read **CODEBASE_MAP.md section 2** - Detailed equation descriptions

---

## Quick Navigation

### Key File Locations
| Purpose | File |
|---------|------|
| Main solver | `/home/user/compInterFoam/compInterFoam.C` |
| Phase-change | `/home/user/compInterFoam/twoPhaseMixtureThermo.C` |
| Electron-lattice | `/home/user/compInterFoam/twoTemperatureModel.C` |
| Laser heating | `/home/user/compInterFoam/femtosecondLaserModel.C` |
| Recoil pressure | `/home/user/compInterFoam/advancedInterfaceCapturing.C` |
| Simulation control | `/home/user/compInterFoam/TestCase/system/controlDict` |
| Discretization | `/home/user/compInterFoam/TestCase/system/fvSchemes` |
| Solvers | `/home/user/compInterFoam/TestCase/system/fvSolution` |
| Material properties | `/home/user/compInterFoam/TestCase/constant/thermophysicalProperties.metal` |
| Laser parameters | `/home/user/compInterFoam/TestCase/constant/laserProperties` |

### Key Equations Summary
- **Momentum:** `∂(ρU)/∂t + ∇·(ρU⊗U) + ∇p = ∇·τ + ∇(P_recoil)`
- **Phase-fraction:** `∂α/∂t + ∇·(αU) = -∇·(αU_interface) + source`
- **Electron temp:** `Ce·∂Te/∂t + ... = G(Tl-Te) + Q_laser + h_gm(Tgas-Te)`
- **Lattice temp:** `Cl·∂Tl/∂t + ... = G(Te-Tl) - Γ_evap·hf + h_gm(Tgas-Tl)`
- **Evaporation:** `Γ ~ √(p_sat/(2πMkT))` (Hertz-Knudsen)
- **Recoil:** `P_recoil ~ √(2πnMkT)·Γ` (kinetic theory)

### Critical Parameters for LIFT on Titanium
- **Laser:** 343 nm, 200 fs, ~20 nJ, 3.2 µm spot
- **Material:** Cp 650 J/(kg·K), ρ 4515 kg/m³, Tvap 3560 K, hf 9.1e6 J/kg
- **Time-step:** 0.1–5 fs (adaptive)
- **Phase-change:** Tsol 1941 K, alphaMin 0.001, alphaMax 0.999
- **Two-temp:** Ce~210, Cl~2.3e6, G~5e17, De~1e-4 (J/m³/K)
- **Recoil:** 50–100 MPa max, temperature-ramped

---

## Building the Solver

```bash
# Source OpenFOAM environment
source /opt/openfoam/etc/bashrc

# Navigate to solver directory
cd /home/user/compInterFoam

# Compile
wmake

# Result: executable at $(FOAM_USER_APPBIN)/compInterFoam
```

## Running a Simulation

```bash
# Set up case (copy from TestCase or TEST1/TEST2)
cd /path/to/MyCase

# Generate mesh
blockMesh

# Initialize fields
setFields

# Run solver
compInterFoam > log.out 2>&1 &

# Monitor progress
tail -f log.out

# Post-process with ParaView
paraFoam -builtin
```

---

## Document Statistics

- **Total lines:** 1,276 documentation lines (not counting source code)
- **Coverage:** All major components and configuration files documented
- **Generated:** 2025-11-16
- **Format:** Plain text + markdown for git-friendly version control

---

## Related Resources

1. **OpenFOAM Documentation:** https://www.openfoam.com/documentation
2. **VOF Method:** OpenFOAM Two-Phase Interface Capturing (MULES, CICSAM)
3. **PIMPLE Algorithm:** Semi-implicit Method for Pressure-Linked Equations
4. **Femtosecond LIFT:** See thesis "Laser-Induced Forward Transfer" in root directory
5. **Two-Temperature Model:** Electron-phonon coupling in femtosecond laser-material interaction

---

For questions about specific implementations, consult the corresponding `.md` report in the `Models/` directory.

**Document Generated:** 2025-11-16  
**Codebase Repository:** `/home/user/compInterFoam`
