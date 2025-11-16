# Implementation Guide for Enhanced Physics Models

This guide explains how to integrate the newly created physics models into the compInterFoam solver.

## New Physics Models Created

1. **Phase Explosion Model**
   - `phaseExplosionModel.H`
   - `phaseExplosionModel.C`

2. **Plasma Ionization Model**
   - `plasmaIonizationModel.H`
   - `plasmaIonizationModel.C`

3. **Droplet Breakup Model**
   - `dropletBreakupModel.H`
   - `dropletBreakupModel.C`

---

## Integration Steps

### Step 1: Update Make/files

Add the new source files to the compilation list:

**File:** `Make/files` (create if it doesn't exist)

```makefile
compInterFoam.C
twoPhaseMixtureThermo.C
twoTemperatureModel.C
femtosecondLaserModel.C
advancedInterfaceCapturing.C
phaseExplosionModel.C
plasmaIonizationModel.C
dropletBreakupModel.C

EXE = $(FOAM_USER_APPBIN)/compInterFoam
```

### Step 2: Update Make/options

Ensure proper include paths and libraries:

**File:** `Make/options` (create if it doesn't exist)

```makefile
EXE_INC = \
    -I$(LIB_SRC)/finiteVolume/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude \
    -I$(LIB_SRC)/sampling/lnInclude \
    -I$(LIB_SRC)/transportModels/twoPhaseMixture/lnInclude \
    -I$(LIB_SRC)/transportModels/interfaceProperties/lnInclude \
    -I$(LIB_SRC)/transportModels/incompressible/lnInclude \
    -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
    -I$(LIB_SRC)/thermophysicalModels/specie/lnInclude \
    -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
    -I$(LIB_SRC)/TurbulenceModels/incompressible/lnInclude \
    -I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
    -I$(LIB_SRC)/dynamicFvMesh/lnInclude

EXE_LIBS = \
    -lfiniteVolume \
    -lfvOptions \
    -lmeshTools \
    -lsampling \
    -ltwoPhaseMixture \
    -linterfaceProperties \
    -lincompressibleTransportModels \
    -lfluidThermophysicalModels \
    -lspecie \
    -lturbulenceModels \
    -lincompressibleTurbulenceModels \
    -lcompressibleTurbulenceModels \
    -ldynamicFvMesh
```

### Step 3: Integrate Models into Main Solver

**File:** `compInterFoam.C`

Add the following includes near the top of the file (after existing model includes):

```cpp
// Enhanced physics models
#include "phaseExplosionModel.H"
#include "plasmaIonizationModel.H"
#include "dropletBreakupModel.H"
```

**Initialize models after existing model initialization (around line 200-300):**

```cpp
// Initialize enhanced physics models
autoPtr<phaseExplosionModel> phaseExplosion;
autoPtr<plasmaIonizationModel> plasmaModel;
autoPtr<dropletBreakupModel> breakupModel;

const dictionary& controlDict = runTime.controlDict();

if (controlDict.found("phaseExplosionCoeffs"))
{
    phaseExplosion.reset
    (
        new phaseExplosionModel
        (
            mesh,
            controlDict.subDict("phaseExplosionCoeffs")
        )
    );
}

if (controlDict.found("plasmaIonizationCoeffs"))
{
    plasmaModel.reset
    (
        new plasmaIonizationModel
        (
            mesh,
            controlDict.subDict("plasmaIonizationCoeffs")
        )
    );
}

if (controlDict.found("dropletBreakupCoeffs"))
{
    breakupModel.reset
    (
        new dropletBreakupModel
        (
            mesh,
            mixture,
            controlDict.subDict("dropletBreakupCoeffs")
        )
    );
}
```

**In the time loop, after two-temperature model update and before PIMPLE loop:**

```cpp
// Update enhanced physics models
if (phaseExplosion.valid() && phaseExplosion->enabled())
{
    phaseExplosion->update(Tl, alpha1, rho);
}

if (plasmaModel.valid() && plasmaModel->enabled())
{
    plasmaModel->update(Tl, rho, alpha1);
}

if (breakupModel.valid() && breakupModel->enabled())
{
    breakupModel->update(alpha1, U, rho);
}
```

**After the PIMPLE loop, apply droplet breakup:**

```cpp
// Apply droplet breakup to alpha field
if (breakupModel.valid() && breakupModel->enabled())
{
    breakupModel->applyBreakup(alpha1, U);
}
```

**Integrate phase explosion and plasma pressures into momentum equation:**

In the section where recoil pressure is applied to the momentum equation, add:

```cpp
// Get recoil pressure from advanced interface capturing
volScalarField recoilP = advancedInterface->recoilPressure();

// Add phase explosion pressure
if (phaseExplosion.valid() && phaseExplosion->enabled())
{
    recoilP += phaseExplosion->explosivePressure();
}

// Add plasma pressure
if (plasmaModel.valid() && plasmaModel->enabled())
{
    recoilP += plasmaModel->plasmaPressure();
}

// Apply total pressure to momentum equation
fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  - fvm::laplacian(mixture.mu(), U)
  ==
  - fvc::grad(recoilP)
  + ... // other terms
);
```

**Enhance evaporation with phase explosion multiplier:**

In `twoPhaseMixtureThermo.C`, in the `computePhaseChange()` function, after computing `j_net`:

```cpp
// Check for phase explosion enhancement
if (mesh.foundObject<volScalarField>("explosionIndicator"))
{
    const volScalarField& explosionInd =
        mesh.lookupObject<volScalarField>("explosionIndicator");

    if (explosionInd[cellI] > 0.01)
    {
        // Enhance mass flux by explosion multiplier
        const scalar multiplier = 1.0 + 99.0 * explosionInd[cellI]; // 1-100x
        j_net *= multiplier;
    }
}
```

**Reduce laser absorption due to plasma shielding:**

In `femtosecondLaserModel.C`, in the laser absorption calculation:

```cpp
// Check for plasma shielding
scalar shieldingFactor = 1.0;
if (mesh.foundObject<volScalarField>("plasmaShielding"))
{
    const volScalarField& plasmaShield =
        mesh.lookupObject<volScalarField>("plasmaShielding");

    shieldingFactor = 1.0 - plasmaShield[cellI];
}

// Apply shielding to laser intensity
scalar I_absorbed = I_laser * absorptionCoeff * shieldingFactor;
```

---

## Step 4: Compilation

```bash
cd /home/user/compInterFoam
wmake
```

If compilation errors occur, check:
1. All header files are in the same directory
2. Make/files includes all .C files
3. Make/options has correct library links

---

## Step 5: Testing

Run the enhanced test case:

```bash
cd TestCase
./Allrun
```

Monitor the output for:
- "Phase explosion model enabled"
- "Plasma ionization model enabled"
- "Droplet breakup model enabled"

Check log files for:
- "Phase explosion active: Explosive cells = ..."
- "Plasma formation active: Plasma cells = ..."
- "Droplet breakup active: Breakup cells = ..."

---

## Step 6: Validation

**Check field outputs:**

```bash
ls -lh TestCase/0.5e-09/
```

Expected new fields:
- `explosiveMassSource`
- `explosivePressure`
- `explosionIndicator`
- `ionizationDegree`
- `n_electron`
- `plasmaPressure`
- `plasmaShielding`
- `WeberNumber`
- `breakupIndicator`
- `dropletDiameter`

**Post-process results:**

```bash
paraFoam -case TestCase
```

Visualize:
- Temperature contours (check for >8550 K for phase explosion)
- Pressure fields (recoil + explosive + plasma)
- Velocity magnitude (target: >30 m/s)
- Weber number (check for breakup regions)
- Ionization degree (plasma formation)

**Quantitative checks:**

```bash
postProcess -func 'mag(U)' -latestTime -case TestCase
```

Expected outcomes:
- Maximum velocity > 30 m/s
- Peak temperature ~6000-10000 K
- Recoil pressure ~100 MPa - 3 GPa
- Phase explosion above 8550 K
- Plasma formation above 30000 K
- Weber number > 10 in jet regions

---

## Troubleshooting

**Issue: Models not recognized**
- Check controlDict has the configuration sections
- Verify `enablePhaseExplosion`, `enablePlasmaModel`, `enableDropletBreakup` are set to `true`

**Issue: Compilation errors**
- Ensure OpenFOAM environment is loaded: `source $WM_PROJECT_DIR/etc/bashrc`
- Check C++11 compatibility: add `-std=c++11` to Make/options if needed

**Issue: Simulation crashes**
- Reduce time step: `maxDeltaT 1e-15;`
- Increase nOuterCorrectors: `nOuterCorrectors 30;`
- Enable pressure clamping: `pressureClamp true;`

**Issue: No phase explosion detected**
- Check temperature field reaches T_spinodal (8550 K)
- Increase laser fluence: higher energy or smaller spot size
- Verify `T_spinodal` setting in controlDict

**Issue: No plasma formation**
- Check temperature exceeds T_ionization (30000 K)
- Lower T_ionization threshold for testing: `T_ionization 20000;`
- Ensure metal vapor present (alpha1 > 0.1)

---

## Performance Notes

**Memory usage:**
- Each new model adds ~10 fields
- Expect +20-30% memory overhead

**Computational cost:**
- Phase explosion: +2-5% (only in explosive cells)
- Plasma model: +3-7% (Saha equation iterative)
- Droplet breakup: +1-3% (Weber number calculation)
- Total overhead: ~10-15%

**Scaling:**
- Models are cell-local (embarrassingly parallel)
- MPI scaling unchanged from base solver

---

## Next Steps

After successful integration and testing:

1. **Adaptive Mesh Refinement (AMR):**
   - Enable dynamicFvMesh
   - Refine based on temperature gradients, Weber number, explosion indicator

2. **Non-ideal Equation of State:**
   - Implement stiffened gas or virial EOS for high-T vapor
   - Replace ideal gas assumption

3. **Temperature-dependent Properties:**
   - Add polynomials for Cp(T), k(T), μ(T)
   - Update material property tables

4. **Advanced Breakup:**
   - Lagrangian droplet tracking
   - Sub-grid fragmentation model

---

## Support

For questions or contributions:
- Documentation: `PHYSICS_ENHANCEMENTS.md`
- Codebase map: `CODEBASE_MAP.md`
- Architecture: `ARCHITECTURE_DIAGRAM.txt`

---

**Last Updated:** 2025-01-16
**Solver Version:** compInterFoam + Enhanced LIFT Physics
