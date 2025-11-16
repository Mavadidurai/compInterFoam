# Quick Integration Steps for Enhanced LIFT Physics

**Single Combined File:** All three physics models (phase explosion, plasma ionization, droplet breakup) are now in one file for easier integration.

---

## Step 1: Add to Make/files

Create or edit `/home/user/compInterFoam/Make/files`:

```makefile
compInterFoam.C
twoPhaseMixtureThermo.C
twoTemperatureModel.C
femtosecondLaserModel.C
advancedInterfaceCapturing.C
enhancedLIFTPhysics.C

EXE = $(FOAM_USER_APPBIN)/compInterFoam
```

---

## Step 2: Update Make/options

Create or edit `/home/user/compInterFoam/Make/options`:

```makefile
EXE_INC = \
    -I$(LIB_SRC)/finiteVolume/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude \
    -I$(LIB_SRC)/transportModels/twoPhaseMixture/lnInclude \
    -I$(LIB_SRC)/transportModels/interfaceProperties/lnInclude \
    -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
    -I$(LIB_SRC)/thermophysicalModels/specie/lnInclude \
    -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
    -I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
    -I$(LIB_SRC)/dynamicFvMesh/lnInclude

EXE_LIBS = \
    -lfiniteVolume \
    -lfvOptions \
    -lmeshTools \
    -ltwoPhaseMixture \
    -linterfaceProperties \
    -lfluidThermophysicalModels \
    -lspecie \
    -lturbulenceModels \
    -lcompressibleTurbulenceModels \
    -ldynamicFvMesh
```

---

## Step 3: Integrate into compInterFoam.C

### 3.1 Add Include (after existing includes):

```cpp
// Enhanced LIFT physics models
#include "enhancedLIFTPhysics.H"
```

### 3.2 Initialize Model (after mixture initialization):

```cpp
// Initialize enhanced LIFT physics
autoPtr<enhancedLIFTPhysics> liftPhysics;

if (runTime.controlDict().found("enhancedLIFTPhysics"))
{
    liftPhysics.reset
    (
        new enhancedLIFTPhysics
        (
            mesh,
            mixture,
            runTime.controlDict()
        )
    );
}
```

### 3.3 Update in Time Loop (after two-temperature model update):

```cpp
// Update enhanced LIFT physics
if (liftPhysics.valid())
{
    liftPhysics->updateAll(Tl, alpha1, rho, U);
}
```

### 3.4 Add Pressures to Momentum Equation:

Find where recoil pressure is applied and modify:

```cpp
// Get total pressure contributions
volScalarField totalRecoilP = advancedInterface->recoilPressure();

// Add phase explosion pressure
if (liftPhysics.valid() && liftPhysics->phaseExplosionEnabled())
{
    totalRecoilP += liftPhysics->explosivePressure();
}

// Add plasma pressure
if (liftPhysics.valid() && liftPhysics->plasmaEnabled())
{
    totalRecoilP += liftPhysics->plasmaPressure();
}

// Apply to momentum equation
fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  - fvm::laplacian(mixture.mu(), U)
  ==
  - fvc::grad(totalRecoilP)
  + ... // other terms
);
```

### 3.5 Apply Droplet Breakup (after PIMPLE loop):

```cpp
// Apply droplet breakup
if (liftPhysics.valid() && liftPhysics->breakupEnabled())
{
    liftPhysics->applyBreakup(alpha1, U);
}
```

### 3.6 Enhance Evaporation (in twoPhaseMixtureThermo::computePhaseChange()):

After computing j_net in the evaporation calculation:

```cpp
// Apply phase explosion enhancement if active
if (mesh.foundObject<volScalarField>("explosionIndicator"))
{
    const volScalarField& explosionInd =
        mesh.lookupObject<volScalarField>("explosionIndicator");

    if (explosionInd[cellI] > 0.01)
    {
        // Get multiplier (1 to 100x based on superheat)
        scalar T_local = TlField[cellI];

        // Simple enhancement: 1 + 99*indicator
        scalar multiplier = 1.0 + 99.0 * explosionInd[cellI];
        j_net *= multiplier;
    }
}
```

### 3.7 Add Plasma Shielding (in femtosecondLaserModel.C):

In the laser absorption calculation:

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

## Step 4: Configuration (controlDict)

Your controlDict already has the configuration blocks. Make sure they're under a parent dictionary:

```cpp
// Wrap all three model configs in one dictionary
enhancedLIFTPhysics
{
    phaseExplosionCoeffs
    {
        enablePhaseExplosion   true;
        T_critical             9500;
        T_spinodal             8550;
        p_critical             1.2e8;
        tau_explosion          1e-12;
        explosionMultiplier    100.0;
    }

    plasmaIonizationCoeffs
    {
        enablePlasmaModel      true;
        ionizationEnergy       1.093e-18;
        T_ionization           30000;
        atomicMass             7.95e-26;
        atomicNumberDensity    5.68e28;
    }

    dropletBreakupCoeffs
    {
        enableDropletBreakup   true;
        We_critical            10.0;
        minJetDiameter         1e-7;
        maxJetDiameter         10e-6;
        breakupTimeCoeff       1.0;
    }
}
```

---

## Step 5: Compile

```bash
cd /home/user/compInterFoam
wmake
```

---

## Step 6: Test

```bash
cd TestCase
./Allrun
```

Watch for startup messages:
```
Enhanced LIFT physics initialized:
    Phase explosion: ON
    Plasma ionization: ON
    Droplet breakup: ON

Phase explosion model enabled:
    T_critical   = 9500 K
    T_spinodal   = 8550 K
    ...

Plasma ionization model enabled:
    Ionization energy = 6.82 eV
    ...

Droplet breakup model enabled:
    We_critical = 10
    ...
```

During simulation:
```
Phase explosion active:
    Explosive cells = 1247
    Max T = 8923 K
    ...

Plasma formation active:
    Plasma cells = 347
    ...

Droplet breakup active:
    Breakup cells = 523
    ...
```

---

## Step 7: Post-Process

Check for new fields:
```bash
ls TestCase/0.5e-09/
```

Should see:
- `explosivePressure`
- `explosiveMassSource`
- `explosionIndicator`
- `plasmaPressure`
- `ionizationDegree`
- `n_electron`
- `plasmaShielding`
- `WeberNumber`
- `breakupIndicator`
- `dropletDiameter`

Visualize:
```bash
paraFoam -case TestCase
```

---

## Simplified API

With the combined file, you only need to:

```cpp
// 1. Include
#include "enhancedLIFTPhysics.H"

// 2. Initialize
autoPtr<enhancedLIFTPhysics> liftPhysics(new enhancedLIFTPhysics(mesh, mixture, dict));

// 3. Update (one call updates all three models)
liftPhysics->updateAll(T, alpha1, rho, U);

// 4. Access results
if (liftPhysics->phaseExplosionEnabled())
    const volScalarField& p_expl = liftPhysics->explosivePressure();

if (liftPhysics->plasmaEnabled())
    const volScalarField& p_plasma = liftPhysics->plasmaPressure();

if (liftPhysics->breakupEnabled())
    liftPhysics->applyBreakup(alpha1, U);
```

Much simpler than managing three separate models!

---

## Troubleshooting

**Compilation errors:**
- Make sure `enhancedLIFTPhysics.H` and `enhancedLIFTPhysics.C` are in solver directory
- Check Make/files includes `enhancedLIFTPhysics.C`
- Verify OpenFOAM environment: `source $WM_PROJECT_DIR/etc/bashrc`

**Models not initializing:**
- Check controlDict has `enhancedLIFTPhysics { ... }` wrapper dictionary
- Verify `enable*` switches are `true`

**No active cells:**
- Phase explosion requires T > 8550 K (increase laser fluence)
- Plasma requires T > 30000 K (very high, may not occur)
- Breakup requires We > 10 (high velocity jets)

---

## Summary

**Files Created:**
1. `enhancedLIFTPhysics.H` - Combined header (all 3 models)
2. `enhancedLIFTPhysics.C` - Combined implementation (all 3 models)

**Advantages:**
- ✓ Single include statement
- ✓ One initialization call
- ✓ Unified update: `updateAll()`
- ✓ Easier to manage and maintain
- ✓ Reduced file count: 6 → 2 files

**Integration Points:**
1. Make/files - add `enhancedLIFTPhysics.C`
2. compInterFoam.C - include, initialize, update
3. Momentum equation - add pressures
4. Evaporation - enhance with explosion multiplier
5. Laser absorption - reduce with plasma shielding
6. Alpha field - apply breakup

Follow this guide and you'll have all three physics models integrated and working!
