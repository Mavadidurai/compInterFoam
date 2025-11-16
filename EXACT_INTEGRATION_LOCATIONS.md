# Complete Integration Guide for enhancedLIFTPhysics.C

This guide shows you **exactly** where to add code in `compInterFoam.C` to integrate the combined physics models.

---

## Step 1: Add Include Header

**Location:** After line 70 (after the other LIFT model includes)

**Find this section:**
```cpp
// LIFT-specific model headers
#include "femtosecondLaserModel.H"
#include "twoTemperatureModel.H"
#include "advancedInterfaceCapturing.H"
#include "OFstream.H"
```

**Add this line:**
```cpp
// LIFT-specific model headers
#include "femtosecondLaserModel.H"
#include "twoTemperatureModel.H"
#include "advancedInterfaceCapturing.H"
#include "enhancedLIFTPhysics.H"        // ← ADD THIS LINE
#include "OFstream.H"
```

---

## Step 2: Initialize the Model

**Location:** In `createFields.H` after line 732 (right after advancedInterfaceCapturing initialization)

**Find this section (around line 720-732 in createFields.H):**
```cpp
    pInterfaceCapturing.reset
    (
        new advancedInterfaceCapturing
        (
            mesh,
            alpha1,
            mixture,
            T,
            &ttm.Tl()
        )
    );
    recoilPressurePtr = &pInterfaceCapturing->recoilPressure();
}
```

**Add this code RIGHT AFTER the closing brace:**
```cpp
    recoilPressurePtr = &pInterfaceCapturing->recoilPressure();
}

// ========== ADD THIS SECTION ==========
// Initialize enhanced LIFT physics (phase explosion, plasma, droplet breakup)
autoPtr<enhancedLIFTPhysics> liftPhysics;

if (runTime.controlDict().found("enhancedLIFTPhysics"))
{
    if (verbose && Foam::Pstream::master())
    {
        Info<< "Creating enhanced LIFT physics models\n" << endl;
    }

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
// ========== END OF ADDITION ==========
```

---

## Step 3: Update Models in Time Loop

**Location:** In `compInterFoam.C` around line 1027 (right after `ttm.solve(...)` and `mixture.setClTTM(...)`)

**Find this section (around line 1018-1030):**
```cpp
        ttm.solve
        (
            QLaser,
            phaseChangeSource,
            phaseChangeRelaxCoeff,
            gasMetalHeatFlux
        );

        mixture.setClTTM(ttm.Cl());
    }

    // ✅ UPDATE PHASE CHANGE WITH LATEST TEMPERATURES
    mixture.correct();
```

**Add this code BETWEEN `mixture.setClTTM(...)` and `mixture.correct()`:**
```cpp
        mixture.setClTTM(ttm.Cl());
    }

    // ========== ADD THIS SECTION ==========
    // Update enhanced LIFT physics (all three models updated in one call)
    if (liftPhysics.valid())
    {
        liftPhysics->updateAll(ttm.Tl(), alpha1, rho, U);

        if (verbose && master && runTime.timeIndex() % 100 == 0)
        {
            Info<< "Enhanced physics updated:" << nl;
            if (liftPhysics->phaseExplosionEnabled())
                Info<< "    Phase explosion active" << nl;
            if (liftPhysics->plasmaEnabled())
                Info<< "    Plasma ionization active" << nl;
            if (liftPhysics->breakupEnabled())
                Info<< "    Droplet breakup active" << nl;
        }
    }
    // ========== END OF ADDITION ==========

    // ✅ UPDATE PHASE CHANGE WITH LATEST TEMPERATURES
    mixture.correct();
```

---

## Step 4: Add Pressures to Momentum Equation

**Location:** In `UEqn.H`, find where recoil pressure is applied

**You need to edit:** `/home/user/compInterFoam/UEqn.H`

**Find the section that looks like:**
```cpp
// Apply recoil pressure
fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  - fvm::laplacian(mixture.mu(), U)
  ==
  - fvc::grad(p_rgh)
  - ghf*fvc::snGrad(rho)*mesh.magSf()
  + ... other terms
);
```

**Change it to:**
```cpp
// Calculate total pressure including enhanced physics
volScalarField totalPressure(p_rgh);

// Add recoil pressure from advanced interface capturing
if (useAdvancedCapturing && pInterfaceCapturing.valid())
{
    totalPressure += pInterfaceCapturing->recoilPressure();
}

// ========== ADD THESE CONTRIBUTIONS ==========
// Add phase explosion pressure
if (liftPhysics.valid() && liftPhysics->phaseExplosionEnabled())
{
    totalPressure += liftPhysics->explosivePressure();
}

// Add plasma pressure
if (liftPhysics.valid() && liftPhysics->plasmaEnabled())
{
    totalPressure += liftPhysics->plasmaPressure();
}
// ========== END OF ADDITIONS ==========

// Apply total pressure to momentum equation
fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  - fvm::laplacian(mixture.mu(), U)
  ==
  - fvc::grad(totalPressure)    // ← Use totalPressure instead of p_rgh
  - ghf*fvc::snGrad(rho)*mesh.magSf()
  + ... other terms
);
```

---

## Step 5: Apply Droplet Breakup

**Location:** In `compInterFoam.C` around line 1066 (right after the PIMPLE loop closes)

**Find this section (around line 1062-1068):**
```cpp
    if (pimple.turbCorr())
    {
        transportModel.correct();
    }
}

        dimensionedScalar Ek = fvc::domainIntegrate(0.5*rho*magSqr(U));
```

**Add this code RIGHT AFTER the closing brace of the PIMPLE loop:**
```cpp
    if (pimple.turbCorr())
    {
        transportModel.correct();
    }
}

// ========== ADD THIS SECTION ==========
// Apply droplet breakup to alpha field (interface remapping)
if (liftPhysics.valid() && liftPhysics->breakupEnabled())
{
    liftPhysics->applyBreakup(alpha1, U);
}
// ========== END OF ADDITION ==========

        dimensionedScalar Ek = fvc::domainIntegrate(0.5*rho*magSqr(U));
```

---

## Step 6: Enhance Evaporation with Phase Explosion

**Location:** In `twoPhaseMixtureThermo.C` around line 1188 (in the `computePhaseChange` function)

**Find this section (around line 1180-1200):**
```cpp
        const scalar j_evap = evaporationCoeff_*p_vapor/(sqrt_2piR*sqrt_T);
        const scalar j_cond = evaporationCoeff_*p_metalVapor/(sqrt_2piR*sqrt_T);

        scalar j_net = j_evap - j_cond;
```

**Add this code RIGHT AFTER `scalar j_net = j_evap - j_cond;`:**
```cpp
        scalar j_net = j_evap - j_cond;

        // ========== ADD THIS SECTION ==========
        // Apply phase explosion enhancement if active
        if (mesh.foundObject<volScalarField>("explosionIndicator"))
        {
            const volScalarField& explosionInd =
                mesh.lookupObject<volScalarField>("explosionIndicator");

            if (explosionInd[cellI] > 0.01)
            {
                // Enhanced mass flux: 1 to 100x based on superheat indicator
                const scalar multiplier = 1.0 + 99.0 * explosionInd[cellI];
                j_net *= multiplier;

                if (Foam::Pstream::master() && cellI == 0 && mesh.time().timeIndex() % 100 == 0)
                {
                    Info<< "  Phase explosion enhancement at cell 0: "
                        << multiplier << "x (indicator=" << explosionInd[cellI] << ")" << nl;
                }
            }
        }
        // ========== END OF ADDITION ==========
```

---

## Step 7: Add Plasma Shielding to Laser

**Location:** In `femtosecondLaserModel.C` in the laser absorption calculation

**Find the absorption calculation (should be in the `correct()` or `update()` function):**
```cpp
// Calculate absorbed laser intensity
scalar I_absorbed = I_laser * absorptionCoeff;
```

**Change it to:**
```cpp
// ========== MODIFY THIS SECTION ==========
// Check for plasma shielding
scalar shieldingFactor = 1.0;
if (mesh.foundObject<volScalarField>("plasmaShielding"))
{
    const volScalarField& plasmaShield =
        mesh.lookupObject<volScalarField>("plasmaShielding");

    shieldingFactor = max(1.0 - plasmaShield[cellI], 0.0);  // 0 to 1

    if (shieldingFactor < 0.99 && mesh.time().timeIndex() % 100 == 0)
    {
        Info<< "  Plasma shielding active: reducing laser by "
            << (1.0 - shieldingFactor)*100.0 << "%" << nl;
    }
}

// Calculate absorbed laser intensity with plasma shielding
scalar I_absorbed = I_laser * absorptionCoeff * shieldingFactor;
// ========== END OF MODIFICATION ==========
```

---

## Step 8: Update Make/files

**File:** `/home/user/compInterFoam/Make/files`

**Replace or create with:**
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

## Step 9: Create or Update Make/options

**File:** `/home/user/compInterFoam/Make/options`

**Create or replace with:**
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

---

## Step 10: Compile

```bash
cd /home/user/compInterFoam
wmake
```

**Expected output:**
```
Making dependency list for source file compInterFoam.C
...
Making dependency list for source file enhancedLIFTPhysics.C
...
g++ ... -o compInterFoam
```

---

## Step 11: Test

```bash
cd TestCase
./Allrun
```

**Watch for startup messages:**
```
Creating enhanced LIFT physics models

Enhanced LIFT physics initialized:
    Phase explosion: ON
    Plasma ionization: ON
    Droplet breakup: ON

Phase explosion model enabled:
    T_critical   = 9500 K
    T_spinodal   = 8550 K
    p_critical   = 120 MPa
    ...

Plasma ionization model enabled:
    Ionization energy = 6.82 eV
    T_ionization = 30000 K
    ...

Droplet breakup model enabled:
    We_critical = 10
    ...
```

**During simulation:**
```
Enhanced physics updated:
    Phase explosion active
    Plasma ionization active
    Droplet breakup active

Phase explosion active:
    Explosive cells = 1247
    Max T = 8923 K
    Max explosion pressure = 45.3 MPa
    ...
```

---

## Summary of Changes

| File | Changes | Lines Added |
|------|---------|-------------|
| `compInterFoam.C` | Add include, update time loop | ~30 |
| `createFields.H` | Initialize model | ~20 |
| `UEqn.H` | Add pressures to momentum | ~15 |
| `twoPhaseMixtureThermo.C` | Enhance evaporation | ~15 |
| `femtosecondLaserModel.C` | Add plasma shielding | ~15 |
| `Make/files` | Add source file | 1 |
| `Make/options` | Already exists | 0 |
| **Total** | | **~96 lines** |

---

## Quick Checklist

Before compiling:
- [ ] Added `#include "enhancedLIFTPhysics.H"` to compInterFoam.C
- [ ] Initialized `liftPhysics` in createFields.H
- [ ] Added `liftPhysics->updateAll()` in time loop
- [ ] Modified UEqn.H to include enhanced pressures
- [ ] Added phase explosion multiplier to evaporation
- [ ] Added plasma shielding to laser absorption
- [ ] Updated Make/files
- [ ] Verified Make/options exists

After compiling:
- [ ] No compilation errors
- [ ] Solver executable created
- [ ] Test case runs without errors
- [ ] Startup messages show models enabled
- [ ] New fields appear in output (explosivePressure, etc.)
- [ ] Enhanced physics reports activity during simulation

---

## Troubleshooting

**"undefined reference to enhancedLIFTPhysics"**
→ Make sure `enhancedLIFTPhysics.C` is in Make/files

**"No such file: enhancedLIFTPhysics.H"**
→ Check files are in /home/user/compInterFoam/ directory

**Models don't initialize**
→ Check controlDict has `enhancedLIFTPhysics { ... }` dictionary

**No active cells reported**
→ Temperatures too low; increase laser fluence

---

## File Locations

All modifications are in:
- `/home/user/compInterFoam/compInterFoam.C`
- `/home/user/compInterFoam/createFields.H`
- `/home/user/compInterFoam/UEqn.H`
- `/home/user/compInterFoam/twoPhaseMixtureThermo.C`
- `/home/user/compInterFoam/femtosecondLaserModel.C`
- `/home/user/compInterFoam/Make/files`

Configuration:
- `/home/user/compInterFoam/TestCase/system/controlDict` (already updated)

---

**Ready to integrate!** Follow the steps in order and you'll have all three physics models working together.
