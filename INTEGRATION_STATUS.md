# Integration Status - Enhanced LIFT Physics

## ✅ COMPLETED: Plasma Shielding Fix

**Problem:** Plasma shielding code was incorrectly placed at line 369 (outside function scope), causing:
```
femtosecondLaserModel.C:369:1: error: expected unqualified-id before 'if'
femtosecondLaserModel.C:376:21: error: 'I_laser' was not declared in this scope
```

**Solution:** Code has been moved to **correct location** at lines 1593-1607 inside the cell loop where:
- Variables `cellI`, `baseIntensity` are in scope
- Laser intensity is calculated and ready to be modified
- Plasma shielding reduces effective laser intensity before Beer-Lambert absorption

**Modified File:**
- `femtosecondLaserModel.C` (lines 1593-1607)

**Code Added:**
```cpp
// Check for plasma shielding (reduces laser absorption)
scalar shieldingFactor = 1.0;
if (mesh_.foundObject<volScalarField>("plasmaShielding"))
{
    const volScalarField& plasmaShield =
        mesh_.lookupObject<volScalarField>("plasmaShielding");

    shieldingFactor = Foam::max(1.0 - plasmaShield[cellI], scalar(0.0));
}

// Apply plasma shielding to laser intensity
const scalar effectiveIntensity = baseIntensity * shieldingFactor;
```

Then uses `effectiveIntensity` instead of `baseIntensity` for Ein/Eout calculations.

---

## 🔧 REMAINING INTEGRATION STEPS

The enhanced LIFT physics models are ready but **NOT YET INTEGRATED** into the main solver. Here's what still needs to be done:

### Step 1: Add to Make/files ⏳

**File:** `/home/user/compInterFoam/Make/files`

**Current status:** Missing `enhancedLIFTPhysics.C`

**Required change:**
```makefile
twoPhaseMixtureThermo.C
femtosecondLaserModel.C
twoTemperatureModel.C
advancedInterfaceCapturing.C
enhancedLIFTPhysics.C          # ← ADD THIS LINE
compressibleInterPhaseTransportModel.C
compInterFoam.C

EXE = $(FOAM_USER_APPBIN)/compInterFoam
```

---

### Step 2: Add Include Header ⏳

**File:** `/home/user/compInterFoam/compInterFoam.C`

**Location:** Around line 70 (after other LIFT model includes)

**Required change:**
```cpp
#include "femtosecondLaserModel.H"
#include "twoTemperatureModel.H"
#include "advancedInterfaceCapturing.H"
#include "enhancedLIFTPhysics.H"        // ← ADD THIS LINE
#include "OFstream.H"
```

---

### Step 3: Initialize Enhanced LIFT Physics ⏳

**File:** `/home/user/compInterFoam/createFields.H`

**Location:** Around line 733 (after advancedInterfaceCapturing initialization)

**Required change:**
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
        Info<< "Creating enhanced LIFT physics models\\n" << endl;
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

### Step 4: Update Models in Time Loop ⏳

**File:** `/home/user/compInterFoam/compInterFoam.C`

**Location:** Around line 1027 (after `ttm.solve(...)` and `mixture.setClTTM(...)`)

**Required change:**
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

### Step 5: Add Pressures to Momentum Equation ⏳

**File:** `/home/user/compInterFoam/UEqn.H`

**Location:** Where recoil pressure is applied to momentum equation

**Required change:**
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

### Step 6: Apply Droplet Breakup ⏳

**File:** `/home/user/compInterFoam/compInterFoam.C`

**Location:** Around line 1066 (after PIMPLE loop closes)

**Required change:**
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

### Step 7: Enhance Evaporation with Phase Explosion ⏳

**File:** `/home/user/compInterFoam/twoPhaseMixtureThermo.C`

**Location:** Around line 1188 (in `computePhaseChange` function, after `j_net = j_evap - j_cond`)

**Required change:**
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

## 📋 INTEGRATION CHECKLIST

Before compiling:
- [x] Plasma shielding added to femtosecondLaserModel.C ✅
- [ ] `enhancedLIFTPhysics.C` added to Make/files
- [ ] `#include "enhancedLIFTPhysics.H"` added to compInterFoam.C
- [ ] `liftPhysics` initialized in createFields.H
- [ ] `liftPhysics->updateAll()` added in time loop
- [ ] Enhanced pressures added to UEqn.H
- [ ] Phase explosion multiplier added to evaporation
- [ ] Droplet breakup applied after PIMPLE loop

After compiling:
- [ ] No compilation errors
- [ ] Solver executable created
- [ ] Test case runs without errors
- [ ] Startup messages show models enabled
- [ ] New fields appear in output (explosivePressure, plasmaPressure, WeberNumber)
- [ ] Enhanced physics reports activity during simulation

---

## 🎯 NEXT STEPS

1. **Complete remaining integrations** (Steps 1-7 above)
2. **Compile the solver:**
   ```bash
   cd /home/user/compInterFoam
   source $WM_PROJECT_DIR/etc/bashrc  # or wherever OpenFOAM is installed
   wmake
   ```
3. **Test the solver:**
   ```bash
   cd TestCase
   ./Allrun
   ```
4. **Verify enhanced physics:**
   - Check startup messages for "Enhanced LIFT physics initialized"
   - Monitor for "Phase explosion active", "Plasma ionization active", etc.
   - Check for new output fields in time directories

---

## 📚 REFERENCE DOCUMENTS

- **EXACT_INTEGRATION_LOCATIONS.md** - Detailed step-by-step with exact line numbers
- **INTEGRATION_FLOW_DIAGRAM.txt** - Visual flow diagram showing execution order
- **INTEGRATION_STEPS.md** - Quick integration guide
- **README_COMBINED_PHYSICS.md** - Overview of combined physics model

---

**Status:** Plasma shielding fix complete ✅ | Main integration pending ⏳

**Modified Files:** 1 (femtosecondLaserModel.C)

**Remaining Files to Modify:** 5 (Make/files, compInterFoam.C, createFields.H, UEqn.H, twoPhaseMixtureThermo.C)
