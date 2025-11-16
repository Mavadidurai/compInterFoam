# Enhanced LIFT Physics - Combined Single File

**Simplified Integration:** All three physics models are now in a single file!

---

## Files

**Combined Implementation (Use These):**
- ✅ `enhancedLIFTPhysics.H` - Single header with all 3 models
- ✅ `enhancedLIFTPhysics.C` - Single implementation with all 3 models

**Individual Files (Legacy - Can be deleted):**
- ~~`phaseExplosionModel.H/C`~~ - Superseded by combined file
- ~~`plasmaIonizationModel.H/C`~~ - Superseded by combined file
- ~~`dropletBreakupModel.H/C`~~ - Superseded by combined file

---

## Quick Start

### 1. Add One Include

```cpp
#include "enhancedLIFTPhysics.H"
```

### 2. Initialize Once

```cpp
autoPtr<enhancedLIFTPhysics> liftPhysics
(
    new enhancedLIFTPhysics
    (
        mesh,
        mixture,
        runTime.controlDict()
    )
);
```

### 3. Update All Models (One Call)

```cpp
// Updates all three models in one call
liftPhysics->updateAll(T, alpha1, rho, U);
```

### 4. Access Results

```cpp
// Phase explosion
if (liftPhysics->phaseExplosionEnabled())
{
    const volScalarField& p_expl = liftPhysics->explosivePressure();
    totalPressure += p_expl;
}

// Plasma ionization
if (liftPhysics->plasmaEnabled())
{
    const volScalarField& p_plasma = liftPhysics->plasmaPressure();
    totalPressure += p_plasma;
}

// Droplet breakup
if (liftPhysics->breakupEnabled())
{
    liftPhysics->applyBreakup(alpha1, U);
}
```

---

## Configuration

**Single dictionary wraps all three models:**

```cpp
enhancedLIFTPhysics
{
    phaseExplosionCoeffs { ... }
    plasmaIonizationCoeffs { ... }
    dropletBreakupCoeffs { ... }
}
```

See `TestCase/system/controlDict` for complete example.

---

## Compilation

**Make/files:**
```makefile
compInterFoam.C
twoPhaseMixtureThermo.C
twoTemperatureModel.C
femtosecondLaserModel.C
advancedInterfaceCapturing.C
enhancedLIFTPhysics.C        # Add this line

EXE = $(FOAM_USER_APPBIN)/compInterFoam
```

Then:
```bash
wmake
```

---

## Advantages

✅ **Single file** - easier to manage
✅ **Unified API** - `updateAll()` updates all models
✅ **Cleaner includes** - only one header
✅ **Better organization** - all physics in one place
✅ **Same functionality** - all 3 models included

---

## Output Fields

When enabled, you'll get these fields automatically:

**Phase Explosion:**
- `explosivePressure`
- `explosiveMassSource`
- `explosionIndicator`

**Plasma Ionization:**
- `plasmaPressure`
- `ionizationDegree`
- `n_electron`
- `plasmaShielding`

**Droplet Breakup:**
- `WeberNumber`
- `breakupIndicator`
- `dropletDiameter`
- `breakupRate`

---

## Documentation

- **INTEGRATION_STEPS.md** - Step-by-step integration guide
- **PHYSICS_ENHANCEMENTS.md** - Detailed physics documentation
- **IMPLEMENTATION_GUIDE.md** - Full implementation details
- **SUMMARY.md** - Overview and status

---

## Support

For detailed integration instructions, see **INTEGRATION_STEPS.md**

For physics details and equations, see **PHYSICS_ENHANCEMENTS.md**

---

**Version:** 2.0 (Combined File)
**Date:** 2025-01-16
**Status:** Ready for integration
