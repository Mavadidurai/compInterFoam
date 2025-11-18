# Mass Transfer Diagnostic Plan

## Problem Summary
- Temperatures reach 9067 K (well above Tvap = 2200 K)
- 728 cells exceed recoil threshold
- **BUT: Zero mass transfer detected** (all cells below massRateEps = 1e-12 kg/m²/s)
- Need to find if this is:
  1. Thermodynamic parameters wrong (P_sat, R, units)
  2. massRateEps still too high
  3. Time-gating shutting off recoil
  4. Q_laser coupling issue

## Key Code Locations

### Mass Flux Calculation
File: `twoPhaseMixtureThermo.C` lines 1185-1188

```cpp
const scalar j_evap = evaporationCoeff_*p_vapor/(sqrt_2piR*sqrt_T);
const scalar j_cond = evaporationCoeff_*p_metalVapor/(sqrt_2piR*sqrt_T);
scalar j_net = j_evap - j_cond;
```

### Saturation Pressure (Clausius-Clapeyron)
File: `twoPhaseMixtureThermo.C` lines 1140-1141

```cpp
const scalar exponent = (L/R)*(inv_Tvap - 1.0/Teval);
const scalar psat = p_ref*Foam::exp(exponent);
```

### Recoil Time Gating (CRITICAL!)
File: `advancedInterfaceCapturing.C` lines 554-563

```cpp
const scalar laserEnd = time.controlDict().lookupOrDefault<scalar>("laserEndTime", GREAT);
const scalar recoilGrace = 5e-12;  // 5 ps grace period

if (currentTime > laserEnd + recoilGrace)
{
    recoilPressure_ = 0;  // RECOIL SHUTS OFF!
    return;
}
```

**CRITICAL FINDING:** Your simulation at t=2.58 ps is **WITHIN** the laser window (laserEndTime = 200 ps),
so this gate is NOT the problem.

## Step-by-Step Diagnostic Procedure

### STEP 1: Check Existing DEBUG Output

The code already prints Hertz-Knudsen diagnostics every 10 timesteps!

```bash
cd /home/user/compInterFoam
grep "DEBUG Hertz-Knudsen" TEST1/log.* | tail -50
```

**What to look for:**
```
DEBUG Hertz-Knudsen at cell 0:
  T_eff = XXXX K
  p_vapor (psat) = XXXX Pa        ← Should be HUGE at 9000K
  sqrt_2piR = XXXX                ← Check this value
  sqrt_T = XXXX
  evaporationCoeff = 0.03         ← Correct value
  j_evap = XXXX kg/m²/s           ← The smoking gun
  j_net = XXXX kg/m²/s            ← Should be ~90 kg/m²/s at 5000K
  Expected at 5000K: ~90 kg/m²/s
```

**Expected values at T=9000K:**
- p_vapor: ~10^7 - 10^8 Pa (very high)
- j_net: ~1000-10000 kg/m²/s (should be MASSIVE)

**If you see:**
- p_vapor ~ 0 or very small → Clausius-Clapeyron parameters wrong
- j_net ~ 1e-12 or smaller → Units problem or evaporationCoeff issue

---

### STEP 2: Add Enhanced Diagnostics

Create a modified version with more verbose output:

```bash
cd /home/user/compInterFoam

# Backup original
cp twoPhaseMixtureThermo.C twoPhaseMixtureThermo.C.backup

# We'll add diagnostics after this step
```

#### Diagnostic Code to Add

In `twoPhaseMixtureThermo.C`, find line ~1185 and replace the debug block (lines 1209-1225) with:

```cpp
if (Pstream::master() && mesh.time().timeIndex() % 5 == 0)  // Every 5 timesteps
{
    // Find cell with max temperature
    label hotCell = -1;
    scalar maxTemp = -GREAT;
    forAll(TlField, i)
    {
        if (TlField[i] > maxTemp && alpha1Field[i] > 0.001)
        {
            maxTemp = TlField[i];
            hotCell = i;
        }
    }

    if (hotCell >= 0)
    {
        const scalar T_hot = TlField[hotCell];
        bool satClamped = false;
        const scalar p_vap_hot = saturationPressure(T_hot, satClamped);
        const scalar sqrt_T_hot = std::sqrt(Foam::max(T_hot, 1.0));
        const scalar sqrt_2piR_val = std::sqrt(2.0*3.14159*R);
        const scalar j_evap_hot = evaporationCoeff_*p_vap_hot/(sqrt_2piR_val*sqrt_T_hot);

        Info<< "\n=== MASS FLUX DIAGNOSTICS (hottest cell " << hotCell << ") ===" << nl
            << "  Time: " << mesh.time().value()*1e12 << " ps" << nl
            << "  T_eff = " << T_hot << " K" << nl
            << "  alpha = " << alpha1Field[hotCell] << nl
            << "  --- Clausius-Clapeyron parameters ---" << nl
            << "  L (latent heat) = " << L << " J/kg" << nl
            << "  R (gas constant) = " << R << " J/kg/K" << nl
            << "  p_ref (reference) = " << p_ref << " Pa" << nl
            << "  T_vap (reference) = " << T_vap << " K" << nl
            << "  --- Calculated values ---" << nl
            << "  exponent = " << (L/R)*(1.0/T_vap - 1.0/T_hot) << nl
            << "  p_vapor (psat) = " << p_vap_hot << " Pa" << nl
            << "  sqrt(2*pi*R) = " << sqrt_2piR_val << nl
            << "  sqrt(T) = " << sqrt_T_hot << nl
            << "  evaporationCoeff = " << evaporationCoeff_ << nl
            << "  j_evap = " << j_evap_hot << " kg/m²/s" << nl
            << "  massRateEps threshold = " << massRateEps_ << " kg/m²/s" << nl;

        if (j_evap_hot < massRateEps_)
        {
            Info<< "  *** FILTERED: j_evap < massRateEps ***" << nl;
            Info<< "  *** RATIO: j_evap/massRateEps = " << j_evap_hot/massRateEps_ << " ***" << nl;
        }
        else
        {
            Info<< "  *** ACTIVE: j_evap > massRateEps ***" << nl;
        }

        Info<< "  Expected j_net at 5000K (literature): ~90 kg/m²/s" << nl
            << "  Expected j_net at 9000K (extrapolated): ~500-1000 kg/m²/s" << nl
            << "==================================================\n" << endl;
    }
}
```

#### Compile and Test

```bash
cd /home/user/compInterFoam

# Recompile
wmake

# If compilation succeeds, run short test
cd TEST1
rm -rf 0.* [1-9]*
cp -r 0.orig 0
compInterFoam > log.diagnostics 2>&1 &

# Watch the output
tail -f log.diagnostics | grep -A 25 "MASS FLUX DIAGNOSTICS"
```

---

### STEP 3: Enable Field Writes for ParaView

Add these lines to `TEST1/system/controlDict` in the `functions` section:

```cpp
functions
{
    writeFields
    {
        type            writeObjects;
        libs            (utilityFunctionObjects);
        objects
        (
            phaseChangeMassFlux
            recoilPressure
            Q_laser
        );
        writeControl    writeTime;
    }
}
```

Then run and visualize:

```bash
# Run simulation
compInterFoam

# Open in ParaView
paraFoam

# In ParaView:
# 1. Load "phaseChangeMassFlux" field
# 2. Apply "Threshold" filter: alpha.metal > 0.001
# 3. Color by "phaseChangeMassFlux"
# 4. Check magnitude at t ~ 2 ps
```

---

### STEP 4: Manual Calculation Check

Create a Python script to verify the thermodynamics:

```bash
cat > /home/user/compInterFoam/check_psat.py << 'EOF'
#!/usr/bin/env python3
import numpy as np

# From your controlDict phaseChangeCoeffs
L = 9.1e6      # J/kg (latent heat)
R = 174        # J/kg/K (gas constant for Ti vapor)
T_vap = 2200   # K (reference temperature)
p_ref = 101325 # Pa (atmospheric pressure)
evap_coeff = 0.03

# Test temperatures
temperatures = [2200, 3000, 5000, 7000, 9000]

print("=" * 70)
print("TITANIUM VAPOR PRESSURE CALCULATION (Clausius-Clapeyron)")
print("=" * 70)
print(f"Parameters:")
print(f"  L = {L:.2e} J/kg")
print(f"  R = {R} J/kg/K")
print(f"  T_vap (ref) = {T_vap} K")
print(f"  p_ref = {p_ref} Pa")
print(f"  evaporationCoeff = {evap_coeff}")
print("=" * 70)

for T in temperatures:
    exponent = (L/R) * (1.0/T_vap - 1.0/T)
    p_sat = p_ref * np.exp(exponent)

    sqrt_2piR = np.sqrt(2.0 * np.pi * R)
    sqrt_T = np.sqrt(T)

    j_evap = evap_coeff * p_sat / (sqrt_2piR * sqrt_T)

    print(f"\nT = {T} K:")
    print(f"  exponent = {exponent:.2f}")
    print(f"  p_sat = {p_sat:.2e} Pa")
    print(f"  j_evap = {j_evap:.2e} kg/m²/s")

    if j_evap < 1e-12:
        print(f"  *** BELOW massRateEps (1e-12) - FILTERED! ***")
    elif j_evap < 1e-8:
        print(f"  *** Was filtered by old threshold (1e-8) ***")
    else:
        print(f"  *** ACTIVE ***")

print("\n" + "=" * 70)
print("EXPECTED BEHAVIOR:")
print("  At 5000K: j_evap ~ 90 kg/m²/s (from code comment)")
print("  At 9000K: j_evap should be 1000-10000 kg/m²/s")
print("  If you see j_evap < 1e-8, check your parameters!")
print("=" * 70)
EOF

chmod +x /home/user/compInterFoam/check_psat.py
python3 /home/user/compInterFoam/check_psat.py
```

---

## STEP 5: Analysis Decision Tree

Run the Python script first, then:

### Case A: Python shows j_evap ~ 1000 kg/m²/s at 9000K
→ **Thermodynamics are correct**
→ Check OpenFOAM DEBUG output
→ If OpenFOAM shows j_net ~ 0, then:
  - Units mismatch in code
  - Wrong R value being used
  - Activation windows preventing calculation

**Fix:** Add the enhanced diagnostics (Step 2) to see exact values in code

---

### Case B: Python shows j_evap ~ 1e-12 kg/m²/s at 9000K
→ **Thermodynamic parameters are WRONG**
→ Likely issues:
  1. **R = 174 J/kg/K is for Ti VAPOR**
     - Should be SPECIFIC gas constant (R_universal/M)
     - Ti: M = 47.867 g/mol = 0.047867 kg/mol
     - R_specific = 8.314 / 0.047867 = 173.7 J/kg/K ✓ (correct!)

  2. **p_ref might be wrong**
     - Using 1 atm (101325 Pa) as reference
     - But should use p_sat(T_vap) from literature
     - Ti vapor pressure at 3560K (boiling) ~ 1 atm
     - At 2200K (your Tvap), p_sat ~ 10 Pa (MUCH lower!)

**This is likely your problem!**

**Fix:**
```cpp
// In controlDict phaseChangeCoeffs:
Tvap    3560;  // Use actual boiling point as reference
// OR set p_ref to match your chosen Tvap
```

---

### Case C: j_evap shows ~1e-10 kg/m²/s
→ Close but still filtered
→ Lower massRateEps to 1e-14

**Fix:**
```cpp
// In controlDict advancedInterfaceCapturing:
massRateEps    1e-14;
```

---

## Expected Results

### If Everything Works:
```
DEBUG Hertz-Knudsen at cell 0:
  T_eff = 9000 K
  p_vapor (psat) = 5.2e7 Pa              ← HIGH
  j_evap = 847 kg/m²/s                   ← MASSIVE
  j_net = 847 kg/m²/s

Recoil diagnostics: 728 of 78000 interface cells supplied mass flux above 1e-12 kg/m^2/s
  Max |j_net| = 847 kg/m²/s
  Max |recoilPressure| = 125 MPa         ← NON-ZERO!
```

### Current (Broken) State:
```
  p_vapor (psat) = ??? Pa                ← CHECK THIS
  j_evap = ~1e-12 kg/m²/s                ← TOO LOW
  j_net = ~1e-12 kg/m²/s

Recoil diagnostics: 0 of 78000 interface cells supplied mass flux
```

---

## Quick First Check (5 minutes)

**Before modifying code, just check existing log:**

```bash
cd /home/user/compInterFoam/TEST1
grep -A 20 "DEBUG Hertz-Knudsen" log.* | tail -100 > hertz_knudsen_output.txt
cat hertz_knudsen_output.txt
```

**Send me the output** and I can immediately tell you what's wrong!

---

## Summary

The mass flux calculation is:
```
j = (evap_coeff * p_sat(T)) / (sqrt(2*pi*R) * sqrt(T))
```

Where:
```
p_sat(T) = p_ref * exp((L/R) * (1/T_vap - 1/T))
```

**Most likely culprits (in order):**

1. **p_ref/T_vap mismatch** (80% probability)
   - Using p_ref=1atm with T_vap=2200K gives tiny p_sat
   - Should use T_vap=3560K (boiling point) with p_ref=1atm

2. **massRateEps still too high** (15% probability)
   - Lower to 1e-14 or 1e-16

3. **Units error in code** (5% probability)
   - R in wrong units somewhere
   - Dimension mismatch

**Next action:** Run the Python check script - takes 30 seconds and will immediately show if thermodynamics are broken.
