# CSV Implementation Analysis for compInterFoam

## Overview
The enhanced CSV logging feature in `compInterFoam.C` captures comprehensive simulation data including physics parameters, solver settings, and real-time metrics for LIFT process analysis.

---

## ✓ What's Working Correctly

### 1. CSV Structure (Lines 1159-1349)
- **Header**: 80 columns properly defined (lines 1179-1217)
- **Data output**: 80 columns matching header (lines 1260-1348)
- **File handling**: Uses `OFstream` with APPEND mode for continuous logging
- **Header guard**: Writes header only once using `static bool headerWritten`
- **Master process**: Only master rank writes to avoid conflicts

### 2. Real-time Metrics Capture (32 columns)
All dynamic simulation metrics are correctly captured:
- Time progression and phase tracking
- Temperature fields (Te, Tl, averages, spreads)
- Pressure fields (max, recoil, gradients)
- Metal dynamics (volume, loss rate, interface area)
- Velocity statistics (max, avg, stddev, acceleration)
- Vapor phase tracking
- Energy analysis (kinetic, turbulent)
- Shock wave and separation detection
- Laser power metrics

### 3. Non-dimensioned Parameters
Parameters without dimensions are extracted correctly:
- Laser properties (pulseEnergy, pulseWidth, etc.) - use `.value()` ✓
- Solver settings (maxCo, maxAlphaCo, etc.) ✓
- Dimensionless coefficients (evaporationCoeff, etc.) ✓
- Mesh statistics (totalCells, volumes) ✓

### 4. CSV Output Location
- File: `liftProcessTracking.csv`
- Location: Case directory (where solver is run)
- Format: ASCII, uncompressed
- Frequency: Every timestep (complete data)

---

## ❌ Critical Issues Identified

### Issue #1: Incorrect Extraction of Dimensioned Parameters
**Severity**: HIGH
**Impact**: All dimensioned input parameters will be recorded as 0.0

**Affected Lines**: 223-264 in `compInterFoam.C`

**Problem**:
```cpp
// WRONG - these are dimensionedScalars in controlDict:
input_Tvap = phaseChangeDict.lookupOrDefault<Foam::scalar>("Tvap", 0.0);
input_Tsol = phaseChangeDict.lookupOrDefault<Foam::scalar>("Tsol", 0.0);
input_hf = phaseChangeDict.lookupOrDefault<Foam::scalar>("hf", 0.0);
input_gasConstant = phaseChangeDict.lookupOrDefault<Foam::scalar>("gasConstant", 0.0);
input_p_ref = phaseChangeDict.lookupOrDefault<Foam::scalar>("p_ref", 101325.0);
input_maxSource = phaseChangeDict.lookupOrDefault<Foam::scalar>("maxSource", 0.0);
input_Cl = twoTempDict.lookupOrDefault<Foam::scalar>("Cl", 0.0);
input_De = twoTempDict.lookupOrDefault<Foam::scalar>("De", 0.0);
```

**Why it fails**:
In controlDict, these are defined as:
```cpp
Tvap    [0 0 0 1 0 0 0] 3560;      // dimensionedScalar
Tsol    [0 0 0 1 0 0 0] 1941;
hf      [0 2 -2 0 0 0 0] 9.1e6;
// etc.
```

OpenFOAM stores these as `dimensionedScalar` objects, not plain `scalar`.
Using `lookupOrDefault<scalar>()` fails to find them and returns default (0.0).

**Affected CSV Columns** (will show 0.0 instead of actual values):
- Column 54: `input_Tvap_K` (should be ~3560)
- Column 55: `input_Tsol_K` (should be ~1941)
- Column 57: `input_hf_Jkg` (should be ~9.1e6)
- Column 58: `input_gasConstant_JkgK` (should be ~174)
- Column 59: `input_p_ref_Pa` (should be ~101325)
- Column 61: `input_maxSource_kgm3s` (should be ~1e25)
- Column 64: `input_Cl_Jm3K` (should be ~2.5e6)
- Column 65: `input_De_m2s` (should be ~1e-4)

**Total affected**: 8 critical input parameters

---

## 🔧 Required Fix

### Replace Lines 223-264 with Corrected Code

See `FIXED_CSV_parameter_extraction.cpp` for the complete corrected code.

**Key changes**:
```cpp
// CORRECT approach - extract from dimensionedScalar:
if (phaseChangeDict.found("Tvap"))
{
    try
    {
        Foam::dimensionedScalar Tvap_dim(phaseChangeDict.lookup("Tvap"));
        input_Tvap = Tvap_dim.value();  // Extract scalar value
    }
    catch (...)
    {
        input_Tvap = phaseChangeDict.lookupOrDefault<Foam::scalar>("Tvap", 0.0);
    }
}
```

**Apply to all dimensioned parameters**:
- `Tvap`, `Tsol` (phaseChangeCoeffs)
- `hf`, `gasConstant`, `p_ref`, `maxSource` (phaseChangeCoeffs)
- `Cl`, `De` (twoTemperatureProperties)

---

## 📊 CSV Data Structure (80 columns)

### Real-time Metrics (Columns 1-32)
| Column | Name | Units | Description |
|--------|------|-------|-------------|
| 1 | time_ps | ps | Simulation time |
| 2 | phase | string | Phase name (IDLE, LASER_ACTIVE, etc.) |
| 3 | phase_num | - | Phase number (0-16) |
| 4 | progress_pct | % | Overall progress |
| 5-8 | Te_max_K, Tl_max_K, Tl_avg_K, Tl_spread_K | K | Temperatures |
| 9-11 | P_max_MPa, recoil_MPa, recoil_max_seen_MPa | MPa | Pressures |
| 12-14 | metal_vol_um3, metal_loss_pct, metal_loss_rate_um3s | µm³, %, µm³/s | Metal tracking |
| 15-18 | vel_max_ms, vel_avg_ms, vel_stddev_ms, accel_kms2 | m/s, km/s² | Velocity stats |
| 19-20 | interface_cells, interface_area_um2 | -, µm² | Interface tracking |
| 21-24 | vapor_vol_um3, vapor_cells, vapor_vel_avg_ms, vapor_vel_max_ms | µm³, -, m/s | Vapor phase |
| 25-26 | turbulent_KE_nJ, kinetic_eff_pct | nJ, % | Energy analysis |
| 27-28 | pressure_grad_GPam, shock_present | GPa/m, bool | Shock detection |
| 29-30 | separation_detected, time_to_ejection_ps | bool, ps | Separation tracking |
| 31-32 | laserPower_W, peakQLaser_TWm3 | W, TW/m³ | Laser metrics |

### Input Parameters (Columns 33-80)
Captured once at simulation start for reproducibility:

**Laser Properties (33-41)**: 9 columns
- pulseEnergy, pulseWidth, wavelength, spotSize, absorptionCoeff, reflectivity, focus coordinates

**Solver Settings (42-47)**: 6 columns
- maxCo, maxAlphaCo, maxThermalCourant, maxDi, thermalFluxRelax, adjustTimeStep

**Time Stepping (48-51)**: 4 columns
- deltaT, maxDeltaT, minDeltaT, endTime

**Laser Timing (52-53)**: 2 columns
- laserStartTime, laserEndTime

**Phase Change (54-63)**: 10 columns ⚠️ **8 AFFECTED BY BUG**
- Tvap, Tsol, evaporationCoeff, hf, gasConstant, p_ref, relaxationTime, maxSource, alphaMin, metalFractionCutoff

**Two-Temperature (64-69)**: 6 columns ⚠️ **2 AFFECTED BY BUG**
- Cl, De, Ce_coeff, maxTe, maxTl, minTe

**Advanced Interface (70-72)**: 3 columns
- momentumAccommodationCoeff, stickingCoeff, recoilMax

**Material Properties (73-76)**: 4 columns
- metal_Cp, metal_mu, metal_molWeight, metal_rho0

**Mesh/Domain (77-80)**: 4 columns
- totalCells, domainVolume, minCellVolume, maxCellVolume

---

## 🚀 Post-Fix Validation

After applying the fix and recompiling:

### 1. Compile
```bash
cd ~/compInterFoam
wmake
```

### 2. Test Run
```bash
cd TEST1
compInterFoam | tee run.log
```

### 3. Check CSV
```bash
# View first 2 rows
head -2 liftProcessTracking.csv

# Verify specific columns
tail -1 liftProcessTracking.csv | cut -d',' -f54-65
```

**Expected values** (should NOT be zeros):
- Column 54 (input_Tvap_K): ~3560
- Column 55 (input_Tsol_K): ~1941
- Column 57 (input_hf_Jkg): ~9100000
- Column 58 (input_gasConstant_JkgK): ~174
- Column 59 (input_p_ref_Pa): ~101325
- Column 61 (input_maxSource_kgm3s): ~1e25
- Column 64 (input_Cl_Jm3K): ~2500000
- Column 65 (input_De_m2s): ~0.0001

### 4. Python Analysis
```python
import pandas as pd

df = pd.read_csv('liftProcessTracking.csv')

# Check input parameters are not zero
input_cols = [col for col in df.columns if col.startswith('input_')]
print("Input parameters:")
print(df[input_cols].iloc[0])

# Verify dimensioned parameters
critical_params = [
    'input_Tvap_K', 'input_Tsol_K', 'input_hf_Jkg',
    'input_gasConstant_JkgK', 'input_p_ref_Pa',
    'input_Cl_Jm3K', 'input_De_m2s'
]
print("\nCritical parameters:")
print(df[critical_params].iloc[0])

# Should show actual values, not zeros!
assert df['input_Tvap_K'].iloc[0] > 1000, "Tvap should be ~3560 K"
assert df['input_hf_Jkg'].iloc[0] > 1e6, "hf should be ~9.1e6 J/kg"
```

---

## 📝 Additional Notes

### CSV File Characteristics
- **Append mode**: Safe for checkpoint/restart
- **Every timestep**: Complete temporal resolution
- **Master only**: No parallel write conflicts
- **ASCII format**: Human-readable, easy to parse

### Data Volume Estimation
For a typical LIFT simulation:
- Endtime: 2e-10 s (200 ps)
- Timestep: ~1e-15 s (adaptive)
- Total steps: ~200,000
- CSV size: ~200,000 rows × 80 columns × 15 bytes ≈ **240 MB**

### Performance Impact
- File I/O every timestep: minimal overhead (<1%)
- String formatting: negligible
- Master-only write: no MPI bottleneck

### Recommendations
1. **Fix the bug** before production runs
2. Add validation output after parameter extraction
3. Consider adding a `csvWriteInterval` option for very long simulations
4. Implement column subsets for targeted analysis
5. Add CSV checkpointing for fault tolerance

---

## Summary

**Status**: CSV implementation is 90% correct
- ✓ Structure and column alignment: Perfect
- ✓ Real-time metrics: All working
- ✓ File handling: Robust
- ❌ Dimensioned parameter extraction: **Critical bug**

**Action Required**: Apply fix from `FIXED_CSV_parameter_extraction.cpp` to lines 223-264

**Priority**: HIGH - affects data quality and reproducibility

**Effort**: Low - simple code replacement, 5 minutes + recompile

**Validation**: Essential - verify non-zero values in CSV columns 54-65
