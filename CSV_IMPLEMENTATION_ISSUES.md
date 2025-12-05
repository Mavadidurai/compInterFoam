# Critical Issues in CSV Implementation for compInterFoam

## Issue 1: Incorrect Extraction of Dimensioned Parameters ⚠️

### Location: `compInterFoam.C` lines 227-236

### Problem:
The code attempts to extract dimensioned parameters (like `Tvap`, `Tsol`, `hf`, etc.) using `lookupOrDefault<Foam::scalar>()`, but these values are defined as **dimensioned scalars** in the controlDict:

```cpp
// In controlDict:
Tvap    [0 0 0 1 0 0 0] 3560;  // This is a dimensionedScalar
Tsol    [0 0 0 1 0 0 0] 1941;
hf      [0 2 -2 0 0 0 0] 9.1e6;
gasConstant [0 2 -2 -1 0 0 0] 174;
p_ref   [1 -1 -2 0 0 0 0] 101325;

// But the code tries to read them as plain scalars:
input_Tvap = phaseChangeDict.lookupOrDefault<Foam::scalar>("Tvap", 0.0);  // WRONG!
input_Tsol = phaseChangeDict.lookupOrDefault<Foam::scalar>("Tsol", 0.0);  // WRONG!
```

### Why This Fails:
- OpenFOAM stores dimensioned values as `dimensionedScalar` objects
- Reading with `lookupOrDefault<Foam::scalar>()` will either:
  1. Fail to find the entry (since it's stored with dimensions)
  2. Return the default value (0.0) instead of the actual value
  3. Throw a type conversion error

### Affected Parameters in `phaseChangeCoeffs`:
- `input_Tvap` (line 227)
- `input_Tsol` (line 228)
- `input_hf` (line 230)
- `input_gasConstant` (line 231)
- `input_p_ref` (line 232)
- `input_maxSource` (line 234)

### Affected Parameters in `twoTemperatureProperties`:
- `input_Cl` (line 243)
- `input_De` (line 244)

### Impact:
- **All dimensioned input parameters will be recorded as 0.0 in the CSV file**
- Your CSV data will be missing critical configuration information
- Analysis and reproducibility will be compromised

## Issue 2: Potential Header/Data Column Mismatch

### Location: `compInterFoam.C` lines 1179-1348

### Verification Needed:
Count the header columns vs. data columns to ensure they match exactly.

**Header columns (lines 1179-1217):**
- Basic metrics: ~30 columns
- Laser properties: 9 columns
- Solver settings: 6 columns
- Time stepping: 4 columns
- Laser timing: 2 columns
- Phase change: 10 columns
- Two-temperature: 6 columns
- Advanced interface: 3 columns
- Material properties: 4 columns
- Mesh/Domain: 4 columns
**Total: ~78 columns**

**Data columns (lines 1260-1348):**
Must match exactly!

---

## Solution: Proper Extraction of Dimensioned Values

### Method 1: Read as dimensionedScalar, then extract value
```cpp
// CORRECT approach:
const dimensionedScalar Tvap_dim =
    phaseChangeDict.lookupOrDefault<dimensionedScalar>(
        "Tvap",
        dimensionedScalar("Tvap", dimTemperature, 0.0)
    );
input_Tvap = Tvap_dim.value();
```

### Method 2: Use lookupOrDefaultBackwardsCompatible (if available)
```cpp
input_Tvap = phaseChangeDict.lookupOrDefaultBackwardsCompatible<scalar>(
    {"Tvap"},
    0.0
);
```

### Method 3: Try-catch with dimension stripping
```cpp
try
{
    dimensionedScalar Tvap_dim(phaseChangeDict.lookup("Tvap"));
    input_Tvap = Tvap_dim.value();
}
catch (...)
{
    input_Tvap = phaseChangeDict.lookupOrDefault<scalar>("Tvap", 0.0);
}
```

---

## Recommended Fix

Replace lines 223-237 in `compInterFoam.C`:

```cpp
// ===== Extract phase change parameters =====
if (controlDict.found("phaseChangeCoeffs"))
{
    const Foam::dictionary& phaseChangeDict = controlDict.subDict("phaseChangeCoeffs");

    // Dimensioned parameters - extract the scalar value from dimensionedScalar
    if (phaseChangeDict.found("Tvap"))
    {
        dimensionedScalar Tvap_dim(phaseChangeDict.lookup("Tvap"));
        input_Tvap = Tvap_dim.value();
    }
    if (phaseChangeDict.found("Tsol"))
    {
        dimensionedScalar Tsol_dim(phaseChangeDict.lookup("Tsol"));
        input_Tsol = Tsol_dim.value();
    }
    if (phaseChangeDict.found("hf"))
    {
        dimensionedScalar hf_dim(phaseChangeDict.lookup("hf"));
        input_hf = hf_dim.value();
    }
    if (phaseChangeDict.found("gasConstant"))
    {
        dimensionedScalar gasConstant_dim(phaseChangeDict.lookup("gasConstant"));
        input_gasConstant = gasConstant_dim.value();
    }
    if (phaseChangeDict.found("p_ref"))
    {
        dimensionedScalar p_ref_dim(phaseChangeDict.lookup("p_ref"));
        input_p_ref = p_ref_dim.value();
    }
    if (phaseChangeDict.found("maxSource"))
    {
        dimensionedScalar maxSource_dim(phaseChangeDict.lookup("maxSource"));
        input_maxSource = maxSource_dim.value();
    }

    // Dimensionless parameters - read directly as scalar
    input_evaporationCoeff = phaseChangeDict.lookupOrDefault<Foam::scalar>("evaporationCoeff", 0.0);
    input_relaxationTime = phaseChangeDict.lookupOrDefault<Foam::scalar>("relaxationTime", 1e-11);
    input_alphaMin = phaseChangeDict.lookupOrDefault<Foam::scalar>("alphaMin", 0.001);
    input_metalFractionCutoff = phaseChangeDict.lookupOrDefault<Foam::scalar>("metalFractionCutoff", 1e-6);
}
```

Similarly for `twoTemperatureProperties` (lines 239-264):

```cpp
// ===== Extract two-temperature model parameters =====
if (controlDict.found("twoTemperatureProperties"))
{
    const Foam::dictionary& twoTempDict = controlDict.subDict("twoTemperatureProperties");

    // Dimensioned parameters
    if (twoTempDict.found("Cl"))
    {
        dimensionedScalar Cl_dim(twoTempDict.lookup("Cl"));
        input_Cl = Cl_dim.value();
    }
    if (twoTempDict.found("De"))
    {
        dimensionedScalar De_dim(twoTempDict.lookup("De"));
        input_De = De_dim.value();
    }

    // Ce coefficient extraction (already correct)
    if (twoTempDict.found("Ce"))
    {
        const Foam::dictionary& CeDict = twoTempDict.subDict("Ce");
        if (CeDict.found("coeffs"))
        {
            Foam::List<Foam::Tuple2<Foam::scalar, Foam::scalar>> coeffs =
                CeDict.lookup("coeffs");
            if (coeffs.size() > 1)
            {
                input_Ce_coeff = coeffs[1].first();
            }
        }
    }

    // Dimensionless parameters
    input_maxTe = twoTempDict.lookupOrDefault<Foam::scalar>("maxTe", 20000.0);
    input_maxTl = twoTempDict.lookupOrDefault<Foam::scalar>("maxTl", 10000.0);
    input_minTe = twoTempDict.lookupOrDefault<Foam::scalar>("minTe", 200.0);
}
```

---

## Testing the Fix

After applying the fix:

1. **Recompile:**
   ```bash
   cd ~/compInterFoam
   wmake
   ```

2. **Run a test case:**
   ```bash
   cd TEST1
   compInterFoam
   ```

3. **Verify CSV output:**
   ```bash
   head -2 liftProcessTracking.csv
   ```

   Check that the input parameter columns contain actual values, not zeros:
   - `input_Tvap_K` should be ~3560, not 0
   - `input_Tsol_K` should be ~1941, not 0
   - `input_hf_Jkg` should be ~9.1e6, not 0
   - `input_gasConstant_JkgK` should be ~174, not 0
   - `input_p_ref_Pa` should be ~101325, not 0

4. **Validate with Python/Excel:**
   ```python
   import pandas as pd
   df = pd.read_csv('liftProcessTracking.csv')
   print(df[['input_Tvap_K', 'input_Tsol_K', 'input_hf_Jkg']].head())
   ```

---

## Additional Recommendations

1. **Add validation logging:**
   After extraction, print the values to verify they're correct:
   ```cpp
   Info<< "CSV: Extracted Tvap = " << input_Tvap << " K" << nl
       << "CSV: Extracted hf = " << input_hf << " J/kg" << nl
       << "CSV: Extracted p_ref = " << input_p_ref << " Pa" << endl;
   ```

2. **Error handling:**
   Wrap dimension extractions in try-catch blocks to handle missing entries gracefully

3. **Unit testing:**
   Create a simple test case to verify all CSV columns are populated correctly

---

## Summary

**Current Status:** ❌ CSV implementation has critical bugs
- Dimensioned parameters are being read incorrectly
- All dimensioned input values will be 0.0 in CSV output

**Required Action:** Fix parameter extraction using proper `dimensionedScalar` handling

**Priority:** HIGH - affects all CSV data quality and reproducibility
