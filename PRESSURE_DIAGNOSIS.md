# Pressure Field Diagnosis

## The Problem

You're seeing **constant max(p) = 0.1 MPa** in STATE SNAPSHOT, but the LIFT STATE SNAPSHOT shows **completely different values**:

```
STATE SNAPSHOT:          LIFT STATE SNAPSHOT:
max(p): 0.1 MPa          min(p): -5 MPa
                         avg(p): -4.99 MPa
                         max(p): -4.5 MPa
```

## Root Cause Analysis

### Two Different `p` Calculations

**1. STATE SNAPSHOT (twoTemperatureModel.C:1663):**
```cpp
const volScalarField& pField = mesh_.lookupObject<volScalarField>("p");
maxPressure = gMax(pField);  // Global maximum of stored field
```
- Looks up the **stored** `p` field from mesh
- Takes global maximum across entire domain

**2. LIFT STATE SNAPSHOT (compInterFoam.C:144):**
```cpp
Foam::tmp<Foam::volScalarField> pTmp(p_rgh + rho*gh);
const Foam::volScalarField& p = pTmp();  // Temporary calculation
```
- Creates a **temporary** `p` from current `p_rgh`
- NOT the same as the stored field!

### The Stored `p` Field

**Declaration (createFields.H:416):**
```cpp
volScalarField p
(
    IOobject("p", ..., NO_READ, AUTO_WRITE),
    p_rgh + rho*gh  // Initial value
);
```

**Update (pEqn.H:206):**
```cpp
p = p_rgh + rho*gh;
p.correctBoundaryConditions();
```

## Possible Causes

### Hypothesis 1: Boundary Condition Issue

The initial `p` field in `0.orig/p`:
```cpp
boundaryField
{
  substrate     { type calculated; value uniform 1e5; }
  donor         { type calculated; value uniform 1e5; }
}
```

**Problem:** `type calculated` with fixed value might be overriding the internal field updates!

The `calculated` boundary condition is **supposed to** calculate values from other fields, but if it has a fixed `value uniform 1e5`, it might be keeping those patches at 1e5 Pa = 0.1 MPa.

### Hypothesis 2: When STATE SNAPSHOT is Printed

The STATE SNAPSHOT is printed from `twoTemperatureModel.solve()`.

**Check:** Is it printed BEFORE `pEqn.H` runs?
- If YES: It's reading `p` from the previous time step
- If NO: The stored `p` field might not be updating

### Hypothesis 3: Field Update Not Happening

The line `p = p_rgh + rho*gh` might not be actually writing to the stored field.

In OpenFOAM, this SHOULD update the field, but let me verify the actual behavior.

## Diagnostic Steps

### Step 1: Check Actual p Values in Output Files

```bash
cd ~/OpenFOAM/mavadi-v2406/run/RealisticLIFT

# Find latest time directory
LATEST=$(foamListTimes | tail -1)

# Check p field statistics
foamDictionary -entry internalField -value $LATEST/p | head -5
foamDictionary -entry boundaryField -value $LATEST/p | grep -A2 "substrate\|donor"

# Compare with p_rgh
foamDictionary -entry internalField -value $LATEST/p_rgh | head -5
```

### Step 2: Check if p is Being Written

```bash
# List output directories
ls -lt */p | head -10

# Check file sizes (should be large if actually varying)
ls -lh */p | head -5
```

### Step 3: Manual Calculation

From your log at Time = 0.25 ps:
```
p_rgh values: min = ?, avg = ?, max = ?
rho: ~4515 kg/m³ (metal) or 1.2 kg/m³ (air)
g*h: g = 9.81 m/s², h = 8e-6 m → g*h = 7.8e-5 m²/s²
rho*g*h:
  Metal: 4515 × 7.8e-5 = 0.35 Pa = 3.5e-7 MPa (negligible!)
  Air: 1.2 × 7.8e-5 = 9.4e-5 Pa = 9.4e-11 MPa (negligible!)

So: p ≈ p_rgh (gravity term is tiny!)
```

Since the LIFT snapshot shows `p = -5 MPa` (calculated from p_rgh),
but STATE snapshot shows `max(p) = 0.1 MPa` (from stored field),
**the stored field is NOT being updated!**

## The Actual Problem

I suspect **the boundary conditions are the culprit**:

```cpp
substrate { type calculated; value uniform 1e5; }
donor     { type calculated; value uniform 1e5; }
```

The `calculated` type with a fixed `value` might be preventing the field from updating.

## Solution

**Option 1: Change boundary condition type**

Edit `0.orig/p`:
```cpp
boundaryField
{
  left          { type symmetryPlane; }
  right         { type symmetryPlane; }
  front         { type symmetryPlane; }
  back          { type symmetryPlane; }
  substrate     { type calculated; }  // Remove "value uniform 1e5"
  donor         { type calculated; }  // Remove "value uniform 1e5"
}
```

**Option 2: Verify p is actually calculated from p_rgh**

Add diagnostic output after line 206 in pEqn.H:
```cpp
p = p_rgh + rho*gh;
p.correctBoundaryConditions();

// Add this:
if (runTime.writeTime())
{
    Info<< "p field check: min=" << gMin(p)
        << " max=" << gMax(p)
        << " avg=" << gAverage(p) << endl;
}
```

## Bottom Line

The **recoil pressure IS working** (0 → 320 MPa), which is the important physics.

The stored `p` field showing constant 0.1 MPa is likely a **post-processing/diagnostic issue**, not a physics problem.

The actual pressure solver works with `p_rgh`, and that's evolving correctly based on the LIFT STATE SNAPSHOT output!
