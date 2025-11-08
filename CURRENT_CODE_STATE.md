# Current Code State - Complete Analysis

## 📁 File: alphaSuSp.H

### Current Implementation (Lines 60-65)
```cpp
// Store diagonal magnitude explicitly positive for stable implicit treatment
// Following Hardt & Wondra (J. Comput. Phys. 227, 2008)
const scalar diagMag = evapRate*evapWeight + condRate*condWeight;

SpField[celli] = diagMag;  // Positive diagonal strengthens matrix
SuField[celli] = -rate + diagMag*a1;  // Explicit source term
```

**Status:** ✅ CORRECT - Using explicit positive `diagMag`

### How It Works
```cpp
evapRate = 0.5*(rate + |rate|)  // ≥ 0 (evaporation component)
condRate = 0.5*(|rate| - rate)  // ≥ 0 (condensation component)

evapWeight = α₂·ρ₂ / ρₘᵢₓ  // > 0
condWeight = α₁·ρ₁ / ρₘᵢₓ  // > 0

diagMag = evapRate*evapWeight + condRate*condWeight  // Always ≥ 0
```

**Result:** Sp is always positive → strengthens diagonal → stable matrix ✓

---

## 📁 File: alphaEqn.H

### Section 1: Implicit Solver (Lines 171-259)

#### Phase Change Source Terms (Lines 188-199)
```cpp
if (MULESCorr)
{
    fvScalarMatrix alpha1Eqn(...);

    // Apply phase change diagonal (with compressibility if available)
    if (divU.valid())
    {
        alpha1Eqn -= fvm::Sp(Sp() + divU(), alpha1);
    }
    else
    {
        alpha1Eqn -= fvm::Sp(Sp(), alpha1);  // ✅ FIXED: was += (wrong!)
    }

    // Apply phase change explicit source (common to both cases)
    alpha1Eqn -= Su();  // ✅ FIXED: moved outside, was duplicated
```

**Status:** ✅ CORRECT - Consistent -= signs in both branches

#### Diagonal Stabilization (Lines 210-226)
```cpp
// DISABLED: Stabilizing diagonal was freezing interface motion at femtosecond timescales
// The diagThreshold = SMALL/dt becomes huge (0.3) with dt = 1e-14 s, preventing
// alpha.metal from evolving despite sufficient driving forces (recoil, velocity).
// With the MULES limiter and other safeguards, the matrix should remain stable.
/*
    if (minDiag.value() < diagThreshold.value())
    {
        alpha1Eqn += fvm::Sp(diagThreshold, alpha1);
        ...
    }
*/
```

**Status:** ✅ DISABLED (commented out) - Was freezing alpha evolution

#### Singular Matrix Fallback (Lines 229-239)
```cpp
SolverPerformance<scalar> solverPerf = alpha1Eqn.solve(alpha1.name());

if (solverPerf.singular())
{
    alpha1 = alpha1.oldTime();  // ⚠️ REVERTS TO OLD VALUE!
    alpha1.correctBoundaryConditions();

    if (verbose && master)
    {
        Info<< "alpha1 solve detected singular matrix; reverting to previous"
            << " alpha field" << endl;
    }
}
```

**Status:** ⚠️ ACTIVE - This could explain constant alpha if triggering!

**CRITICAL QUESTION:** Is this fallback triggering in your simulation?

---

### Section 2: MULES Correction Loop (Lines 261-360)

#### MULES Settings from fvSolution
```cpp
MULESCorr       true;      // ✅ ENABLED
nAlphaCorr      5;         // ✅ 5 correction iterations
nAlphaSubCycles 5;         // ✅ 5 subcycles
```

#### MULES::correct Implementation (Lines 278-306)

**With divU (compressible):**
```cpp
if (divU.valid())
{
    MULES::correct
    (
        geometricOneField(),
        alpha1,
        talphaPhi1Un(),
        talphaPhi1Corr.ref(),
        Sp(),           // ✅ HAS phase change
        Su(),           // ✅ HAS phase change
        oneField(),
        zeroField()
    );
}
```

**Without divU (incompressible):**
```cpp
else
{
    // FIX: Include phase change source terms in incompressible case
    MULES::correct
    (
        geometricOneField(),
        alpha1,
        talphaPhi1Un(),
        talphaPhi1Corr.ref(),
        Sp(),           // ✅ FIXED: was oneField()
        Su(),           // ✅ FIXED: was zeroField()
        oneField(),
        zeroField()
    );
}
```

**Status:** ✅ FIXED - Both branches now include phase change

#### MULES::explicitSolve Implementation (Lines 322-350)

**With divU (compressible):**
```cpp
if (divU.valid())
{
    MULES::explicitSolve
    (
        geometricOneField(),
        alpha1,
        phiAlphaRef,
        alphaPhi10,
        Sp(),
        (Su() + divU()*min(alpha1(), scalar(1)))(),  // Combined source
        oneField(),
        zeroField()
    );
}
```

**Without divU (incompressible):**
```cpp
else
{
    // FIX: Include phase change source terms in incompressible case
    MULES::explicitSolve
    (
        geometricOneField(),
        alpha1,
        phiAlphaRef,
        alphaPhi10,
        Sp(),           // ✅ FIXED: was oneField()
        Su(),           // ✅ FIXED: was zeroField()
        oneField(),
        zeroField()
    );
}
```

**Status:** ✅ FIXED - Both branches now include phase change

---

## 🔍 Code Flow Analysis

### What Happens Each Time Step

```
1. Calculate phase change source terms (alphaSuSp.H)
   → Sp() = diagMag (positive) ✓
   → Su() = -rate + diagMag*a1 ✓

2. IF MULESCorr = true (IT IS):

   a) Build implicit matrix equation:
      ∂α/∂t + ∇·(φα) - Sp()*α - Su() = 0

   b) Solve implicit equation:
      solverPerf = alpha1Eqn.solve()

   c) CHECK IF SINGULAR: ⚠️ THIS IS THE KEY!
      if (solverPerf.singular())
      {
          alpha1 = alpha1.oldTime();  // DISCARD SOLUTION!
      }

   d) Run MULES correction loop (5 times):
      - MULES::correct with Sp(), Su()
      - Apply limiters
      - Correct boundaries

3. Update alpha2 = 1 - alpha1
```

---

## ⚠️ THE CRITICAL QUESTION

**Why is alpha constant despite solver iterating?**

Three possible reasons:

### Hypothesis 1: Singular Matrix Fallback Triggering
**Lines 229-239:** If `solverPerf.singular()` returns true, alpha is reverted to oldTime every step.

**How to check:**
- Look for message: "alpha1 solve detected singular matrix; reverting to previous alpha field"
- This requires `verbose = true` in the solver

**If this is happening:**
- The implicit solve finds a solution (2 iterations)
- But OpenFOAM thinks matrix is singular
- Discards the solution
- Alpha stays frozen

### Hypothesis 2: Changes Below Output Precision
**Expected change per step:** ~1e-8

**Output precision:** Typically 8 significant figures

If alpha = 0.40213438, change of 1e-8 would give:
```
0.40213438 → 0.40213437 (might not show in output)
```

### Hypothesis 3: MULES Limiting Too Strong
MULES applies bounds [0, 1] and limiters. With very small changes, limiters might prevent evolution.

From line 353:
```cpp
alpha1 = Foam::max(Foam::min(alpha1, scalar(1)), scalar(0));
```

---

## 🎯 Diagnostic Steps

### Step 1: Check if Singular Fallback is Triggering

**In your log output, look for:**
```
"alpha1 solve detected singular matrix; reverting to previous alpha field"
```

**If present:** The fallback is the problem, not MULES being disabled

**If absent:** Then alpha might be changing but too small to see

### Step 2: Check Verbose Output

Your simulation needs to run with `verbose = true` to see detailed messages.

Check if you see:
```
"alpha1 solve: initial residual = ..., final residual = ..., iterations = ..."
```

### Step 3: Check Matrix Diagonal

The code calculates (lines 208):
```cpp
const dimensionedScalar minDiag = gMin(alpha1Eqn.diag());
```

This is computed but not printed. At dt=1e-14, the diagonal should be:
```
diag ≈ ρ/dt + Sp() + divU()
     ≈ 4515/(1e-14) + diagMag + divU()
     ≈ 4.5e17 + small terms
```

Matrix should be VERY non-singular!

---

## 📊 Current Status Summary

| Component | Status | Evidence |
|-----------|--------|----------|
| **alphaSuSp.H Sp sign** | ✅ FIXED | Using positive diagMag |
| **alphaEqn.H implicit signs** | ✅ FIXED | Consistent -= operators |
| **MULES::correct phase change** | ✅ FIXED | Sp(), Su() in both branches |
| **MULES::explicitSolve phase change** | ✅ FIXED | Sp(), Su() in both branches |
| **MULES enabled** | ✅ YES | MULESCorr = true |
| **Diagonal stabilization** | ✅ DISABLED | Commented out (was freezing) |
| **Singular fallback** | ⚠️ ACTIVE | Could be triggering! |
| **Solver iterations** | ✅ WORKING | 2 iterations (not 0) |
| **Residual convergence** | ✅ WORKING | 1e-14 → 1e-78 |

---

## 🔬 Bottom Line

**All phase change fixes are correctly in place:**
- ✅ Positive Sp diagonal
- ✅ Correct -= signs
- ✅ MULES has phase change source terms
- ✅ MULES is enabled (not disabled)

**But there's a FALLBACK mechanism that could freeze alpha:**

Lines 229-239 will **discard the solution** if the matrix is detected as singular, even if the solver successfully iterates!

**Next step:** Check your log output for the singular matrix warning message. This will tell us if the fallback is the culprit.

---

**Files reviewed:**
- `/home/user/compInterFoam/alphaSuSp.H` - Lines 1-74
- `/home/user/compInterFoam/alphaEqn.H` - Lines 1-401

**Last commit:** 39457c1 - Clean up redundant Su() application in alphaEqn.H
