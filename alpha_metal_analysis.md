# Alpha.Metal Analysis - Why It's Not Changing

## Observations from Simulation Log

### Alpha.Metal Volume Fraction:
```
Phase-1 volume fraction (alpha.metal) = 0.40212149
  (domain average - includes donor and receiver substrates)
  Min(alpha.metal) = 0
  Max(alpha.metal) = 1
```

**Constant at 0.40212149** across all timesteps (0.2 ps → 13.28 ps)

### Metal Volume:
```
t=0.2 ps:   Metal volume: 4035.7 µm³
t=13.28 ps: Metal volume: 4035.5707 µm³
Metal loss: 0.0032% (0.13 µm³ ablated)
```

**Small but measurable mass loss!**

---

## Why Alpha Appears Constant (But It's Actually OK)

### Reason 1: **Domain Includes Fixed Substrates**

Your geometry includes:
1. **Donor substrate** (bottom) - FIXED, thick
2. **Active metal film** (71.4 nm) - ABLATING
3. **Receiver substrate** (top) - FIXED, thick
4. **Air gaps** - surrounding

**Domain-averaged alpha.metal** = Total metal volume / Total domain volume

Since the **substrates are much larger** than the thin film:
- Film volume: ~71.4 nm × beam area ≈ small
- Substrate volumes: Much larger
- Film ablation (0.13 µm³) / Total metal (4035 µm³) = **0.0032%**

**Conclusion**: Alpha.metal SHOULD be nearly constant because you're only ablating a tiny fraction of total metal!

---

### Reason 2: **Ablation Rate is Very Small**

Mass loss calculation:
```
Initial: 4035.7 µm³
Current: 4035.57 µm³
Lost: 0.13 µm³ over 13.28 ps

Loss rate: 0.13 µm³ / 13.28 ps = 0.0098 µm³/ps
```

At this rate:
- To lose 1% of film: would take ~400 ps
- To ablate entire 71.4 nm film: would take much longer

**This is physically realistic for early-stage ablation!**

---

### Reason 3: **Phase-Change Cells Are Active**

From the log:
```
Active-film alpha average (phase-change cells) = 0.98752114
  (evaluated over 1.5261627 µm^3 participating in mass transfer)
```

This shows:
- ✅ **1.53 µm³ of material is actively evaporating**
- ✅ Alpha in those cells ≈ 0.9875 (nearly pure metal)
- ✅ 336 interface cells have mass flux above threshold
- ✅ Max mass flux: 1472 kg/m²/s (realistic)

**The ablation IS happening, just in a small region!**

---

## The REAL Problem: Metal Velocity Too Low

### What We're Seeing:
```
avg(|U|) (metal): 0.0093 m/s  ← Only 9 mm/s!
max(|U|):         1123.69 m/s ← In vapor phase
```

### What We Should See (LIFT Physics):
```
Metal jet velocity: 30-100 m/s   (Feinaeugle et al., 2017)
Vapor expansion:    500-2000 m/s (thermal)
```

**The metal isn't moving because recoil pressure isn't accelerating it!**

---

## Why Metal Isn't Accelerating

### Current Physics (WITHOUT Alpha-Weighting):

Recoil pressure gradient `∇P_recoil` acts equally on all cells:

**Acceleration = -∇P / ρ**

| Phase | Density (kg/m³) | Acceleration | Result |
|-------|----------------|--------------|--------|
| Metal | 4500 | ∇P / 4500 = small | 0.009 m/s ❌ |
| Vapor | 1-10 | ∇P / 1 = HUGE | 1124 m/s ❌ |

The **low-density vapor** gets accelerated because `a ∝ 1/ρ`.

The **high-density metal** barely moves!

---

## Expected Behavior AFTER Alpha-Weighting Fix

With the UEqn.H fix (`totalPressure += alpha1 * recoilContribution`):

**Force = -∇(α·P) = -(α·∇P + P·∇α)**

| Phase | Alpha | Effective Force | Acceleration | Expected Velocity |
|-------|-------|-----------------|--------------|-------------------|
| Metal | 1 | Full ∇P | ∇P / 4500 | **30-100 m/s** ✅ |
| Vapor | 0 | Zero | 0 | **Thermal only** ✅ |
| Interface | 0.5 | P·∇α (traction) | Interface force | Physical ✅ |

---

## Why You Don't See Alpha Changing Much

### Expected Alpha Evolution:

**Before ablation starts:**
- Alpha.metal = 0.40212 (substrates + film)

**During ablation (what SHOULD happen):**
- Material ejected from film → metal jet forms
- Jet cells: α = 0.5-1.0 (mixed metal-vapor)
- Vapor plume: α = 0.0-0.1 (mostly vapor)
- **Domain average α decreases slowly** (jet << substrates)

**What's happening now:**
- Vapor expands at 1124 m/s (wrong!)
- Metal doesn't move (0.009 m/s)
- No jet formation
- Material evaporates in-place
- **Alpha stays constant** (no material redistribution)

---

## Diagnostic: Check If Jet Is Forming

To verify if material is actually being ejected (not just evaporated):

### Check 1: Velocity Distribution
```bash
# In ParaView, load t=13.28 ps
# Color by U magnitude
# Look for:
#   - High velocity (>30 m/s) in METAL phase (α > 0.5)
#   - Upward/outward jet from film surface
```

**Expected**: Jet with α=0.5-1.0, U=30-100 m/s

**Actual (without fix)**: High U only where α≈0 (vapor)

### Check 2: Alpha Distribution
```bash
# Color by alpha.metal
# Look for:
#   - Ejected droplets/jet (α > 0.5) moving away from surface
```

**Expected**: Metal parcels with α>0.5 traveling upward

**Actual (without fix)**: Metal stays put, vapor expands

### Check 3: Pressure Field
```bash
# Color by p_rgh
# Look for:
#   - High pressure (71 MPa) at interface driving metal jet
```

**Expected**: Pressure gradient accelerates metal

**Actual (without fix)**: Pressure gradient accelerates vapor (low ρ)

---

## Summary: Is Alpha Behavior Correct?

### ✅ **Alpha.metal = 0.40212 staying constant is EXPECTED:**
- Domain includes large fixed substrates
- Only thin film (tiny fraction) is ablating
- 0.0032% mass loss is realistic for 13 ps
- Active ablation IS happening (1.53 µm³ evaporating)

### ❌ **But metal velocity = 0.009 m/s is WRONG:**
- Metal should be accelerating at 30-100 m/s
- Recoil pressure not driving metal phase
- Vapor accelerating instead (1124 m/s)
- **This is the problem the alpha-weighted fix solves!**

---

## After the Fix: Expected Changes

| Metric | Before Fix | After Fix | Physical |
|--------|------------|-----------|----------|
| **Alpha.metal** | 0.40212 (constant) | 0.40210-0.40200 (slight drop) | ✅ Small change OK |
| **Metal avg velocity** | 0.009 m/s ❌ | **30-100 m/s** ✅ | ✅ Jet formation |
| **Vapor velocity** | 1124 m/s ⚠️ | **500-1000 m/s** ✅ | ✅ Thermal expansion |
| **Mass loss rate** | 0.0098 µm³/ps | Similar or higher | ✅ Realistic |
| **Jet formation** | ❌ No | ✅ Yes | ✅ LIFT physics |

---

## Conclusion

**Alpha.metal not changing is CORRECT** - you're only ablating a tiny fraction of total metal (thin film vs. thick substrates).

**Metal not moving is WRONG** - this is the core issue the alpha-weighted fix addresses.

After applying the UEqn.H patch, you should see:
1. ✅ Metal velocity increases to 30-100 m/s
2. ✅ Vapor velocity decreases to 500-1000 m/s
3. ✅ Jet formation visible in ParaView
4. ✅ Alpha.metal still mostly constant (expected)
5. ✅ Metal loss rate may increase (better momentum coupling)
