# Laser Ablation Simulation - Recoil Pressure Fix

## Problem Diagnosis

### Critical Issue
Your simulation has a **4x mismatch** between kinetic theory predictions and pressure limits:

- **Kinetic Theory Request**: 320 MPa (at T_lattice ≈ 8000 K)
- **Pressure Limits**: 80 MPa
- **Result**: Flow velocity only 0.009% of expected (severely throttled)

### Symptoms
```
Recoil clamp requests: peak 320.02703 MPa vs limit 80 MPa
|U|_metal_max = 0.017 m/s (0.009% × recoil-limited 188 m/s)
PIMPLE: not converged within 3 iterations
```

## Root Cause

At lattice temperatures of ~8000 K, the **kinetic theory recoil model** calculates:

```
P_recoil ≈ ρ_vapor * k_B * T / m_atom ≈ 320 MPa
```

But your configuration limits recoil to **80 MPa** based on Feinaeugle et al. experimental data.

**This creates two problems:**
1. Physical: Flow cannot develop properly (artificial throttling)
2. Numerical: PIMPLE cannot converge due to inconsistent forcing

---

## Solution Options

### ✅ Option 1: Increase Pressure Limits (APPLIED)

**When to use:** If high temperatures (8000 K) are physically expected in your simulation.

**Changes Made:**
- `maxRecoilPressure`: 80 MPa → **350 MPa**
- `maxPressure`: 80 MPa → **350 MPa**
- `pressureClamp`: false → **true**

**Files Modified:**
- `system/controlDict` (lines 92, 102)
- `system/fvSolution` (lines 181-182)

**Rationale:** Allows kinetic theory to operate without artificial throttling.

**Trade-offs:**
- May produce pressures higher than experimental observations (Feinaeugle: ~80 MPa)
- Check if your material properties justify these extreme conditions

---

### 🔄 Option 2: Reduce Temperature Cap (ALTERNATIVE)

**When to use:** If you want to match experimental observations (~80 MPa) more closely.

**Changes Needed:**

1. **Reduce maximum lattice temperature** (controlDict):
   ```
   maxT     5000;   // was 8000 K
   maxTl    5000;   // was 8000 K
   ```

2. **Keep original pressure limits**:
   ```
   maxRecoilPressure   8e7;    // 80 MPa
   maxPressure         8e7;    // 80 MPa
   ```

**Rationale:**
- At T_l ≈ 5000 K, kinetic theory predicts ~100 MPa (closer to 80 MPa limit)
- Better matches Feinaeugle experimental observations
- Ti vaporization temperature is 3560 K, so 5000 K allows reasonable superheat

**How to apply Option 2:**
```bash
cd /home/user/compInterFoam/TEST1

# Revert Option 1 changes
sed -i 's/3.5e8/8e7/g' system/controlDict
sed -i 's/3.5e8/8e7/g' system/fvSolution

# Apply temperature reduction
sed -i 's/maxT            8000/maxT            5000/' system/controlDict
sed -i 's/maxTl           8000/maxTl           5000/' system/controlDict
```

---

### 🔬 Option 3: Tune Kinetic Theory Parameters (ADVANCED)

**When to use:** If you want to keep T=8000K but reduce recoil via physics.

**Changes Needed** (controlDict, advancedInterfaceCapturing):
```
stickingCoeff       0.18  →  0.50   // Less aggressive vaporization
momentumAccom       0.18  →  0.40   // Reduced momentum transfer
```

**Rationale:**
- Higher sticking coefficient = less vapor flux = lower recoil
- Typical values: 0.1-0.5 for metals
- Current 0.18 may be too aggressive for your laser parameters

---

## Recommended Next Steps

1. **Resume simulation** with Option 1 (already applied):
   ```bash
   cd /home/user/compInterFoam/TEST1
   fg  # Resume stopped job
   # OR restart:
   compInterFoam > log.compInterFoam 2>&1 &
   ```

2. **Monitor diagnostics:**
   ```bash
   tail -f log.compInterFoam | grep -E "Recoil|Velocity diagnostic|max\(Tl\)"
   ```

3. **Check for improvements:**
   - Recoil clamp warnings should disappear
   - Metal velocity should increase from 0.017 m/s to ~50-100 m/s
   - PIMPLE should converge within 3 iterations

4. **Validate results:**
   - Check if flow velocities are physically reasonable
   - Compare ablation depths/transfer with literature
   - If pressures seem too high, try Option 2 or 3

---

## Physics Considerations

### Feinaeugle et al. Observations
- Ti fs-LIFT recoil: **~80 MPa plateau**
- At T ≈ 6600 K
- 200 fs, 343 nm laser

### Your Current Simulation
- T_lattice reaching **8000 K**
- Kinetic theory: **320 MPa**
- This suggests either:
  - **Temperature is too high** → reduce maxTl (Option 2)
  - **Kinetic coefficients too aggressive** → increase sticking/accom (Option 3)
  - **Expected for your parameters** → keep Option 1

### Physical Intuition
- Recoil pressure scales as: P ∝ T^(3/2) * exp(-E_vap/kT)
- Going from 6600K to 8000K increases P by factor of ~2-3
- So 80 MPa → 200-240 MPa is reasonable, but 320 MPa suggests:
  - Over-prediction by kinetic theory, OR
  - Legitimate extreme conditions in your case

---

## Troubleshooting

### If simulation still doesn't converge:
1. Reduce time step: `maxDeltaT 1e-13;` → `5e-14;`
2. Tighten Courant numbers: `maxCo 0.1;` → `maxCo 0.05;`
3. Add more PIMPLE correctors: `nOuterCorrectors 3;` → `5;`

### If velocities seem unrealistic:
1. Try Option 2 (temperature reduction)
2. Check laser power: 30 nJ might be too high for this geometry
3. Verify material properties (absorption coefficient, thermal properties)

---

## Files Modified

### Option 1 (Applied):
- `system/controlDict` - maxRecoilPressure: 80→350 MPa
- `system/fvSolution` - maxPressure: 80→350 MPa, pressureClamp: false→true

### Backups:
Original files are in git history. To revert:
```bash
git checkout HEAD -- system/controlDict system/fvSolution
```

---

## Summary

**Current State:** Option 1 applied - pressure limits increased to 350 MPa

**Expected Outcome:**
- Recoil warnings should stop
- Flow velocity should increase dramatically
- Simulation should continue progressing

**Decision Point:** After running for a few timesteps, evaluate if:
- ✅ Results are physically reasonable → continue
- ❌ Pressures/velocities too extreme → try Option 2 or 3

---

Generated: 2025-11-07
Case: /home/user/compInterFoam/TEST1
