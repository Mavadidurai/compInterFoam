# LIFT Process Tracker v2.0 - Diagnostics Guide

## Overview

The enhanced LIFT progress tracker provides real-time physics-based diagnostics to monitor and understand the laser-induced forward transfer process. This guide explains what factors determine progress and how to diagnose issues.

---

## Progress Determination Factors

### 1. **Phase-Based Progress (16 Phases)**

The tracker monitors 17 distinct phases (0-16) of the LIFT process:

| Phase | Name | Trigger Conditions | Physical Meaning |
|-------|------|-------------------|------------------|
| 0 | IDLE | Initial state | Simulation started, no laser |
| 1 | LASER_ACTIVE | peakQLaser > 1e6 W/m³ | Laser is ON and depositing energy |
| 2 | ABSORPTION | Te > 500 K | Electrons absorbing laser energy |
| 3 | HEATING | Tl > 500 K | Lattice heating from e-l coupling |
| 4 | PRE-MELT | Tl > 1500 K | Approaching melting point |
| 5 | MELTING | Tl > 1941 K | Iron melting (phase transition) |
| 6 | VAPORIZING | Tl > 3560 K | Boiling/vaporization begins |
| 7 | PLASMA | Tl > 5000 K AND p > 10 MPa | High-density plasma formation |
| 8 | RECOIL | Recoil > 10 MPa | Significant recoil pressure |
| 9 | EJECTION | \|U\| > 10 m/s | Film starts moving |
| 10 | ACCELERATION | \|U\| > 100 m/s | Rapid acceleration phase ⬅️ **YOU ARE HERE** |
| 11 | TRANSFER | Metal loss > 2% | Significant material transfer |
| 12 | COOLING | Tl decreasing + t > 50 ps | Thermal decay |
| 13 | CONDENSING | Tl < 3560 K (cooling) | Below vaporization |
| 14 | SOLIDIFYING | Tl < 1941 K (cooling) | Re-solidification |
| 15 | NEAR_AMBIENT | Tl < 500 K | Approaching initial temp |
| 16 | EQUILIBRIUM | ΔT < 200 K | Thermally stable |

**Current Progress Calculation:**
```
Progress (%) = (Current Phase / 16) × 100
```

At Phase 10 (ACCELERATION), you're at **62.5% complete**.

---

## Key Diagnostic Indicators

### 🌡️ **Thermal Diagnostics**

| Parameter | Current Value | Significance | Diagnostic Interpretation |
|-----------|---------------|--------------|--------------------------|
| **Te max** | 10,600 K | Electron temperature peak | Normal for femtosecond laser |
| **Tl max** | 10,000 K | Lattice temperature (clamped) | At safety limit - may indicate intense heating |
| **Tl avg** | 351 K | Domain average | Most material still cool |
| **Temp spread** | 9,700 K | Max - Min temperature | Large gradient → localized heating |

**Diagnosis:**
- ✅ High Te: Laser absorption working correctly
- ⚠️ Tl at clamp (10,000 K): May need to check if physical or numerical
- ✅ Low Tl avg: Heating localized to laser spot (correct)

---

### 💨 **Pressure & Recoil Diagnostics**

| Parameter | Current Value | Physical Meaning | What It Tells You |
|-----------|---------------|------------------|-------------------|
| **Max pressure** | 41,193 MPa (41.2 GPa) | Total pressure | Extreme compression - typical for LIFT |
| **Recoil pressure** | 37,448 MPa | Vapor recoil force | Driving film acceleration |
| **Peak recoil seen** | Tracked automatically | Historical maximum | Used to detect separation |
| **Pressure gradient** | Calculated | Shock wave indicator | High values → shock formation |

**Diagnosis for Your Simulation:**
- ✅ **Recoil = 37.4 GPa**: This is **enormous** and physically realistic for femtosecond LIFT
- ✅ **Pressure ratio** (recoil/total = 91%): Most pressure from vapor, not just thermal expansion
- ℹ️ **What this means**: Vapor layer is **actively pushing** your film

**Key Threshold:**
- When `maxRecoil < 0.5 × maxRecoilSeen` AND `avgVel > 50 m/s` → **SEPARATION DETECTED**

---

### 🎯 **Metal Film Status**

| Metric | Current | What It Means | Action/Diagnosis |
|--------|---------|---------------|------------------|
| **Volume** | 4,035.62 µm³ | Film volume | Track for conservation |
| **Volume loss** | 0.0018% | Material removed | **Still attached!** |
| **Loss rate** | Calculated real-time | Rate of removal | Increases at separation |
| **Interface cells** | Counted | Interface complexity | More cells → fragmentation |
| **Interface area** | Calculated | Surface area | Increases if breakup occurs |

**Your Current Status:**
- **0.0018% loss**: This is **NOT bulk ejection** - just surface evaporation
- **Film is still attached** to donor substrate
- **Acceleration in progress** but no separation yet

**Expected at Ejection:**
- Loss will **jump** to 0.1% → 1% → 10%+ rapidly
- Interface cells may drop (detachment) or spike (fragmentation)

---

### ⚡ **Velocity Analysis**

| Parameter | Current | Interpretation | Physics |
|-----------|---------|----------------|---------|
| **Film avg velocity** | 0.195 m/s | Bulk motion | Film moving as a unit |
| **Film max velocity** | 249.7 m/s | Fastest point | Surface/edge velocities |
| **Velocity std dev** | NEW | Velocity spread | Coherent (low σ) vs turbulent (high σ) |
| **Acceleration** | km/s² | Rate of speed increase | Tells when ejection will occur |
| **Turbulent KE** | NEW | Flow disorder | Energy not going into bulk motion |

**Your Current State:**
```
Avg velocity:  0.195 m/s  ← Bulk film is barely moving
Max velocity:  249.7 m/s  ← Surface is moving fast (6.1% of recoil-limited 4,073 m/s)
```

**Diagnosis:**
- ⚠️ **Large velocity discrepancy**: Surface moving 1000× faster than bulk
- **Physical Meaning**:
  - Recoil is accelerating the **surface layer**
  - Bulk film hasn't fully responded yet
  - Vapor layer building between film and donor

**Expected Evolution:**
1. Avg velocity will **rapidly increase** as bulk film accelerates
2. When avg → 100-500 m/s: **EJECTION IMMINENT**
3. Max/avg ratio decreases as film moves coherently

---

### ☁️ **Vapor Plume Tracking** (NEW)

| Metric | What It Tracks | Why It Matters |
|--------|----------------|----------------|
| **Plume volume** | Volume of high-density gas | Measures vapor production |
| **Plume cells** | Number of vapor cells | Plume extent |
| **Plume avg velocity** | Average vapor expansion | Typical plume speed |
| **Plume max velocity** | Fastest vapor | Shock/explosion velocity |

**Current Values:**
- **Max plume velocity**: 2,047 m/s (from your log)
- This is **supersonic** (>>343 m/s in air)
- Indicates **explosive vapor expansion**

**Diagnosis:**
- ✅ Fast vapor plume confirms strong vaporization
- Vapor is moving **8× faster than max film velocity**
- Normal: vapor expands faster than condensed material

---

### ⚡ **Energy Partitioning** (NEW)

Shows where laser energy is going:

| Energy Path | Tracked Value | Efficiency Metric |
|-------------|---------------|-------------------|
| **Laser input** | 0.0092 W (current) | Active laser power |
| **Peak volumetric source** | 5.65×10¹⁶ W/m³ | Intensity in material |
| **Metal kinetic energy** | Calculated | Energy in motion |
| **Laser→kinetic efficiency** | % | How much laser energy → motion |

**Diagnosis:**
- Low kinetic efficiency early is **normal**
- Most energy goes into: heating, vaporization, phase change
- Efficiency **increases** as film accelerates

---

## How to Diagnose Your Current State

### Question: "Am I near ejection?"

**Check these indicators:**

1. **Velocity Evolution:**
   ```
   Current:  avg = 0.2 m/s, max = 250 m/s
   Ejection: avg > 100 m/s (threshold)
   ```
   → **You are NOT at ejection yet** (avg too low)

2. **Volume Loss:**
   ```
   Current:  0.0018%
   Ejection: >0.01% and rapidly increasing
   ```
   → **Film still attached**

3. **Recoil Pressure Trend:**
   ```
   If: maxRecoil drops to <50% of peak
   AND: avgVel still high
   → SEPARATION
   ```
   → **Not yet** (recoil stable at peak)

4. **Phase:**
   ```
   Current: ACCELERATION (phase 10)
   Ejection: Will transition to TRANSFER (phase 11)
   ```
   → **Still accelerating**

### Question: "What should I watch for?"

**Critical Indicators (in order of occurrence):**

1. **Avg velocity crossing 50 m/s** → Bulk film mobilizing
   - Milestone alert will trigger

2. **Avg velocity crossing 100 m/s** → Ejection threshold
   - "⚡ MILESTONE" message appears

3. **Recoil pressure dropping** → Film separating from donor
   - Vapor layer lost confinement

4. **Volume loss spiking** → Material leaving domain
   - Transfer phase beginning

5. **Big alert box** → Separation confirmed
   ```
   🚀🚀🚀  FILM SEPARATION EVENT DETECTED!  🚀🚀🚀
   ```

---

## Time-to-Ejection Prediction

**How it works:**

```cpp
If (recoil > 10 MPa AND acceleration > 1e10 m/s²):
    time_to_ejection = (100 m/s - current_velocity) / acceleration
```

**Shows in tracker:**
```
╠═══════════════════════════════╣
║     🔮 EJECTION PREDICTION     ║
║  Est. time to ejection: X ps  ║
║  Predicted ejection vel: Y m/s║
║  Time at ejection: Z ps       ║
╚═══════════════════════════════╝
```

**Your situation:**
- Currently NOT showing → acceleration not high enough yet
- Will appear when bulk acceleration picks up

---

## Advanced Diagnostics

### Shock Wave Detection

**Indicator:**
```
Shock wave: DETECTED  (appears in tracker)
```

**Criterion:**
```
Pressure gradient > 1e15 Pa/m  (1 million GPa/m)
```

**Meaning:**
- Supersonic pressure front
- Often accompanies explosive vaporization
- Can affect film integrity

### Film Separation Algorithm

**Logical detection:**
```cpp
IF:
  1. maxRecoilSeen > 10 MPa        (strong recoil occurred)
  2. currentRecoil < 0.5 × peak    (pressure dropped 50%+)
  3. avgVel > 50 m/s               (film still moving fast)
THEN:
  → SEPARATION DETECTED
```

**Why these criteria:**
- **Recoil drop**: Vapor layer collapsed or film moved away
- **Velocity maintained**: Film has inertia (detached, not stopped)
- **Combined**: Film is ballistic (free-flying)

---

## CSV Data Analysis

**Enhanced CSV now includes 30+ parameters** for post-processing:

### Key Columns for Analysis:

1. **Time series:**
   - `time_ps`
   - `phase`, `phase_num`

2. **Thermal evolution:**
   - `Te_max_K`, `Tl_max_K`, `Tl_avg_K`, `Tl_spread_K`

3. **Pressure dynamics:**
   - `P_max_MPa`, `recoil_MPa`, `recoil_max_seen_MPa`

4. **Material tracking:**
   - `metal_vol_um3`, `metal_loss_pct`, `metal_loss_rate_um3s`

5. **Velocity analysis:**
   - `vel_max_ms`, `vel_avg_ms`, `vel_stddev_ms`, `accel_kms2`

6. **Interface topology:**
   - `interface_cells`, `interface_area_um2`

7. **Vapor plume:**
   - `vapor_vol_um3`, `vapor_cells`, `vapor_vel_avg_ms`, `vapor_vel_max_ms`

8. **Energy:**
   - `turbulent_KE_nJ`, `kinetic_eff_pct`

9. **Events:**
   - `shock_present` (0/1)
   - `separation_detected` (0/1)
   - `time_to_ejection_ps` (-1 if not applicable)

### Recommended Plots:

```python
# 1. Velocity evolution
plt.plot(time_ps, vel_avg_ms, label='Avg')
plt.plot(time_ps, vel_max_ms, label='Max')
plt.axhline(100, color='red', label='Ejection threshold')

# 2. Pressure dynamics
plt.plot(time_ps, recoil_MPa)
plt.plot(time_ps, recoil_max_seen_MPa, '--', label='Peak seen')

# 3. Material loss
plt.plot(time_ps, metal_loss_pct)
plt.yscale('log')

# 4. Vapor plume
plt.plot(time_ps, vapor_vel_max_ms)

# 5. Phase progression
plt.plot(time_ps, phase_num)
```

---

## Troubleshooting Guide

### Issue: "Progress stuck at X%"

**Diagnosis:**
1. Check which phase is current
2. Look at trigger conditions for next phase
3. Verify physics is evolving:
   - Is temperature increasing?
   - Is pressure building?
   - Is velocity changing?

**Example:**
```
Stuck at HEATING (phase 3, ~19%)?
→ Check: Is Tl reaching 1500 K (pre-melt threshold)?
→ If not: Laser power too low OR simulation time too short
```

### Issue: "Velocity too high/unrealistic"

**Check:**
1. `maxReasonableVelocity` in controlDict (default: 5000 m/s)
2. Recoil pressure magnitude
3. Time step size

**Warning triggers:**
```
⚠  WARNING: Velocity exceeds realistic LIFT range!
   Observed: XXXX m/s
   Limit: YYYY m/s
```

### Issue: "No separation detected"

**Possible reasons:**
1. **Film hasn't reached ejection velocity** (need avg > 100 m/s)
2. **Recoil still building** (not dropping yet)
3. **Not enough simulation time** (need to continue)

**Check:**
- Current avg velocity
- Recoil trend (increasing/stable/decreasing)
- Metal loss rate (should be minimal before separation)

### Issue: "Too much material loss early"

**Diagnosis:**
```
If metal_loss > 1% before ACCELERATION phase:
  → May have:
    - Too intense laser
    - Ablation instead of LIFT
    - Numerical issues
```

---

## Real-Time Monitoring Checklist

While simulation runs, watch for:

- [ ] **LASER_ACTIVE** phase transition (laser working)
- [ ] **Temperature rise** to melting/boiling (Te, Tl increasing)
- [ ] **Recoil pressure** building (>1 GPa realistic)
- [ ] **Velocity milestones**:
  - [ ] Max > 100 m/s
  - [ ] Avg > 10 m/s  ← *Ejection starts*
  - [ ] Avg > 100 m/s ← *Ejection threshold*
  - [ ] Avg > 500 m/s ← *High-velocity regime*
- [ ] **Volume loss** starts increasing
- [ ] **Separation alert** (big box with 🚀)
- [ ] **Phase transitions** through 11-16 (post-ejection)

---

## Interpreting Your Current Status

### Your Simulation (at 0.76 ps):

| Metric | Value | Status |
|--------|-------|--------|
| **Phase** | ACCELERATION (10/16) | 62.5% ✅ |
| **Time** | 0.76 ps | Early stage |
| **Laser** | OFF (pulse ended at 0.2 ps) | ℹ️ |
| **Temperature** | Te=10.6 kK, Tl=10 kK | Very hot ✅ |
| **Recoil** | 37.4 GPa | Extremely high ✅ |
| **Avg velocity** | 0.195 m/s | Still low ⏳ |
| **Max velocity** | 250 m/s | Surface accelerating ✅ |
| **Volume loss** | 0.0018% | Negligible ✅ |
| **Vapor plume** | 2047 m/s max | Supersonic ✅ |

### Physical Interpretation:

**What's happening:**
1. Laser pulse finished (0.2 ps) but energy still redistributing
2. Extreme temperatures creating **massive vapor pressure**
3. **Vapor layer** between film and donor with 37 GPa recoil
4. Film **surface** accelerating (250 m/s) but bulk **lagging** (0.2 m/s)
5. Vapor plume expanding supersonically (2 km/s)

**What comes next (predictions):**
1. **0.8-1.5 ps**: Bulk velocity rapidly increases as film responds to recoil
2. **1-3 ps**: Avg velocity → 100+ m/s, **ejection event**
3. **1-5 ps**: Recoil pressure drops, **separation detected**
4. **5-50 ps**: Ballistic flight, cooling, transfer phase
5. **50-500 ps**: Film reaches receiver, solidifies

**You are:** ✅ **Progressing normally toward ejection**

---

## Summary: What Determines Progress

### Primary Factors:

1. **Temperature Evolution**
   - Electron heating → Lattice heating → Melting → Vaporization
   - Drives phase transitions 1-6

2. **Pressure Build-up**
   - Thermal + vapor pressure → Recoil
   - Enables phases 7-8

3. **Velocity Development**
   - Recoil acceleration → Film motion → Ejection
   - Defines phases 9-11

4. **Time After Laser**
   - Energy redistribution → Cooling → Solidification
   - Governs phases 12-16

### Secondary Factors:

- Material properties (Fe melting/boiling points)
- Laser parameters (power, duration, spot size)
- Domain size and boundary conditions
- Time step and solver convergence

---

## Quick Reference: "Where Am I?"

**Use this decision tree:**

```
Is laser ON? → Phase 1-2
Is Tl < 1500K? → Phase 3-4
Is Tl < 3560K? → Phase 5
Is Tl < 5000K? → Phase 6
Is recoil < 10 MPa? → Phase 7
Is velocity < 10 m/s? → Phase 8
Is velocity < 100 m/s? → Phase 9
Is velocity < high regime? → Phase 10 ← YOU ARE HERE
Is metal loss < 2%? → Stay in 10
Is Tl still high? → Phase 11
Is Tl decreasing? → Phase 12-16
```

**Your Position:**
- ✅ Tl = 10,000 K (very high)
- ✅ Recoil = 37,400 MPa (>>10 MPa)
- ✅ Velocity max = 250 m/s (>100 m/s)
- ⏳ Velocity avg = 0.2 m/s (<100 m/s) ← **Rate-limiting step**
- ⏳ Metal loss = 0.002% (<2%)

**Next milestone:** Bulk velocity reaching 100 m/s

---

## File Output

**File:** `liftProcessTracking.csv` (updated every 25 timesteps)

**Location:** Case directory (same as log file)

**Usage:**
```bash
# Plot progress
python3 plot_lift_progress.py liftProcessTracking.csv

# Analyze specific column
awk -F',' '{print $1, $16}' liftProcessTracking.csv  # time vs vel_avg
```

---

## Conclusion

You now have **comprehensive real-time diagnostics** with:

- ✅ **16-phase tracking** with clear criteria
- ✅ **8 diagnostic categories** covering all physics
- ✅ **Event notifications** for milestones
- ✅ **Separation detection** algorithm
- ✅ **Ejection prediction** (when conditions met)
- ✅ **30+ CSV columns** for post-processing
- ✅ **Visual progress bar** and formatted output
- ✅ **Shock wave detection**
- ✅ **Vapor plume tracking**
- ✅ **Energy efficiency** analysis

**Your simulation is progressing normally** - continue running to see ejection!

---

*LIFT Process Tracker v2.0 - Enhanced Physics Diagnostics*
