# 3-Shot Bracket Search for Experimental Validation

## Objective
Find pulse energy that produces **6-8 kK peak Tl** with **sustained recoil pressure 5-80 MPa**.

## Known Bounds
- **18.6 nJ** → 1.08 kK (too low, no vaporization)
- **30.0 nJ** → 17.0 kK (too high, unstable/clamped)
- **Transition zone**: 20-26 nJ

## 3-Shot Bracket Plan

### Run 1: 20.0 nJ (lower bracket)
```bash
# Edit constant/laserProperties line 57:
pulseEnergy              2.00e-8;   // 20.0 nJ

# Run
rm -rf 0 [0-9]* postProcessing *.log
cp -r 0.orig 0 && setFields
compInterFoam 2>&1 | tee run_20nJ.log

# Extract key metrics:
grep "max(Tl):" run_20nJ.log | tail -20 > metrics_20nJ.txt
grep "Max |recoilPressure|" run_20nJ.log | tail -20 >> metrics_20nJ.txt
```

### Run 2: 22.5 nJ (center bracket) ← CURRENT
```bash
# Already set in laserProperties
pulseEnergy              2.25e-8;   // 22.5 nJ

# Run
rm -rf 0 [0-9]* postProcessing *.log
cp -r 0.orig 0 && setFields
compInterFoam 2>&1 | tee run_22.5nJ.log

# Extract metrics:
grep "max(Tl):" run_22.5nJ.log | tail -20 > metrics_22.5nJ.txt
grep "Max |recoilPressure|" run_22.5nJ.log | tail -20 >> metrics_22.5nJ.txt
```

### Run 3: 25.0 nJ (upper bracket)
```bash
# Edit constant/laserProperties line 57:
pulseEnergy              2.50e-8;   // 25.0 nJ

# Run
rm -rf 0 [0-9]* postProcessing *.log
cp -r 0.orig 0 && setFields
compInterFoam 2>&1 | tee run_25nJ.log

# Extract metrics:
grep "max(Tl):" run_25nJ.log | tail -20 > metrics_25nJ.txt
grep "Max |recoilPressure|" run_25nJ.log | tail -20 >> metrics_25nJ.txt
```

---

## Decision Tree (After 3 Shots)

### Case A: All 3 below target (max Tl < 5 kK)
→ **Step up to 26-28 nJ bracket**, repeat

### Case B: Straddle found (e.g., 22.5 too low, 25.0 good)
→ **Bisect**: Try (22.5 + 25.0)/2 = 23.75 nJ

### Case C: All 3 above target (max Tl > 9 kK or clamps)
→ **Step down to 18-21 nJ bracket**, repeat

### Case D: Direct hit (one run in 6-8 kK window)
→ **Done!** Use that energy value

---

## Acceptance Criteria

**✓ SUCCESS:**
- Peak Tl: 6,000 - 8,000 K
- Sustained (≥0.5 ps) above 6 kK
- Recoil pressure: 5 - 80 MPa (sustained)
- No pressure/velocity clamp hits
- Co ≤ 0.5, thermal-Co stable

**✗ TOO LOW:**
- Peak Tl < 3,560 K (never reaches Tvap)
- Recoil pressure = 0 MPa
- Phase: MELTING/SOLIDIFYING (no ejection)

**✗ TOO HIGH:**
- Peak Tl > 10,000 K
- Hitting maxPressure (200 MPa) clamp
- Velocity warnings (> 800 m/s)
- Solver instabilities

---

## Key Diagnostics to Extract

From each run's log, get:

1. **Peak temperatures** (last 20 timesteps):
   ```
   max(Te): XXXX K
   max(Tl): XXXX K
   ```

2. **Recoil pressure evolution**:
   ```
   Max |recoilPressure| = XX MPa
   ```

3. **Phase transitions**:
   ```
   Phase: ABSORPTION / MELTING / VAPORIZING / EJECTING / SOLIDIFYING
   ```

4. **Clamp messages** (any present = too high):
   ```
   pressureClamp: enforcing bounds [0, 2e+08] Pa
   maxTl clamp triggered
   ```

5. **Energy audit** (verify consistent delivery):
   ```
   Laser power: XXXX W
   Metal-absorbed power: XXXX W
   ```

---

## Safety Rails

- **Courant limits**: maxCo = 0.5, maxAlphaCo = 0.5 (already set)
- **Thermal Courant**: maxThermalCourant = 0.5 (already set)
- **Adaptive timestep**: enabled, minDeltaT = 1e-14, maxDeltaT = 1e-11
- **Pressure clamp**: 0 - 200 MPa (already set)
- **Temperature clamp**: maxTl = 7000 K (already set)

If instabilities occur:
1. Reduce maxCo to 0.2
2. Reduce maxThermalCourant to 0.1
3. Halve maxDeltaT to 5e-12
4. Rerun that energy point

---

## Expected Timeline

- Each run: ~10-12 hours
- 3 shots: ~30-36 hours total
- 1 bisection: +10-12 hours
- **Total convergence**: ~40-48 hours (2 days)

---

## Notes

- Keep **all other parameters fixed** (mesh, BCs, model coefficients)
- Only change `pulseEnergy` between runs
- Document results in this file after each run
- Compare energy audits to verify consistent absorption
