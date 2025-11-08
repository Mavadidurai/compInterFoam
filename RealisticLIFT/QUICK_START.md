# RealisticLIFT Quick Start Guide

## 🚨 Your Simulation Issue

Your simulation was running **870 days to complete**! Here's why and how to fix it:

**Problem**: Time step too small (0.01 fs) + Serial execution + Excessive output

**Solution**: Apply optimizations → **1-2 days to complete** ✓

---

## ⚡ Quick Fix (3 Simple Steps)

### Step 1: Apply Optimizations

```bash
cd ~/OpenFOAM/mavadi-v2406/run/RealisticLIFT
./applyOptimizations.sh
```

This automatically:
- ✓ Increases time step 5-100x
- ✓ Reduces output frequency 100x
- ✓ Backs up your original settings

### Step 2: Clean & Run

```bash
# Clean old data
./Allclean

# Run simulation
compInterFoam > log.compInterFoam 2>&1 &
```

### Step 3: Monitor Progress

```bash
# Watch log in real-time
tail -f log.compInterFoam

# Check temperatures (in another terminal)
watch -n 60 'grep "max(Te):" log.compInterFoam | tail -1; grep "max(Tl):" log.compInterFoam | tail -1'
```

---

## 📊 What to Expect

After 1-2 hours, you should see:

```
Time = 1e-10 s (100 picoseconds)          ✓ Good progress
max(Te): 8000-12000 K                      ✓ Electrons hot
max(Tl): 2000-4000 K                       ✓ Lattice heating
Max |recoilPressure| = 10-80 MPa          ✓ Phase change active!
```

---

## 🎯 Success Indicators

| Metric | Target | How to Check |
|--------|--------|--------------|
| Progress | >100 ps in 2 hours | `grep "^Time = " log.compInterFoam \| tail -1` |
| Electron Temp | >5000 K | `grep "max(Te):" log.compInterFoam \| tail -1` |
| Lattice Temp | >2000 K | `grep "max(Tl):" log.compInterFoam \| tail -1` |
| Recoil Pressure | >0 MPa | `grep "Max \|recoilPressure\|" log.compInterFoam \| tail -1` |
| No errors | No "FATAL" | `grep -i "fatal\|error" log.compInterFoam` |

---

## 🚀 For Even Faster (Optional)

### Run in Parallel (8 cores → 6-8x speedup)

```bash
# 1. Create decomposition dictionary
cat > system/decomposeParDict << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}
numberOfSubdomains 8;
method          scotch;
EOF

# 2. Decompose mesh
decomposePar

# 3. Run parallel
mpirun -np 8 compInterFoam -parallel > log.compInterFoam 2>&1 &

# 4. Monitor
tail -f log.compInterFoam
```

---

## 📚 Full Documentation

- **OPTIMIZATION_RECOMMENDATIONS.md** - Detailed technical guide
- **PRE_SIMULATION_REPORT.md** - Original validation report
- **applyOptimizations.sh** - Auto-optimization script
- **diagnoseSimulation.sh** - Performance diagnostic tool

---

## ⚠️ Troubleshooting

### Simulation diverges (nan, FOAM FATAL ERROR)
```bash
# Restore original (more conservative) settings
cp system/controlDict.original system/controlDict
./Allclean
compInterFoam > log.compInterFoam 2>&1 &
```

### Still too slow
```bash
# Further reduce output
# Edit system/controlDict, change:
writeInterval   1e-11;  # Write every 10 ps instead of 1 ps
```

### Temperature not rising
```bash
# Check laser is working
grep "Laser input:" log.compInterFoam | tail -5
# Should show ~5000-9000 W

# If zero, check laser focus and timing in constant/laserProperties
```

---

## 🆘 Need Help?

Run diagnostics:
```bash
./diagnoseSimulation.sh
```

This will tell you:
- Current progress and ETA
- Temperature evolution
- Performance bottlenecks
- Specific recommendations

---

**Ready? Run `./applyOptimizations.sh` and start your simulation!** 🚀
