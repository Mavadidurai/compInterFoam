

# SMART LIFT Simulation for PoliMi HPC

## The Pragmatic Approach

This is the **realistic** ultra-fine LIFT setup for HPC clusters. Unlike the "nano-everywhere" approach (102M cells = suicide), this puts high resolution **where it matters** and saves computational resources everywhere else.

---

## Key Principle: Smart Grading

**Don't waste cells on bulk regions far from the action.**

- **Film region**: 2 nm cells (captures nano-scale physics)
- **Near-film zone**: 10-20 nm cells (graded approach)
- **Bulk regions**: 80-100 nm cells (adequate for thermal bulk)

Result: **~5 Million cells** vs 102 Million, but **same physics quality in critical regions**.

---

## Mesh Specifications

### Total: **4.86 Million cells**

| Region | Y Range | Cells | Grading | Cell Size |
|--------|---------|-------|---------|-----------|
| **Substrate** | 0 → 8 μm | 200 | 0.2 | 80nm (bottom) → **16nm (top)** |
| **Air gap** | 8 → 20 μm | 250 | 0.1 | 96nm (bottom) → **10nm (top)** |
| **Ti film** | 20 → 20.0714 μm | 36 | 1.0 | **2 nm (uniform)** |

**Lateral**: 400 cells × 30 μm = 75 nm/cell
**Depth**: 25 cells × 3 μm = 120 nm/cell

### What This Achieves

✅ **2 nm resolution** in Ti film (perfect for e-ph coupling)
✅ **10 nm resolution** at film interfaces (clean phase change)
✅ **Smooth grading** to bulk (no mesh artifacts)
✅ **75 nm lateral** (captures 6 μm laser spot with 80 cells)

---

## Critical Design Choice: ny=1 Decomposition

### The Problem with Splitting the Film

The Ti film is only **71.4 nm = 36 cells thick**.

**Bad decomposition** (ny > 1):
- Processor boundaries cut through film
- Each processor gets ~9 cells of film
- Heat transfer across boundaries is slow
- Phase change at boundaries is messy
- Interface tracking breaks

**Smart decomposition** (ny=1):
- Each processor owns **full vertical stack**
- Film stays intact per processor
- Clean physics
- No vertical communication for film

### Decomposition Strategy

```
128 cores: (32, 1, 4) = 32 in X × 1 in Y × 4 in Z
```

**Each subdomain:**
- ~13 cells wide (X)
- **486 cells tall (full height)**
- ~6 cells deep (Z)
- **~38,000 cells/core** → excellent balance

**Why this works:**
- Film is small (36 cells) compared to total Y (486 cells)
- Mesh is graded → bulk regions are cheap
- X and Z splits provide good parallelism
- No communication overhead for thin-layer physics

---

## Physics Coverage

### Complete LIFT Cycle

**Phase 1** (0-0.2 ps): Laser Heating
- Electron temp → 20,000 K
- Lattice temp → 6,000 K
- Two-temperature dynamics

**Phase 2** (0.2-5 ps): Melting & Superheating
- T > 1,941 K (melting)
- T > 2,200 K (evaporation onset)
- Recoil pressure builds

**Phase 3** (5-50 ps): Vaporization & Jet
- Material ejection
- Jet/droplet formation
- Transfer initiation

**Phase 4** (50+ ps): Cooling & Solidification
- (Extend endTime as needed)

### Time Control: Adaptive but Realistic

```
Initial dt:  0.1 fs  (for fs-laser)
Max dt:      1 ps    (aggressive to reach 50 ps)
Target:      50 ps in 1-2 days runtime
```

**Courant limits:**
- Flow: 0.5
- Interface: 0.3
- Thermal: 0.3 (critical for 2nm cells!)

---

## Computational Requirements

### Recommended Configuration

| Cores | Nodes | Memory | Cells/Core | Est. Time/Step | To 50 ps |
|-------|-------|--------|------------|----------------|----------|
| 64    | 2     | 192 GB | 76k        | ~0.3 s         | ~2-3 days |
| **128**   | **4**     | **384 GB** | **38k**        | **~0.15 s**    | **~1-2 days** |
| 256   | 8     | 768 GB | 19k        | ~0.08 s        | ~12-24 hrs |

**Storage:**
- Mesh (decomposed): ~10 GB
- Output (50 times): ~50-100 GB
- Total: ~100-150 GB

---

## Quick Start

### 1. Transfer to Cluster

```bash
scp -r TEST_LIFT_POLIMI_SMART username@polimi-hpc.it:~/
ssh username@polimi-hpc.it
cd TEST_LIFT_POLIMI_SMART
```

### 2. Edit SLURM Script

```bash
nano run_lift_smart.slurm
```

**Change:**
- Line 9: `#SBATCH --account=your_account_name`
- Line 11: `#SBATCH --partition=...` (check your cluster)
- Line 29: OpenFOAM path if different

### 3. Run Setup

```bash
chmod +x Allrun
./Allrun
```

This will:
- Generate mesh (~1 min)
- Check quality
- Copy initial conditions
- Decompose for 128 cores (~5 min)

### 4. Submit Job

```bash
sbatch run_lift_smart.slurm
```

### 5. Monitor

```bash
# Check queue
squeue -u $USER

# Monitor output
tail -f slurm_*.out

# Check progress
tail -f log.compInterFoam
grep "^Time = " log.compInterFoam | tail -5

# Check Courant
grep "Courant" log.compInterFoam | tail -10
```

---

## Adjusting Core Count

### For 64 Cores

**decomposeParDict:**
```cpp
numberOfSubdomains 64;
hierarchicalCoeffs
{
    n    (16 1 4);   // 16×1×4 = 64
}
```

**SLURM script:**
```bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
NPROCS=64
```

### For 256 Cores

**decomposeParDict:**
```cpp
numberOfSubdomains 256;
hierarchicalCoeffs
{
    n    (64 1 4);   // 64×1×4 = 256
}
```

**SLURM script:**
```bash
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=32
NPROCS=256
```

Then re-decompose:
```bash
rm -rf processor*
decomposePar -copyZero
```

---

## Extending Simulation Time

To run longer (e.g., 100 ps for solidification):

**Edit system/controlDict:**
```cpp
endTime         1e-10;         // 100 ps
maxDeltaT       2e-12;         // Allow 2 ps steps
writeInterval   2e-12;         // Write every 2 ps
```

**Edit SLURM:**
```bash
#SBATCH --time=72:00:00        # 3 days
```

---

## Troubleshooting

### Simulation Crashes Early

**Check:**
```bash
tail -100 log.compInterFoam | grep -i error
```

**Common fixes:**
- Reduce `maxCo` to 0.3 in controlDict
- Reduce `maxDeltaT` to 1e-13
- Increase `maxElectronSubCycles` to 30

### Timestep Too Small

```bash
grep "deltaT" log.compInterFoam | tail -20
```

**If stuck at minDeltaT:**
- Check temperature extremes (may be unphysical)
- Check mesh quality near film
- Increase `minDeltaT` to 1e-16

### Out of Memory

```bash
sacct -j JOBID --format=JobID,MaxRSS
```

**If exceeding memory:**
- Reduce `purgeWrite` to 25
- Increase `writeInterval`
- Request more memory per node

---

## Expected Results

### Key Observables

**Early phase** (0-1 ps):
- Te spikes to 10-20 kK
- Tl reaches 2-6 kK
- Melting begins

**Mid phase** (1-10 ps):
- Superheating above 2200 K
- Recoil pressure 0.1-3 GPa
- Film bulging

**Late phase** (10-50 ps):
- Material ejection
- Jet formation
- Alpha.metal < 1 in film region

### Data Analysis

```bash
# Reconstruct key times
reconstructPar -time '0.2e-12,1e-12,5e-12,10e-12,50e-12'

# View in ParaView
paraFoam

# Extract probe data
cd postProcessing/probes/
gnuplot plot_probes.gp  # (create your own)
```

---

## File Structure

```
TEST_LIFT_POLIMI_SMART/
├── 0.orig/               Initial conditions
├── constant/             Material properties
│   ├── laserProperties
│   ├── thermophysicalProperties.*
│   └── ...
├── system/
│   ├── blockMeshDict     GRADED mesh with smart resolution
│   ├── controlDict       Adaptive time, realistic targets
│   ├── decomposeParDict  ny=1 (film-preserving decomp)
│   ├── fvSchemes
│   └── fvSolution
├── Allrun                Setup script
└── run_lift_smart.slurm  SLURM submission
```

---

## Comparison: Smart vs Ultra-Fine

| Feature | Ultra-Fine Everywhere | SMART (This Case) |
|---------|----------------------|-------------------|
| **Cells** | 102 Million | 5 Million |
| **Film resolution** | 1 nm | 2 nm |
| **Bulk resolution** | 5 nm (wasted!) | 80 nm (adequate) |
| **Cores needed** | 256-512 | 64-256 |
| **Runtime to 50 ps** | Weeks | 1-2 days |
| **Storage** | ~2 TB | ~100 GB |
| **Practical?** | ❌ Suicide | ✅ Realistic |

---

## Why This Works

1. **Physics is local**: Nano-scale effects matter only in/near film
2. **Grading is smooth**: No mesh artifacts at transitions
3. **Decomposition respects geometry**: Film stays intact
4. **Time stepping is adaptive**: Fast where possible
5. **Output is manageable**: Purging keeps only recent history

**This is how you run LIFT on a cluster without suicide.**

---

## Next Steps

1. ✅ Setup complete
2. ✅ Ready to submit
3. → Run to 50 ps
4. → Analyze jet formation
5. → Extend to 100-200 ps if needed
6. → Publish!

---

## Support

**Documentation:**
- This README (complete guide)
- `system/blockMeshDict` (detailed mesh comments)
- `system/controlDict` (runtime estimates)
- `system/decomposeParDict` (decomposition strategy)

**Common Commands:**
```bash
# Setup
./Allrun

# Submit
sbatch run_lift_smart.slurm

# Monitor
tail -f log.compInterFoam

# Check time
grep "^Time = " log.compInterFoam | tail -1

# Check performance
grep "ExecutionTime" log.compInterFoam | tail -10

# Post-process
reconstructPar -latestTime
paraFoam
```

**Need help?**
- Check log files first
- Verify mesh: `checkMesh`
- Test small: reduce endTime to 1e-12 for quick test
- Contact cluster support for SLURM/module issues

---

**Ready to run realistic LIFT physics? Submit the job!** 🚀
