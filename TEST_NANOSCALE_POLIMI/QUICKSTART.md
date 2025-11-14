# Quick Start Guide - Nanoscale LIFT on PoliMi HPC

## 30-Second Setup

```bash
# 1. Make scripts executable
chmod +x Allrun Allclean Allpost run_nanoscale_lift.slurm

# 2. Edit SLURM script (add your account name)
nano run_nanoscale_lift.slurm
# Change: #SBATCH --account=your_account_name

# 3. Generate mesh and decompose
blockMesh | tee log.blockMesh
cp -r 0.orig 0
decomposePar -copyZero | tee log.decomposePar

# 4. Submit job
sbatch run_nanoscale_lift.slurm

# 5. Monitor
tail -f slurm_*.out
```

## Essential Commands

### Before Running
```bash
# Check mesh
checkMesh -allTopology | tee log.checkMesh

# Verify decomposition
ls -d processor* | wc -l  # Should match numberOfSubdomains
```

### During Run
```bash
# Check job status
squeue -u $USER

# Monitor progress
tail -f log.compInterFoam
grep "^Time = " log.compInterFoam | tail -5

# Check Courant numbers
grep "Courant" log.compInterFoam | tail -10
```

### After Run
```bash
# Reconstruct latest time
reconstructPar -latestTime

# View in ParaView
paraFoam
```

## Key Files to Modify

1. **run_nanoscale_lift.slurm**
   - Line 9: `#SBATCH --account=YOUR_ACCOUNT`
   - Line 6-8: Adjust nodes/cores if needed
   - Line 41: Update OpenFOAM path

2. **system/decomposeParDict**
   - Line 17: `numberOfSubdomains` (must match SLURM cores)

3. **system/controlDict**
   - Line 23: `endTime` (extend/reduce simulation time)
   - Line 28: `writeInterval` (adjust output frequency)

## Common Configurations

### 64 cores (small test)
```bash
# In decomposeParDict:
numberOfSubdomains 64;
hierarchicalCoeffs { n (8 4 2); }

# In SLURM script:
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
NPROCS=64
```

### 128 cores (recommended)
```bash
# Already configured by default!
```

### 256 cores (production)
```bash
# In decomposeParDict:
numberOfSubdomains 256;
hierarchicalCoeffs { n (16 8 2); }

# In SLURM script:
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=32
NPROCS=256
```

## Troubleshooting One-Liners

```bash
# Simulation crashed? Check last error
tail -100 log.compInterFoam | grep -i error

# Timestep too small?
grep "deltaT" log.compInterFoam | tail -20

# Out of memory?
sacct -j JOBID --format=JobID,MaxRSS,Elapsed

# Restart from latest time
# (Edit controlDict: startFrom latestTime)
# Then resubmit

# Clean everything and start over
./Allclean
```

## Expected Output

```
Time directories: 0, 2e-13, 4e-13, ..., 2e-9
Log file: log.compInterFoam
SLURM output: slurm_JOBID.out
Processor dirs: processor0/ ... processor127/
```

## Performance Check

```bash
# Calculate average timestep duration
STEPS=$(grep "^Time = " log.compInterFoam | wc -l)
TIME=$(grep "ExecutionTime" log.compInterFoam | tail -1 | awk '{print $3}')
echo "scale=4; $TIME / $STEPS" | bc
# Should be 0.1-0.5 seconds per step for 128 cores
```

## Critical Parameters

| What | Where | Default | Notes |
|------|-------|---------|-------|
| Cells | blockMeshDict | 102M | Don't change |
| Cores | decomposeParDict | 128 | Match SLURM |
| End time | controlDict | 2e-9 s | 2 nanoseconds |
| Write interval | controlDict | 2e-13 s | 200 femtoseconds |
| Laser energy | laserProperties | 60 nJ | LIFT threshold |

## Help

- Full docs: See README.md
- Clean case: `./Allclean`
- Post-process: `./Allpost`
- Check mesh: `checkMesh`
- Solver help: `compInterFoam -help`

---
**Ready to run? → `sbatch run_nanoscale_lift.slurm`**
