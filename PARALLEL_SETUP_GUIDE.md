# Parallel Processing Setup Guide for compInterFoam

This guide provides instructions for setting up and running parallel simulations with `compInterFoam` on your system with an Intel Core i7-9750H (12 threads).

## System Requirements

Based on your system information:
- **CPU**: Intel Core i7-9750H (6 cores, 12 threads @ 2.60GHz)
- **RAM**: 32 GB
- **GPU**: NVIDIA GeForce RTX 2060
- **OS**: Ubuntu 24.04.2 LTS
- **OpenFOAM**: v2406

## Quick Start

### 1. Generate System Information Report

First, generate a comprehensive system report to verify your setup:

```bash
cd /home/user/compInterFoam
./generate_system_report.sh
```

This will create a `system_info_report.txt` file with detailed information about your system configuration.

### 2. Prepare Your Case

Before running parallel simulations, ensure your case is properly set up:

```bash
# Navigate to your case directory (e.g., TEST1 or RealisticLIFT)
cd TEST1

# Make sure OpenFOAM environment is loaded
source /home/mavadi/OpenFOAM-v2406/etc/bashrc

# Generate the mesh (if not already done)
blockMesh

# Verify mesh quality
checkMesh
```

### 3. Decompose the Domain

The `decomposeParDict` files have been configured for your 12-core system:

```bash
# Decompose the mesh into 12 subdomains
decomposePar

# Optional: Check the decomposition quality
decomposePar -cellDist
```

This will create `processor0` through `processor11` directories containing the decomposed mesh and fields.

### 4. Run the Simulation in Parallel

Run your simulation using MPI:

```bash
# Run on 12 cores
mpirun -np 12 compInterFoam -parallel > log.compInterFoam 2>&1 &

# Monitor progress
tail -f log.compInterFoam
```

### 5. Reconstruct Results

After the simulation completes, reconstruct the fields:

```bash
# Reconstruct all time directories
reconstructPar

# Reconstruct only the latest time
reconstructPar -latestTime

# Reconstruct specific times
reconstructPar -time '0.1:0.5'
```

### 6. Cleanup Processor Directories (Optional)

After reconstruction, you can remove the processor directories to save space:

```bash
# Remove processor directories
rm -rf processor*

# Or use the Allclean script if available
./Allclean
```

## DecomposeParDict Configuration

The `decomposeParDict` files in both `TEST1/system/` and `RealisticLIFT/system/` are configured with:

### Current Settings (Scotch Method)
```c++
numberOfSubdomains 12;
method scotch;
```

**Scotch** is an automatic load-balancing method that works well for complex geometries. It minimizes inter-processor boundaries and ensures balanced workload distribution.

### Alternative Methods

#### Hierarchical Decomposition
If you prefer geometric decomposition, you can modify the `decomposeParDict`:

```c++
method hierarchical;

hierarchicalCoeffs
{
    n           (3 2 2);  // 3*2*2 = 12 subdomains
    delta       0.001;
    order       xyz;
}
```

This divides the domain into 3 parts along x, 2 along y, and 2 along z.

#### Simple Decomposition
```c++
method simple;

simpleCoeffs
{
    n           (3 2 2);
    delta       0.001;
}
```

## Performance Optimization

### CPU Affinity
For better performance, you can pin MPI processes to specific cores:

```bash
mpirun --bind-to core --map-by core -np 12 compInterFoam -parallel
```

### Memory Considerations
With 32 GB RAM, you have approximately 2.67 GB per core. Monitor memory usage:

```bash
# During simulation, in another terminal:
watch -n 5 'free -h'
```

### Choosing the Number of Processors

While your system has 12 threads, consider these factors:

1. **Mesh Size**: Small meshes may not benefit from 12 processors due to communication overhead
2. **Memory**: Ensure each processor has enough memory (at least 1-2 GB)
3. **Efficiency**: Test scalability with different core counts (e.g., 4, 6, 8, 12)

Example scaling test:
```bash
# Test with different core counts
for np in 4 6 8 12; do
    decomposePar -numberOfSubdomains $np
    mpirun -np $np compInterFoam -parallel > log.np${np} 2>&1
    reconstructPar
    rm -rf processor*
done
```

## Troubleshooting

### Issue: decomposePar fails with "inconsistent patch"
**Solution**: Check your mesh with `checkMesh` and ensure all boundaries are properly defined

### Issue: Simulation crashes during parallel run
**Possible causes**:
- Insufficient memory per processor
- Unstable numerical schemes
- CFL number too high

**Solution**:
- Reduce number of processors
- Check `system/fvSchemes` and `system/fvSolution`
- Reduce time step or enable `adjustTimeStep`

### Issue: Poor parallel speedup
**Possible causes**:
- Mesh too small relative to number of processors
- Poor decomposition balance
- High communication overhead

**Solution**:
- Use fewer processors
- Try different decomposition methods
- Check load balance: `decomposePar -cellDist`

## Monitoring Parallel Performance

### Check Load Balance
```bash
# After decomposition
foamJob -l decomposePar -cellDist
```

### Monitor Resource Usage
```bash
# Install and use htop for real-time monitoring
htop

# Or use top with specific process filtering
top -p $(pgrep -d',' compInterFoam)
```

## Best Practices

1. **Always verify mesh quality** before decomposition using `checkMesh`
2. **Test with smaller problems** before running large simulations
3. **Save processor directories** until you verify results are correct
4. **Use version control** for your case setups
5. **Document your runs** in a log book or README
6. **Backup important results** regularly

## Case-Specific Notes

### TEST1 Case
- Configured with standard settings
- Suitable for testing and validation

### RealisticLIFT Case
- Production case with checkpoint preservation
- Uses `purgeWrite` disabled to keep all time directories
- Ensure sufficient disk space before running

## Additional Resources

- OpenFOAM User Guide: https://www.openfoam.com/documentation/user-guide
- Parallel Computing Guide: https://www.openfoam.com/documentation/guides/latest/doc/guide-parallel.html
- OpenFOAM Forums: https://www.cfd-online.com/Forums/openfoam/

## Quick Reference Commands

```bash
# Full workflow
source $FOAM_BASHRC
blockMesh
checkMesh
decomposePar
mpirun -np 12 compInterFoam -parallel > log &
tail -f log
reconstructPar
paraFoam

# Cleanup
rm -rf processor*
foamListTimes -rm
```

## Performance Benchmarks

Consider running benchmarks to understand your system's performance:

```bash
# Simple scaling test
./generate_system_report.sh
# Document results in a spreadsheet for future reference
```

Record:
- Number of processors used
- Wall clock time
- Speedup factor
- Parallel efficiency (%)

---

**Last Updated**: 2025-11-11
**OpenFOAM Version**: v2406
**Target System**: Intel i7-9750H with 12 threads, 32GB RAM
