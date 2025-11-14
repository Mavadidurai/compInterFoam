#!/bin/bash
#SBATCH --job-name=LIFT_smart
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=384G
#SBATCH --partition=parallel
#SBATCH --account=your_account_name

#==============================================================================
# SMART LIFT SIMULATION - PoliMi HPC
#==============================================================================
#
# Configuration: 4 nodes × 32 cores = 128 cores
# Mesh: ~5M cells → ~38k cells/core (excellent balance)
# Memory: ~3 GB/core
# Target: 50 ps in ~1-2 days
#
# BEFORE SUBMITTING:
# 1. Update --account=your_account_name
# 2. Check --partition (may vary by cluster)
# 3. Verify OpenFOAM path below
#==============================================================================

NPROCS=128
CASE_DIR="${SLURM_SUBMIT_DIR}"
FOAM_PATH="/opt/openfoam2406/etc/bashrc"

echo "============================================================================"
echo " SMART LIFT SIMULATION - Job ${SLURM_JOB_ID}"
echo "============================================================================"
echo "Started:   $(date)"
echo "Nodes:     ${SLURM_JOB_NODELIST}"
echo "Cores:     ${SLURM_NTASKS}"
echo "Memory:    ${SLURM_MEM_PER_NODE} MB/node"
echo "============================================================================"
echo ""

# Load modules (ADJUST FOR YOUR CLUSTER)
module purge
module load gcc/11.2.0
module load openmpi/4.1.1

# Source OpenFOAM
if [ -f "${FOAM_PATH}" ]; then
    source "${FOAM_PATH}"
    echo "OpenFOAM loaded: ${WM_PROJECT_DIR}"
else
    echo "ERROR: OpenFOAM not found at ${FOAM_PATH}"
    exit 1
fi

cd "${CASE_DIR}" || exit 1
echo "Case: ${CASE_DIR}"
echo ""

# Pre-flight checks
if [ ! -d "processor0" ]; then
    echo "ERROR: Case not decomposed!"
    exit 1
fi

NPROC_DIRS=$(ls -d processor* | wc -l)
if [ "${NPROC_DIRS}" -ne "${NPROCS}" ]; then
    echo "ERROR: Found ${NPROC_DIRS} processors, expected ${NPROCS}"
    exit 1
fi

echo "✓ Pre-flight checks passed"
echo ""

# Performance settings
export FOAM_SIGFPE=false
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=close

echo "============================================================================"
echo " RUNNING SIMULATION"
echo "============================================================================"
echo "Start: $(date)"
echo ""

# Run with timing
time mpirun -np ${NPROCS} compInterFoam -parallel > log.compInterFoam 2>&1

EXIT_STATUS=$?

echo ""
echo "============================================================================"
echo " SIMULATION COMPLETE"
echo "============================================================================"
echo "End: $(date)"
echo "Exit status: ${EXIT_STATUS}"
echo ""

if [ ${EXIT_STATUS} -eq 0 ]; then
    echo "✓ SUCCESS!"

    if [ -f "log.compInterFoam" ]; then
        FINAL_TIME=$(grep "^Time = " log.compInterFoam | tail -1 | awk '{print $3}')
        TIMESTEPS=$(grep "^Time = " log.compInterFoam | wc -l)
        echo "Final time: ${FINAL_TIME} s"
        echo "Timesteps:  ${TIMESTEPS}"

        # Quick diagnostics
        echo ""
        echo "Final Courant numbers:"
        grep "Courant Number" log.compInterFoam | tail -3

        echo ""
        echo "Temperature extremes:"
        grep "max(Tl)" log.compInterFoam | tail -3
    fi
else
    echo "✗ FAILED (exit ${EXIT_STATUS})"
    echo "Check log.compInterFoam for errors"
fi

echo ""
echo "Post-processing:"
echo "  reconstructPar -latestTime"
echo "  paraFoam"
echo ""
echo "============================================================================"

# EOF
