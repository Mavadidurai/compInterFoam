#!/bin/bash
#SBATCH --job-name=LIFT_nanoscale
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --partition=parallel
#SBATCH --mem=512G
#SBATCH --account=your_account_name

#==============================================================================
# SLURM Job Script for Nanoscale LIFT Simulation - PoliMi HPC
#==============================================================================
#
# Job configuration:
# - 4 nodes × 32 cores = 128 total cores
# - ~512 GB total memory (~4 GB per core)
# - 72 hours walltime
# - Binary parallel execution
#
# MODIFY BEFORE SUBMISSION:
# 1. --account=your_account_name (set your PoliMi HPC account)
# 2. --partition=parallel (check available partitions)
# 3. OpenFOAM environment path (see below)
# 4. Adjust --nodes, --ntasks-per-node for different core counts
#
#==============================================================================

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------

# Total number of MPI processes (must match decomposeParDict)
NPROCS=128  # 4 nodes × 32 cores

# Case directory (where this script is located)
CASE_DIR="${SLURM_SUBMIT_DIR}"

# OpenFOAM version and path (MODIFY FOR POLIMI HPC)
FOAM_VERSION="v2406"
FOAM_PATH="/opt/openfoam${FOAM_VERSION}/etc/bashrc"

#------------------------------------------------------------------------------
# Environment Setup
#------------------------------------------------------------------------------

echo "============================================================================"
echo "NANOSCALE LIFT SIMULATION - SLURM Job ${SLURM_JOB_ID}"
echo "============================================================================"
echo "Job started: $(date)"
echo "Running on: ${SLURM_JOB_NODELIST}"
echo "Number of nodes: ${SLURM_JOB_NUM_NODES}"
echo "Total MPI tasks: ${SLURM_NTASKS}"
echo "CPUs per task: ${SLURM_CPUS_PER_TASK}"
echo "Memory per node: ${SLURM_MEM_PER_NODE}"
echo "============================================================================"
echo ""

# Load modules (MODIFY FOR POLIMI HPC)
module purge
module load gcc/11.2.0
module load openmpi/4.1.1
# module load openfoam/v2406  # If available as module

# Source OpenFOAM environment
if [ -f "${FOAM_PATH}" ]; then
    source "${FOAM_PATH}"
    echo "OpenFOAM environment loaded: ${WM_PROJECT_DIR}"
else
    echo "ERROR: OpenFOAM environment not found at ${FOAM_PATH}"
    echo "Please modify FOAM_PATH in the script"
    exit 1
fi

# Change to case directory
cd "${CASE_DIR}" || exit 1
echo "Case directory: ${CASE_DIR}"
echo ""

#------------------------------------------------------------------------------
# Pre-flight checks
#------------------------------------------------------------------------------

echo "Pre-flight checks:"
echo "------------------------------------------------------------"

# Check if case is decomposed
if [ ! -d "processor0" ]; then
    echo "ERROR: Case not decomposed! Run decomposePar first."
    exit 1
fi

# Check number of processor directories matches NPROCS
NPROC_DIRS=$(find . -maxdepth 1 -type d -name "processor*" | wc -l)
if [ "${NPROC_DIRS}" -ne "${NPROCS}" ]; then
    echo "ERROR: Found ${NPROC_DIRS} processor directories, expected ${NPROCS}"
    echo "Re-run decomposePar with correct numberOfSubdomains"
    exit 1
fi

echo "✓ Case decomposed into ${NPROC_DIRS} subdomains"
echo "✓ MPI processes: ${NPROCS}"
echo ""

#------------------------------------------------------------------------------
# Performance and monitoring setup
#------------------------------------------------------------------------------

# Set OpenFOAM parallel run options
export FOAM_SIGFPE=false           # Disable FPE trapping for performance
export FOAM_SETNAN=false           # Disable NaN initialization
export OMPI_MCA_btl_openib_allow_ib=1  # Enable InfiniBand if available

# Set thread affinity
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "Environment variables set for optimal performance"
echo ""

#------------------------------------------------------------------------------
# Run simulation
#------------------------------------------------------------------------------

echo "============================================================================"
echo "Starting LIFT simulation..."
echo "============================================================================"
echo "Start time: $(date)"
echo ""

# Run with time tracking
time mpirun -np ${NPROCS} compInterFoam -parallel > log.compInterFoam 2>&1

# Capture exit status
EXIT_STATUS=$?

echo ""
echo "============================================================================"
echo "Simulation finished"
echo "============================================================================"
echo "End time: $(date)"
echo "Exit status: ${EXIT_STATUS}"
echo ""

#------------------------------------------------------------------------------
# Post-run summary
#------------------------------------------------------------------------------

if [ ${EXIT_STATUS} -eq 0 ]; then
    echo "✓ Simulation completed successfully!"

    # Extract final simulation time
    if [ -f "log.compInterFoam" ]; then
        FINAL_TIME=$(grep "^Time = " log.compInterFoam | tail -1 | awk '{print $3}')
        echo "Final simulation time: ${FINAL_TIME} seconds"

        # Count timesteps
        TIMESTEPS=$(grep "^Time = " log.compInterFoam | wc -l)
        echo "Total timesteps: ${TIMESTEPS}"

        # Check for warnings/errors
        WARNINGS=$(grep -i "warning" log.compInterFoam | wc -l)
        ERRORS=$(grep -i "error" log.compInterFoam | wc -l)
        echo "Warnings: ${WARNINGS}"
        echo "Errors: ${ERRORS}"
    fi
else
    echo "✗ Simulation failed with exit status ${EXIT_STATUS}"
    echo "Check log.compInterFoam for details"
fi

echo ""
echo "============================================================================"
echo "Job statistics"
echo "============================================================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Walltime used: $(squeue -j ${SLURM_JOB_ID} -h -o %M 2>/dev/null || echo 'N/A')"
echo "Memory used: $(sstat -j ${SLURM_JOB_ID} --format=MaxRSS -n 2>/dev/null || echo 'N/A')"
echo ""

#------------------------------------------------------------------------------
# Automatic post-processing (optional)
#------------------------------------------------------------------------------

if [ ${EXIT_STATUS} -eq 0 ]; then
    echo "Post-processing options:"
    echo "1. Reconstruct parallel data: reconstructPar"
    echo "2. Extract specific times: reconstructPar -time 0:2e-9"
    echo "3. Visualize: paraFoam"
    echo ""
    echo "To automatically reconstruct, uncomment the following line:"
    # mpirun -np ${NPROCS} reconstructPar -parallel
fi

echo "============================================================================"
echo "SLURM job ${SLURM_JOB_ID} finished: $(date)"
echo "============================================================================"

# EOF
