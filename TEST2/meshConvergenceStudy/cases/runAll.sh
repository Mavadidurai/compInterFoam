#!/bin/bash
# Master script to run all mesh convergence cases sequentially

echo "=========================================================================="
echo "Mesh Convergence Study - TEST2 (Graded Mesh)"
echo "Running all refinement levels sequentially"
echo "=========================================================================="
echo ""
echo "WARNING: This will take a LONG time (hours to days)!"
echo "Consider running cases individually or in parallel on HPC"
echo ""
read -p "Continue? (y/N) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled"
    exit 0
fi

# Get start time
start_time=$(date +%s)

# Run each case
for level in coarse medium fine very_fine; do
    echo ""
    echo "=========================================================================="
    echo "Running: $level"
    echo "=========================================================================="

    cd "$level" || exit 1
    ./Allrun
    exit_code=$?
    cd ..

    if [ $exit_code -ne 0 ]; then
        echo ""
        echo "ERROR: $level case failed!"
        echo "Check $level/log.compInterFoam for details"
        exit 1
    fi

    echo "✓ $level completed"
done

# Calculate total time
end_time=$(date +%s)
elapsed=$((end_time - start_time))
hours=$((elapsed / 3600))
minutes=$(((elapsed % 3600) / 60))

echo ""
echo "=========================================================================="
echo "All cases completed successfully!"
echo "Total time: ${hours}h ${minutes}m"
echo "=========================================================================="
echo ""
echo "Next step: Run analyzeConvergence.py to compute GCI values"
