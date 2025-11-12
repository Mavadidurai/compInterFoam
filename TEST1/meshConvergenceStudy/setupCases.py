#!/usr/bin/env python3
"""
Setup case directories for mesh convergence study - TEST1

This script creates complete OpenFOAM case directories for each refinement level.
It copies the TEST1 case structure and installs the appropriate blockMeshDict.
"""

import os
import sys
import shutil
from pathlib import Path

# Refinement levels to setup
REFINEMENT_LEVELS = ['coarse', 'medium', 'fine', 'very_fine']

def copy_directory(src, dst, exclude=None):
    """Recursively copy directory, excluding specified patterns"""
    if exclude is None:
        exclude = []

    if os.path.exists(dst):
        shutil.rmtree(dst)

    shutil.copytree(src, dst, ignore=shutil.ignore_patterns(*exclude))

def setup_case(level, base_case, meshes_dir, cases_dir):
    """Setup a single case directory for a refinement level"""

    # Paths
    case_dir = os.path.join(cases_dir, level)
    mesh_file = os.path.join(meshes_dir, f'blockMeshDict.{level}')

    # Check if mesh file exists
    if not os.path.exists(mesh_file):
        print(f"ERROR: Mesh file not found: {mesh_file}")
        print(f"       Run generateMeshes.py first!")
        return False

    print(f"Setting up case: {level}")
    print(f"  Source: {base_case}")
    print(f"  Target: {case_dir}")

    # Copy base case structure
    # Exclude time directories (0, 0.orig might exist), processor directories, logs, and meshConvergenceStudy
    exclude_patterns = ['0', '0.*', 'processor*', '*.log', 'log.*',
                        'postProcessing', 'dynamicCode', 'PyFoam*', 'meshConvergenceStudy']

    try:
        copy_directory(base_case, case_dir, exclude=exclude_patterns)
    except Exception as e:
        print(f"  ERROR copying case: {e}")
        return False

    # Copy 0.orig to 0 if it exists
    orig_dir = os.path.join(base_case, '0.orig')
    if os.path.exists(orig_dir):
        zero_dir = os.path.join(case_dir, '0')
        shutil.copytree(orig_dir, zero_dir)
        print(f"  Copied initial conditions from 0.orig -> 0")
    else:
        print(f"  WARNING: No 0.orig directory found in base case")

    # Install blockMeshDict for this refinement level
    target_mesh = os.path.join(case_dir, 'system', 'blockMeshDict')
    shutil.copy2(mesh_file, target_mesh)
    print(f"  Installed mesh: blockMeshDict.{level}")

    # Create Allrun script
    allrun_content = f"""#!/bin/sh
cd "${{0%/*}}" || exit 1    # Run from this directory

# Source OpenFOAM environment
if [ -f "/opt/openfoam/etc/bashrc" ]; then
    . /opt/openfoam/etc/bashrc
else
    echo "ERROR: OpenFOAM environment not found"
    exit 1
fi

echo "Running mesh convergence study - TEST1 - {level}"
echo "=================================================="
echo ""

# Clean previous results
echo "Cleaning previous results..."
./Allclean

# Generate mesh
echo "Generating mesh..."
blockMesh > log.blockMesh 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: blockMesh failed! Check log.blockMesh"
    exit 1
fi
echo "Mesh generated successfully"
echo ""

# Check mesh
echo "Checking mesh quality..."
checkMesh > log.checkMesh 2>&1
if [ $? -ne 0 ]; then
    echo "WARNING: checkMesh reported issues. Check log.checkMesh"
fi
echo ""

# Run solver
echo "Running compInterFoam..."
echo "This may take several hours depending on mesh size and hardware..."
echo "Monitor progress with: tail -f log.compInterFoam"
echo ""
compInterFoam > log.compInterFoam 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: compInterFoam failed! Check log.compInterFoam"
    exit 1
fi

echo ""
echo "Simulation completed successfully!"
echo "Results available in time directories"
echo "=================================================="
"""

    allrun_path = os.path.join(case_dir, 'Allrun')
    with open(allrun_path, 'w') as f:
        f.write(allrun_content)
    os.chmod(allrun_path, 0o755)
    print(f"  Created Allrun script")

    # Create Allclean script
    allclean_content = """#!/bin/sh
cd "${0%/*}" || exit 1    # Run from this directory

# Remove time directories (except 0)
echo "Removing time directories..."
for f in [1-9]* processor*; do
    if [ -e "$f" ]; then
        rm -rf "$f"
    fi
done

# Remove logs
echo "Removing logs..."
rm -f log.* *.log

# Remove mesh
echo "Removing mesh..."
rm -rf constant/polyMesh

# Remove post-processing
echo "Removing post-processing data..."
rm -rf postProcessing

# Remove dynamic code
rm -rf dynamicCode

echo "Case cleaned"
"""

    allclean_path = os.path.join(case_dir, 'Allclean')
    with open(allclean_path, 'w') as f:
        f.write(allclean_content)
    os.chmod(allclean_path, 0o755)
    print(f"  Created Allclean script")

    print(f"  ✓ Case setup complete")
    print()

    return True

def create_master_run_script(cases_dir):
    """Create a master script to run all cases sequentially"""

    script_content = """#!/bin/bash
# Master script to run all mesh convergence cases sequentially

echo "=========================================================================="
echo "Mesh Convergence Study - TEST1 (Uniform Mesh)"
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
"""

    script_path = os.path.join(cases_dir, 'runAll.sh')
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    print(f"Created master run script: {script_path}")

def main():
    """Main function to setup all cases"""

    # Determine paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_case = os.path.dirname(script_dir)  # Parent directory (TEST1)
    meshes_dir = os.path.join(script_dir, 'meshes')
    cases_dir = os.path.join(script_dir, 'cases')

    print("="*70)
    print("Mesh Convergence Study - TEST1 Case Setup")
    print("="*70)
    print()
    print(f"Base case: {base_case}")
    print(f"Meshes:    {meshes_dir}")
    print(f"Output:    {cases_dir}")
    print()

    # Check if meshes exist
    if not os.path.exists(meshes_dir):
        print("ERROR: Meshes directory not found!")
        print("       Run generateMeshes.py first to create mesh files")
        sys.exit(1)

    # Create cases directory
    os.makedirs(cases_dir, exist_ok=True)

    # Setup each case
    success_count = 0
    for level in REFINEMENT_LEVELS:
        if setup_case(level, base_case, meshes_dir, cases_dir):
            success_count += 1

    print("="*70)
    if success_count == len(REFINEMENT_LEVELS):
        print(f"✓ All {success_count} cases setup successfully!")
        print()

        # Create master run script
        create_master_run_script(cases_dir)
        print()

        print("Next steps:")
        print(f"  1. cd {cases_dir}")
        print("  2. Run individual cases: cd <level> && ./Allrun")
        print("  3. Or run all: ./runAll.sh")
        print("  4. After completion: cd .. && python3 analyzeConvergence.py")
    else:
        print(f"WARNING: Only {success_count}/{len(REFINEMENT_LEVELS)} cases setup successfully")
        print("         Check errors above")
    print("="*70)

if __name__ == "__main__":
    main()
