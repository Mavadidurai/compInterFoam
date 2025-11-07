#!/usr/bin/env python3
"""
Mesh Convergence Study - Case Setup Script
==========================================
Creates complete case directories for each mesh refinement level.

Author: Generated for compInterFoam mesh convergence study
Date: 2025-11-07
"""

import os
import sys
import shutil
from pathlib import Path


class CaseSetup:
    """Sets up case directories for mesh convergence study"""

    def __init__(self, base_case="LiftTest1", study_dir="meshStudy"):
        self.base_case_path = Path("..") / base_case
        self.study_dir = Path(study_dir)
        self.mesh_dir = Path("meshes")

        if not self.base_case_path.exists():
            raise FileNotFoundError(f"Base case not found: {self.base_case_path}")

    def setup_case(self, level_name):
        """Setup a single case directory for a refinement level"""

        case_dir = self.study_dir / level_name
        mesh_file = self.mesh_dir / f"blockMeshDict.{level_name}"

        if not mesh_file.exists():
            print(f"ERROR: Mesh file not found: {mesh_file}")
            print(f"       Run generateMeshes.py first!")
            return False

        print(f"\nSetting up case: {level_name}")
        print(f"  Target directory: {case_dir}")

        # Create case directory structure
        case_dir.mkdir(parents=True, exist_ok=True)

        # Copy directories from base case
        dirs_to_copy = ['0.orig', 'constant', 'system']

        for dir_name in dirs_to_copy:
            src = self.base_case_path / dir_name
            dst = case_dir / dir_name

            if src.exists():
                if dst.exists():
                    shutil.rmtree(dst)
                shutil.copytree(src, dst)
                print(f"  ✓ Copied {dir_name}/")
            else:
                print(f"  ✗ WARNING: {dir_name}/ not found in base case")

        # Replace blockMeshDict with the appropriate mesh
        system_dir = case_dir / "system"
        if system_dir.exists():
            block_mesh_dict = system_dir / "blockMeshDict"
            shutil.copy(mesh_file, block_mesh_dict)
            print(f"  ✓ Installed blockMeshDict for {level_name}")

        # Modify controlDict for shorter run time (for convergence study)
        self._modify_controlDict(case_dir, level_name)

        # Create Allrun script
        self._create_allrun_script(case_dir, level_name)

        # Create Allclean script
        self._create_allclean_script(case_dir)

        print(f"  ✓ Case {level_name} setup complete")

        return True

    def _modify_controlDict(self, case_dir, level_name):
        """Modify controlDict for convergence study"""

        controlDict_path = case_dir / "system" / "controlDict"

        if not controlDict_path.exists():
            print(f"  ✗ WARNING: controlDict not found")
            return

        # Read original
        with open(controlDict_path, 'r') as f:
            lines = f.readlines()

        # Modify parameters for shorter test run
        # For full convergence study, you may want to run longer
        modifications = {
            'endTime': '1e-10',  # 100 ps instead of 2 ns (for testing)
            'writeInterval': '5e-12',  # Write every 5 ps
        }

        modified = []
        for line in lines:
            modified_line = line
            for key, value in modifications.items():
                if line.strip().startswith(key):
                    # Keep the line structure but change the value
                    indent = len(line) - len(line.lstrip())
                    modified_line = ' ' * indent + f"{key}    {value};\n"
                    break
            modified.append(modified_line)

        # Write modified version
        with open(controlDict_path, 'w') as f:
            f.writelines(modified)

        print(f"  ✓ Modified controlDict (shorter run for convergence study)")

    def _create_allrun_script(self, case_dir, level_name):
        """Create Allrun script for the case"""

        allrun_path = case_dir / "Allrun"

        content = f"""#!/bin/sh
cd "${{0%/*}}" || exit 1    # Run from this directory

# Source OpenFOAM environment
if [ -f /opt/openfoam/etc/bashrc ]; then
    . /opt/openfoam/etc/bashrc
else
    echo "ERROR: OpenFOAM not found"
    exit 1
fi

echo ""
echo "=========================================="
echo "Running LIFT simulation: {level_name}"
echo "=========================================="
echo ""

# Clean previous results
./Allclean

# Generate mesh
echo "Generating mesh..."
blockMesh > log.blockMesh 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: blockMesh failed. Check log.blockMesh"
    exit 1
fi
echo "Mesh generated successfully"

# Check mesh quality
echo "Checking mesh quality..."
checkMesh > log.checkMesh 2>&1

# Initialize fields
echo "Initializing fields..."
cp -r 0.orig 0
topoSet > log.topoSet 2>&1
setFields > log.setFields 2>&1
echo "Fields initialized"

# Run solver
echo "Running compInterFoam..."
echo "This may take a while..."
compInterFoam > log.compInterFoam 2>&1 &
SOLVER_PID=$!

# Monitor progress
echo "Solver PID: $SOLVER_PID"
echo "Monitor with: tail -f log.compInterFoam"
echo ""

# Wait for completion
wait $SOLVER_PID
SOLVER_EXIT=$?

if [ $SOLVER_EXIT -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Simulation completed successfully!"
    echo "=========================================="
else
    echo ""
    echo "=========================================="
    echo "ERROR: Simulation failed (exit code: $SOLVER_EXIT)"
    echo "Check log.compInterFoam for details"
    echo "=========================================="
    exit $SOLVER_EXIT
fi

# Post-process
echo ""
echo "Extracting probe data..."
postProcess -func 'probes' > log.postProcess 2>&1

echo ""
echo "Run complete. Results in postProcessing/"
echo ""
"""

        with open(allrun_path, 'w') as f:
            f.write(content)

        # Make executable
        os.chmod(allrun_path, 0o755)

        print(f"  ✓ Created Allrun script")

    def _create_allclean_script(self, case_dir):
        """Create Allclean script for the case"""

        allclean_path = case_dir / "Allclean"

        content = """#!/bin/sh
cd "${0%/*}" || exit 1    # Run from this directory

echo "Cleaning case..."

# Remove time directories (keep 0.orig)
rm -rf [1-9]* 0 0.* processor* postProcessing

# Remove mesh
rm -rf constant/polyMesh

# Remove logs
rm -f log.*

# Remove dynamicCode
rm -rf dynamicCode

echo "Case cleaned"
"""

        with open(allclean_path, 'w') as f:
            f.write(content)

        # Make executable
        os.chmod(allclean_path, 0o755)

        print(f"  ✓ Created Allclean script")

    def setup_all_cases(self, levels=None):
        """Setup all case directories"""

        if levels is None:
            levels = ['coarse', 'medium', 'fine', 'very_fine']

        print("="*70)
        print("MESH CONVERGENCE STUDY - Setting Up Cases")
        print("="*70)

        success_count = 0
        for level in levels:
            if self.setup_case(level):
                success_count += 1

        print("\n" + "="*70)
        print(f"Setup complete: {success_count}/{len(levels)} cases ready")
        print("="*70)

        # Create master run script
        self._create_master_run_script(levels)

    def _create_master_run_script(self, levels):
        """Create a master script to run all cases"""

        script_path = self.study_dir / "runAll.sh"

        content = """#!/bin/bash
# Master script to run all mesh convergence cases

echo "========================================"
echo "MESH CONVERGENCE STUDY - Running All Cases"
echo "========================================"
echo ""

LEVELS=(""" + " ".join(levels) + """)

for level in "${LEVELS[@]}"; do
    echo ""
    echo "========================================"
    echo "Starting case: $level"
    echo "========================================"

    if [ -d "$level" ]; then
        cd "$level" || exit 1
        ./Allrun
        cd ..
    else
        echo "ERROR: Directory $level not found"
    fi
done

echo ""
echo "========================================"
echo "All simulations complete!"
echo "========================================"
echo ""
echo "Next step: Run ../analyzeConvergence.py to analyze results"
"""

        with open(script_path, 'w') as f:
            f.write(content)

        os.chmod(script_path, 0o755)

        print(f"\n  ✓ Created master run script: {script_path}")


def main():
    """Main execution"""

    levels = ['coarse', 'medium', 'fine', 'very_fine']

    if len(sys.argv) > 1:
        levels = sys.argv[1:]

    try:
        setup = CaseSetup()
        setup.setup_all_cases(levels=levels)

        print("\nNext steps:")
        print("  1. cd meshStudy")
        print("  2. Run individual cases: cd coarse && ./Allrun")
        print("  3. Or run all: ./runAll.sh")
        print("  4. Analyze: cd .. && python3 analyzeConvergence.py")

    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
