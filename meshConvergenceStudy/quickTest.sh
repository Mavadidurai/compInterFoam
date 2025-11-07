#!/bin/bash
###############################################################################
# Quick Test Script for Mesh Convergence Study
# ============================================
# This script performs a quick validation of the mesh convergence study setup
# without running full simulations (which can take hours/days).
#
# Usage: ./quickTest.sh
###############################################################################

set -e  # Exit on error

echo "========================================================================"
echo "MESH CONVERGENCE STUDY - QUICK TEST"
echo "========================================================================"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if we're in the right directory
if [ ! -f "generateMeshes.py" ]; then
    echo -e "${RED}ERROR: Must run from meshConvergenceStudy/ directory${NC}"
    exit 1
fi

# Step 1: Generate meshes
echo -e "${YELLOW}[1/4] Generating mesh definitions...${NC}"
python3 generateMeshes.py coarse medium
echo -e "${GREEN}✓ Mesh definitions generated${NC}"
echo ""

# Step 2: Setup cases (only coarse and medium for quick test)
echo -e "${YELLOW}[2/4] Setting up case directories...${NC}"
python3 setupCases.py coarse medium
echo -e "${GREEN}✓ Cases setup complete${NC}"
echo ""

# Step 3: Validate mesh generation (blockMesh only, no simulation)
echo -e "${YELLOW}[3/4] Validating mesh generation...${NC}"

for level in coarse medium; do
    echo "  Testing $level mesh..."
    cd meshStudy/$level

    # Check if blockMesh command is available
    if ! command -v blockMesh &> /dev/null; then
        echo -e "${YELLOW}    ⚠ blockMesh not found (OpenFOAM not loaded)${NC}"
        echo -e "${YELLOW}    ⚠ Skipping mesh validation${NC}"
        cd ../..
        continue
    fi

    # Run blockMesh
    blockMesh > log.blockMesh 2>&1

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}    ✓ $level mesh generated successfully${NC}"

        # Check mesh quality
        checkMesh > log.checkMesh 2>&1

        # Extract mesh statistics
        CELLS=$(grep "cells:" log.checkMesh | awk '{print $2}')
        NON_ORTH=$(grep "Max non-orthogonality" log.checkMesh | awk '{print $4}')

        echo "      Cells: $CELLS"
        echo "      Max non-orthogonality: $NON_ORTH"
    else
        echo -e "${RED}    ✗ $level mesh generation FAILED${NC}"
        echo "      Check meshStudy/$level/log.blockMesh"
    fi

    cd ../..
done

echo -e "${GREEN}✓ Mesh validation complete${NC}"
echo ""

# Step 4: Summary
echo -e "${YELLOW}[4/4] Test Summary${NC}"
echo "========================================================================"
echo ""
echo "Quick test completed successfully!"
echo ""
echo "Directory structure:"
tree -L 2 -d meshes meshStudy 2>/dev/null || find meshes meshStudy -type d -maxdepth 2
echo ""
echo "Mesh files generated:"
ls -lh meshes/blockMeshDict.* 2>/dev/null || echo "  (No mesh files found)"
echo ""
echo "========================================================================"
echo "NEXT STEPS:"
echo "========================================================================"
echo ""
echo "1. Generate all mesh levels (including fine and very_fine):"
echo "   ./generateMeshes.py"
echo "   ./setupCases.py"
echo ""
echo "2. Run simulations (WARNING: computationally expensive!):"
echo "   cd meshStudy"
echo "   ./runAll.sh"
echo ""
echo "3. After simulations complete, analyze results:"
echo "   python3 analyzeConvergence.py"
echo ""
echo "4. Review convergence report:"
echo "   cat convergenceResults/convergence_report.txt"
echo ""
echo "========================================================================"
echo ""
echo -e "${GREEN}Quick test complete! ✓${NC}"
echo ""
