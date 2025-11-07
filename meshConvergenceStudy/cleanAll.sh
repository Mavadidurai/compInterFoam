#!/bin/bash
###############################################################################
# Clean All - Remove all generated files from mesh convergence study
# Usage: ./cleanAll.sh [--keep-meshes]
###############################################################################

KEEP_MESHES=false

# Parse arguments
if [ "$1" == "--keep-meshes" ]; then
    KEEP_MESHES=true
fi

echo "========================================================================"
echo "CLEANING MESH CONVERGENCE STUDY"
echo "========================================================================"
echo ""

# Remove case directories
if [ -d "meshStudy" ]; then
    echo "Removing case directories..."
    rm -rf meshStudy
    echo "  ✓ meshStudy/ removed"
fi

# Remove mesh definitions (unless --keep-meshes specified)
if [ "$KEEP_MESHES" == "false" ]; then
    if [ -d "meshes" ]; then
        echo "Removing mesh definitions..."
        rm -rf meshes
        echo "  ✓ meshes/ removed"
    fi
else
    echo "Keeping mesh definitions (--keep-meshes specified)"
fi

# Remove convergence results
if [ -d "convergenceResults" ]; then
    echo "Removing convergence results..."
    rm -rf convergenceResults
    echo "  ✓ convergenceResults/ removed"
fi

# Remove Python cache
if [ -d "__pycache__" ]; then
    rm -rf __pycache__
    echo "  ✓ __pycache__/ removed"
fi

# Remove any .pyc files
find . -name "*.pyc" -delete 2>/dev/null
find . -name "*.pyo" -delete 2>/dev/null

echo ""
echo "========================================================================"
echo "Cleanup complete!"
echo "========================================================================"
echo ""

if [ "$KEEP_MESHES" == "false" ]; then
    echo "All generated files have been removed."
    echo "Run ./generateMeshes.py to start fresh."
else
    echo "Case directories and results removed, mesh definitions preserved."
    echo "Run ./setupCases.py to recreate cases."
fi

echo ""
