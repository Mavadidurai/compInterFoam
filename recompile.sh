#!/bin/bash
# Recompile compInterFoam solver after source code changes

echo "========================================="
echo "Recompiling compInterFoam solver..."
echo "========================================="

# Check if OpenFOAM environment is loaded
if [ -z "$WM_PROJECT" ]; then
    echo "ERROR: OpenFOAM environment not loaded!"
    echo ""
    echo "Please source your OpenFOAM environment first:"
    echo "    source /opt/openfoam2406/etc/bashrc"
    echo "    (or wherever your OpenFOAM is installed)"
    echo ""
    exit 1
fi

echo "OpenFOAM version: $WM_PROJECT_VERSION"
echo "Compiler: $WM_COMPILER"
echo ""

# Clean previous build
echo "Cleaning previous build..."
wclean

# Compile
echo ""
echo "Compiling solver..."
wmake

if [ $? -eq 0 ]; then
    echo ""
    echo "========================================="
    echo "✓ Compilation successful!"
    echo "========================================="
    echo ""
    echo "The updated solver is ready to use."
    echo "You can now run your simulation from TEST1:"
    echo "    cd TEST1"
    echo "    compInterFoam"
else
    echo ""
    echo "========================================="
    echo "✗ Compilation failed!"
    echo "========================================="
    echo ""
    echo "Please check the error messages above."
    exit 1
fi
