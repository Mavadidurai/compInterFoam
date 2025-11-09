#!/bin/bash

echo "=== CLEAN RESTART FOR ADAPTIVE TIMESTEPPING ==="

# 1. Kill any running simulations
echo "1. Killing any running compInterFoam processes..."
pkill -9 compInterFoam
sleep 2

# 2. Clean time directories (keep only 0)
echo "2. Cleaning old time directories..."
find . -maxdepth 1 -type d -regex '\./[0-9]+\.?[0-9]*' ! -name "0" ! -name "0.orig" -exec rm -rf {} \; 2>/dev/null

# 3. Show current settings
echo ""
echo "3. Current controlDict settings:"
grep -E "startFrom|deltaT|maxDeltaT|writeInterval|adjustTimeStep" system/controlDict | grep -v "//"

# 4. Start simulation
echo ""
echo "4. Starting simulation (will show deltaT growth)..."
echo "   Watch for 'deltaT =' to verify it's growing!"
echo ""

compInterFoam 2>&1 | tee log.restart | grep --line-buffered -E "Time =|deltaT =|Courant Number mean"
