#!/bin/bash
# Quick fix script for ejection blocking issue
# This increases maxPressure to allow shock waves to develop fully

set -e

echo "================================================"
echo "Applying Ejection Fix to TEST1/system/fvSolution"
echo "================================================"
echo ""
echo "Current maxPressure: 5.0e10 Pa (50 GPa)"
echo "New maxPressure:     1.5e11 Pa (150 GPa)"
echo ""
echo "Reason: The pressure clamp at 40.6 GPa is too close"
echo "        to the 27 GPa recoil pressure, leaving only"
echo "        50% headroom for shock wave dynamics."
echo ""

# Backup original file
cp TEST1/system/fvSolution TEST1/system/fvSolution.backup_$(date +%Y%m%d_%H%M%S)
echo "✓ Backed up original fvSolution"

# Apply fix using sed
sed -i 's/maxPressure        5\.0e10/maxPressure        1.5e11   \/\/ 150 GPa - increased for shock wave headroom/' TEST1/system/fvSolution
sed -i 's/minPressure       -5\.e10/minPressure       -1.5e11  \/\/ Symmetric bounds/' TEST1/system/fvSolution

echo "✓ Updated pressure clamp limits"
echo ""
echo "Changes applied:"
grep -A1 "maxPressure" TEST1/system/fvSolution | head -3
echo ""
echo "================================================"
echo "Fix applied successfully!"
echo "================================================"
echo ""
echo "Next steps:"
echo "  1. Resume your simulation: fg"
echo "  2. Let it run for 10-20 picoseconds"
echo "  3. Watch for:"
echo "     - Volume loss > 1%"
echo "     - Metal velocity > 2000 m/s"
echo "     - Velocity ratio > 0.7"
echo ""
echo "If simulation was stopped, restart with:"
echo "  cd ~/compInterFoam/TEST1"
echo "  compInterFoam"
echo ""
