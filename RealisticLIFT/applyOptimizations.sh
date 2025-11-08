#!/bin/bash
#------------------------------------------------------------------------------
# Apply Performance Optimizations to controlDict
#------------------------------------------------------------------------------

cd "${0%/*}" || exit

echo "========================================="
echo "RealisticLIFT Performance Optimization"
echo "========================================="
echo ""

# Backup original
if [ ! -f "system/controlDict.original" ]; then
    echo "Creating backup: system/controlDict.original"
    cp system/controlDict system/controlDict.original
    echo "✓ Backup created"
else
    echo "✓ Backup already exists"
fi

echo ""
echo "Applying optimizations..."
echo ""

# Apply optimizations using sed
cp system/controlDict system/controlDict.tmp

# 1. Increase base deltaT from 1e-13 to 5e-13
sed -i 's/^deltaT\s\+1e-13;/deltaT          5e-13;         /' system/controlDict.tmp

# 2. Increase minDeltaT from 1e-14 to 1e-13
sed -i 's/^minDeltaT\s\+1e-14;/minDeltaT       1e-13;         /' system/controlDict.tmp

# 3. Increase maxDeltaT from 2e-13 to 1e-12
sed -i 's/^maxDeltaT\s\+2e-13;/maxDeltaT       1e-12;         /' system/controlDict.tmp

# 4. Increase maxCo from 0.1 to 0.5
sed -i 's/^maxCo\s\+0\.1;/maxCo           0.5;           /' system/controlDict.tmp

# 5. Increase maxAlphaCo from 0.02 to 0.1
sed -i 's/^maxAlphaCo\s\+0\.02;/maxAlphaCo      0.1;           /' system/controlDict.tmp

# 6. Reduce write frequency from 1e-14 to 1e-12
sed -i 's/^writeInterval\s\+1e-14;/writeInterval   1e-12;         /' system/controlDict.tmp

# 7. Increase G coefficient at 300K from 1.0e18 to 5.0e18 (optional, commented out)
# sed -i 's/(300\s\+1\.0e18)/(300      5.0e18)   /' system/controlDict.tmp

# Move optimized file to place
mv system/controlDict.tmp system/controlDict

echo "✓ Optimizations applied"
echo ""
echo "CHANGES MADE:"
echo "-------------"
echo "1. deltaT:        1e-13 → 5e-13 s    (5x larger base time step)"
echo "2. minDeltaT:     1e-14 → 1e-13 s    (10x larger minimum)"
echo "3. maxDeltaT:     2e-13 → 1e-12 s    (5x larger maximum)"
echo "4. maxCo:         0.1 → 0.5          (less conservative)"
echo "5. maxAlphaCo:    0.02 → 0.1         (less aggressive compression)"
echo "6. writeInterval: 1e-14 → 1e-12 s    (100x less frequent output)"
echo ""
echo "EXPECTED RESULTS:"
echo "----------------"
echo "• 50-100x faster time stepping"
echo "• 100x less file I/O"
echo "• Overall ~100-500x speedup"
echo "• Estimated completion: 1-2 days (from 870 days)"
echo ""
echo "========================================="
echo "NEXT STEPS:"
echo "========================================="
echo ""
echo "1. Review changes:"
echo "   diff system/controlDict.original system/controlDict"
echo ""
echo "2. Clean previous run:"
echo "   ./Allclean"
echo ""
echo "3. Run simulation:"
echo "   # Serial (for testing):"
echo "   compInterFoam > log.compInterFoam 2>&1 &"
echo ""
echo "   # Parallel (recommended for speed):"
echo "   decomposePar"
echo "   mpirun -np 8 compInterFoam -parallel > log.compInterFoam 2>&1 &"
echo ""
echo "4. Monitor progress:"
echo "   tail -f log.compInterFoam"
echo ""
echo "5. Check temperatures periodically:"
echo "   grep 'max(Te):' log.compInterFoam | tail -1"
echo "   grep 'max(Tl):' log.compInterFoam | tail -1"
echo ""
echo "========================================="
echo ""
echo "To restore original settings:"
echo "  cp system/controlDict.original system/controlDict"
echo ""
