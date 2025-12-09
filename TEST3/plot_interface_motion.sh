#!/bin/bash
# Extract alpha.metal volume over time
echo "Time(ps) MetalVolume(µm³)"
for dir in [0-9]*; do
    if [ -f "$dir/alpha.metal" ]; then
        time=$(echo "$dir" | sed 's/e/*10^/g')
        vol=$(grep "alpha.metal volume" "$dir/alpha.metal" 2>/dev/null | awk '{print $4}')
        echo "$time $vol"
    fi
done
