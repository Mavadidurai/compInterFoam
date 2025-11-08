#!/bin/bash
#------------------------------------------------------------------------------
# Simulation Diagnostics and Performance Analysis
#------------------------------------------------------------------------------

LOG_FILE="log.compInterFoam"

if [ ! -f "$LOG_FILE" ]; then
    echo "Error: log.compInterFoam not found"
    exit 1
fi

echo "========================================="
echo "SIMULATION PERFORMANCE DIAGNOSTICS"
echo "========================================="
echo ""

# Extract time steps
echo "TIME STEP ANALYSIS:"
echo "-------------------"
TIMES=$(grep "^Time = " "$LOG_FILE" | awk '{print $3}' | tail -5)
DELTAT=$(grep "^deltaT = " "$LOG_FILE" | tail -1 | awk '{print $3}')
echo "Last 5 simulation times:"
echo "$TIMES"
echo ""
echo "Current deltaT: $DELTAT s"

# Calculate progress
START_TIME=0
END_TIME=$(grep "endTime" system/controlDict | grep -v "//" | awk '{print $2}' | tr -d ';')
CURRENT_TIME=$(grep "^Time = " "$LOG_FILE" | tail -1 | awk '{print $3}')

if [ -n "$CURRENT_TIME" ] && [ -n "$END_TIME" ]; then
    PROGRESS=$(echo "scale=6; $CURRENT_TIME / $END_TIME * 100" | bc -l)
    REMAINING=$(echo "scale=10; $END_TIME - $CURRENT_TIME" | bc -l)
    echo "Progress: $PROGRESS %"
    echo "Remaining time: $REMAINING s"
fi

echo ""
echo "EXECUTION TIME:"
echo "---------------"
EXEC_TIME=$(grep "ExecutionTime" "$LOG_FILE" | tail -1 | awk '{print $3}')
CLOCK_TIME=$(grep "ClockTime" "$LOG_FILE" | tail -1 | awk '{print $8}')
echo "Execution time: $EXEC_TIME s"
echo "Clock time: $CLOCK_TIME s"

# Time steps completed
N_STEPS=$(grep "^Time = " "$LOG_FILE" | wc -l)
echo "Time steps completed: $N_STEPS"

if [ -n "$EXEC_TIME" ] && [ -n "$N_STEPS" ] && [ "$N_STEPS" -gt 0 ]; then
    TIME_PER_STEP=$(echo "scale=2; $EXEC_TIME / $N_STEPS" | bc -l)
    echo "Average time per step: $TIME_PER_STEP s/step"

    if [ -n "$REMAINING" ] && [ -n "$DELTAT" ]; then
        STEPS_REMAINING=$(echo "scale=0; $REMAINING / $DELTAT" | bc -l)
        ETA=$(echo "scale=0; $STEPS_REMAINING * $TIME_PER_STEP" | bc -l)
        ETA_HOURS=$(echo "scale=1; $ETA / 3600" | bc -l)
        ETA_DAYS=$(echo "scale=1; $ETA / 86400" | bc -l)
        echo "Estimated steps remaining: $STEPS_REMAINING"
        echo "Estimated time to completion: $ETA s (~$ETA_HOURS hours / ~$ETA_DAYS days)"
    fi
fi

echo ""
echo "TEMPERATURE ANALYSIS:"
echo "--------------------"
echo "Last reported temperatures:"
grep "max(Te):" "$LOG_FILE" | tail -1
grep "max(Tl):" "$LOG_FILE" | tail -1

MAX_TE=$(grep "max(Te):" "$LOG_FILE" | tail -1 | awk '{print $2}')
MAX_TL=$(grep "max(Tl):" "$LOG_FILE" | tail -1 | awk '{print $2}')
echo ""
echo "Peak Te reached: $MAX_TE K"
echo "Peak Tl reached: $MAX_TL K"
echo "Vaporization threshold: 3560 K"

echo ""
echo "LASER ENERGY:"
echo "-------------"
grep "Cumulative absorbed:" "$LOG_FILE" | tail -1
grep "Laser input:" "$LOG_FILE" | tail -1

echo ""
echo "PRESSURE ANALYSIS:"
echo "------------------"
grep "Max |recoilPressure|" "$LOG_FILE" | tail -5 | tail -1

echo ""
echo "COURANT NUMBERS:"
echo "----------------"
grep "Courant Number mean:" "$LOG_FILE" | tail -1

echo ""
echo "SOLVER PERFORMANCE:"
echo "-------------------"
echo "Pressure solver iterations:"
grep "Solving for p_rgh" "$LOG_FILE" | tail -20 | head -10

echo ""
echo "MESH INFO:"
echo "----------"
if [ -d "constant/polyMesh" ]; then
    NCELLS=$(grep "nCells:" constant/polyMesh/boundary 2>/dev/null | head -1 | awk '{print $2}')
    if [ -z "$NCELLS" ]; then
        # Try to count from owner file
        if [ -f "constant/polyMesh/owner" ]; then
            NCELLS=$(grep -v "^/\|^$\|^(\|^)" constant/polyMesh/owner | wc -l)
            echo "Estimated cells: $NCELLS"
        fi
    else
        echo "Number of cells: $NCELLS"
    fi
fi

echo ""
echo "RECOMMENDATIONS:"
echo "================"

if [ -n "$ETA_DAYS" ]; then
    if (( $(echo "$ETA_DAYS > 1" | bc -l) )); then
        echo "⚠️  CRITICAL: Simulation will take $ETA_DAYS days at current rate!"
        echo ""
        echo "Suggested actions:"
        echo "1. INCREASE time step (minDeltaT and deltaT in controlDict)"
        echo "2. Run in PARALLEL (decompose mesh, use mpirun)"
        echo "3. Reduce output frequency (writeInterval)"
        echo "4. Consider coarser mesh if resolution allows"
    fi
fi

if [ -n "$MAX_TE" ]; then
    if (( $(echo "$MAX_TE < 1000" | bc -l) )); then
        echo "⚠️  Electron temperature not rising fast enough"
        echo "   Check: laser focus, absorption coefficient, coupling parameters"
    fi
fi

if [ -n "$MAX_TL" ]; then
    if (( $(echo "$MAX_TL < 500" | bc -l) )); then
        echo "⚠️  Lattice temperature very low - heating too slow"
        echo "   Check: electron-lattice coupling coefficient G"
    fi
fi

echo ""
echo "========================================="
