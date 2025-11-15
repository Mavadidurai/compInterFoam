#!/bin/bash
# Quick debugging script for parallel hang in TEST2
# Usage: ./debug_parallel.sh [test_number]

set -e

CASE_DIR="$HOME/OpenFOAM/mavadi-v2406/run/TEST2"
cd "$CASE_DIR" || exit 1

echo "===== OpenFOAM Parallel Debug Helper ====="
echo "Case: $CASE_DIR"
echo ""

# Function to clean processor directories
clean_processors() {
    echo "[1] Cleaning processor directories..."
    rm -rf processor* dynamicCode
    echo "    Done."
}

# Function to check field consistency
check_fields() {
    echo "[2] Checking field consistency..."
    ORIG_COUNT=$(ls -1 0/ | wc -l)
    echo "    Original fields in 0/: $ORIG_COUNT"

    if [ -d "processor0" ]; then
        PROC_COUNT=$(ls -1 processor0/0/ | wc -l)
        echo "    Processor fields in processor0/0/: $PROC_COUNT"

        if [ "$ORIG_COUNT" -ne "$PROC_COUNT" ]; then
            echo "    WARNING: Field count mismatch!"
            echo "    Missing fields in processor directories."
            return 1
        else
            echo "    OK: Field counts match."
        fi
    else
        echo "    No processor directories found."
    fi
}

# Test 1: Simplest possible run
test_1_simple() {
    echo ""
    echo "===== TEST 1: Minimal Configuration ====="
    echo "Disabling complex features for basic parallel test..."

    clean_processors

    # Backup controlDict
    cp system/controlDict system/controlDict.backup

    # Create simplified controlDict
    cat > system/controlDict.simple <<'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}

application     compInterFoam;
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         1e-14;        // Very short test
deltaT          2e-15;
writeInterval   1e-14;
adjustTimeStep  no;           // Disable adaptive timestep
purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   6;
fileHandler     collated;     // Use collated for simplicity

runTimeModifiable no;         // Disable runtime modification

// Disable advanced features
useAdvancedInterfaceCapturing  false;
enableLiftProcessTracker false;
enableTClamp    false;
verbose         false;

// Keep minimal physics
phaseChangeCoeffs
{
    model               clausius_clapeyron;
    hf                  [0 2 -2 0 0 0 0] 9.1e6;
    gasConstant         [0 2 -2 -1 0 0 0] 174;
    Tsol                [0 0 0 1 0 0 0] 1941;
    Tvap                [0 0 0 1 0 0 0] 2200;
    evaporationCoeff    0.03;
    alphaMin            0.001;
    dtFloor             1e-15;
    relaxationTime      1e-11;
    maxSource           [1 -1 -3 0 0 0 0] 1e22;
    minCoefficient      1e6;
    onlyAboveVapor      false;
    activationTime      ((0 2e-10));
}

massTransferCoeffs
{
    tStart  (0);
    tEnd    (2e-10);
    rateMax 3e13;
}

twoTemperatureProperties
{
    Ce
    {
        type        polynomial;
        value       0;
        coeffs      ((0 0) (630 1));
    }
    G
    {
        type table;
        values
        (
            (300   1.0e18)
            (1000  3.0e18)
            (2000  6.0e18)
            (3000  1.0e19)
            (5000  2.0e19)
            (10000 5.0e19)
            (20000 1.0e20)
        );
    }
    Cl  [1 -1 -2 -1 0 0 0] 2.5e6;
    De  [0  2 -1  0 0 0 0] 1e-4;
    gasMetalExchangeCoeff
    {
        type        kapitza;
        Z_metal     2.3e7;
        Z_gas       383;
        rho_metal   4500;
        rho_gas     1.7;
        maxTemperature 0;
    }
    gasMetalCutoffTemperature 2200;
    minMetalFractionForExchange 0.5;
    minGasMetalExchangeCoeff [1 -1 -3 -1 0 0 0] 1e8;
    maxGasMetalExchangeCoeff [1 -1 -3 -1 0 0 0] 5e9;
    maxGasMetalHeatFlux    5e9;
    minTe           200;
    maxTe           20000;
    maxTl           10000;
    ambientTemperature 300;
    temperatureDifferenceFloor 10.0;
    electronSubCycles        1;
    minElectronSubCycles     1;
    maxElectronSubCycles     20;
    maxElectronDeltaT        [0 0 1 0 0 0 0] 1e-15;
    nInnerCouplingSweeps         2;
    innerCouplingResidualTol     1e-3;
    klHighTThreshold 1000;
    klExponent       0.5;
    metalFractionFloor      1e-9;
    metalFractionCutoff     1e-6;
    metalAmbientBlendWidth  1e-3;
    CvolGasFloor            1e6;
    energyTolerance         0.02;
    energyDiagnostics       false;
    temperatureDiagnostics  false;
    energyAudit             false;
}
EOF

    mv system/controlDict.simple system/controlDict

    echo "Decomposing mesh..."
    decomposePar > /dev/null 2>&1

    echo "Running with 2 processors (minimal test)..."
    timeout 30 mpirun --mca btl tcp,self --bind-to none -np 2 compInterFoam -parallel 2>&1 | tee log.test1

    # Restore original
    mv system/controlDict.backup system/controlDict

    echo ""
    echo "Test 1 complete. Check log.test1 for results."
}

# Test 2: Without function objects
test_2_no_functions() {
    echo ""
    echo "===== TEST 2: Disable Function Objects ====="

    clean_processors

    cp system/controlDict system/controlDict.backup

    # Remove function objects section
    sed -i '/^functions/,/^}/d' system/controlDict

    echo "Decomposing mesh..."
    decomposePar > /dev/null 2>&1

    echo "Running with 4 processors..."
    timeout 60 mpirun --mca btl tcp,self -np 4 compInterFoam -parallel 2>&1 | tee log.test2

    mv system/controlDict.backup system/controlDict

    echo ""
    echo "Test 2 complete. Check log.test2 for results."
}

# Test 3: Different file handler
test_3_file_handler() {
    echo ""
    echo "===== TEST 3: Different File Handler ====="

    clean_processors

    cp system/controlDict system/controlDict.backup

    # Change file handler
    sed -i 's/fileHandler.*uncollated.*/fileHandler     collated;/' system/controlDict

    echo "Decomposing mesh..."
    decomposePar > /dev/null 2>&1

    echo "Running with 4 processors (collated file handler)..."
    timeout 60 mpirun -np 4 compInterFoam -parallel 2>&1 | tee log.test3

    mv system/controlDict.backup system/controlDict

    echo ""
    echo "Test 3 complete. Check log.test3 for results."
}

# Test 4: Serial run
test_4_serial() {
    echo ""
    echo "===== TEST 4: Serial Run ====="

    # Backup and remove processor dirs
    if [ -d "processor0" ]; then
        mkdir -p processor_backup
        mv processor* processor_backup/
    fi

    echo "Running in serial mode..."
    timeout 60 compInterFoam 2>&1 | tee log.serial

    echo ""
    echo "Test 4 complete. Check log.serial for results."
}

# Test 5: Check for deadlock patterns
test_5_deadlock_check() {
    echo ""
    echo "===== TEST 5: Deadlock Detection ====="

    clean_processors

    echo "Decomposing mesh..."
    decomposePar > /dev/null 2>&1

    echo "Starting parallel run in background..."
    mpirun -np 4 compInterFoam -parallel > log.deadlock 2>&1 &
    PID=$!

    echo "Monitoring for 20 seconds..."
    for i in {1..20}; do
        sleep 1

        # Check if still running
        if ! kill -0 $PID 2>/dev/null; then
            echo "Process completed."
            break
        fi

        # Check CPU usage
        CPU=$(ps -p $PID -o %cpu= 2>/dev/null || echo "0")
        echo "  [$i/20] CPU usage: $CPU%"

        if (( $(echo "$CPU < 1.0" | bc -l) )); then
            echo "  WARNING: Low CPU usage - possible deadlock!"
        fi
    done

    # Kill if still running
    if kill -0 $PID 2>/dev/null; then
        echo ""
        echo "Process still running after 20s - killing..."
        kill -9 $PID 2>/dev/null
        echo "RESULT: Likely deadlock detected (process hung)"
    else
        echo ""
        echo "RESULT: Process completed normally"
    fi

    echo ""
    echo "Test 5 complete. Check log.deadlock for details."
}

# Main menu
if [ -z "$1" ]; then
    echo "Select a test:"
    echo "  1 - Minimal configuration (safest)"
    echo "  2 - Disable function objects"
    echo "  3 - Different file handler"
    echo "  4 - Serial run (baseline)"
    echo "  5 - Deadlock detection"
    echo "  all - Run all tests"
    echo ""
    echo "Usage: $0 [1|2|3|4|5|all]"
    exit 0
fi

case "$1" in
    1) test_1_simple ;;
    2) test_2_no_functions ;;
    3) test_3_file_handler ;;
    4) test_4_serial ;;
    5) test_5_deadlock_check ;;
    all)
        test_1_simple
        test_2_no_functions
        test_3_file_handler
        test_4_serial
        test_5_deadlock_check
        ;;
    *)
        echo "Invalid test number: $1"
        echo "Use: $0 [1|2|3|4|5|all]"
        exit 1
        ;;
esac

echo ""
echo "===== Debug Complete ====="
echo "Log files generated: log.test* or log.serial"
