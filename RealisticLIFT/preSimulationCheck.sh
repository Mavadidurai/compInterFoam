#!/bin/bash
#------------------------------------------------------------------------------
# Pre-Simulation Validation Script for RealisticLIFT Case
# Checks case structure, mesh, initial conditions, physical properties,
# numerical settings, and simulation parameters
#------------------------------------------------------------------------------

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Counters
ERRORS=0
WARNINGS=0
CHECKS=0

# Helper functions
error() {
    echo -e "${RED}[ERROR]${NC} $1"
    ((ERRORS++))
}

warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
    ((WARNINGS++))
}

success() {
    echo -e "${GREEN}[OK]${NC} $1"
}

info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

check() {
    ((CHECKS++))
}

section() {
    echo ""
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
}

# Get case directory
CASE_DIR="${0%/*}"
cd "$CASE_DIR" || exit 1

section "1. CASE STRUCTURE VALIDATION"

# Check required directories
check
if [ -d "0.orig" ]; then
    success "0.orig directory exists"
else
    error "0.orig directory not found"
fi

check
if [ -d "constant" ]; then
    success "constant directory exists"
else
    error "constant directory not found"
fi

check
if [ -d "system" ]; then
    success "system directory exists"
else
    error "system directory not found"
fi

# Check required files in system
section "2. SYSTEM DIRECTORY FILES"

SYSTEM_FILES=("controlDict" "fvSchemes" "fvSolution" "blockMeshDict")
for file in "${SYSTEM_FILES[@]}"; do
    check
    if [ -f "system/$file" ]; then
        success "system/$file exists"
    else
        error "system/$file not found"
    fi
done

# Check required files in constant
section "3. CONSTANT DIRECTORY FILES"

CONSTANT_FILES=("transportProperties" "thermophysicalProperties" "thermophysicalProperties.metal" "thermophysicalProperties.air" "laserProperties" "g" "turbulenceProperties")
for file in "${CONSTANT_FILES[@]}"; do
    check
    if [ -f "constant/$file" ]; then
        success "constant/$file exists"
    else
        error "constant/$file not found"
    fi
done

# Check required fields in 0.orig
section "4. INITIAL CONDITIONS (0.orig)"

REQUIRED_FIELDS=("alpha.metal" "alpha.air" "T" "Te" "Tl" "U" "p" "p_rgh")
for field in "${REQUIRED_FIELDS[@]}"; do
    check
    if [ -f "0.orig/$field" ]; then
        success "0.orig/$field exists"
    else
        error "0.orig/$field not found"
    fi
done

# Check mesh
section "5. MESH VALIDATION"

check
if [ -d "constant/polyMesh" ]; then
    success "Mesh exists in constant/polyMesh"

    # Check mesh files
    MESH_FILES=("points" "faces" "owner" "neighbour" "boundary")
    for mfile in "${MESH_FILES[@]}"; do
        check
        if [ -f "constant/polyMesh/$mfile" ]; then
            success "constant/polyMesh/$mfile exists"
        else
            error "constant/polyMesh/$mfile not found"
        fi
    done

    # Try to get mesh statistics
    if command -v checkMesh &> /dev/null; then
        info "Running checkMesh for mesh quality validation..."
        checkMesh -case . > checkMesh.log 2>&1
        if [ $? -eq 0 ]; then
            success "checkMesh completed - see checkMesh.log for details"
            # Extract key metrics
            if grep -q "Mesh OK" checkMesh.log; then
                success "Mesh quality: OK"
            else
                warning "Mesh may have quality issues - check checkMesh.log"
            fi
        else
            warning "checkMesh failed - see checkMesh.log"
        fi
    else
        warning "checkMesh not available in PATH - skipping mesh quality check"
    fi
else
    warning "Mesh not found - run blockMesh before simulation"
    info "To generate mesh: blockMesh"
fi

# Validate numerical settings
section "6. NUMERICAL SETTINGS VALIDATION"

check
if [ -f "system/controlDict" ]; then
    # Extract time step settings
    deltaT=$(grep -E "^\s*deltaT\s+" system/controlDict | awk '{print $2}' | tr -d ';')
    maxCo=$(grep -E "^\s*maxCo\s+" system/controlDict | awk '{print $2}' | tr -d ';')
    maxAlphaCo=$(grep -E "^\s*maxAlphaCo\s+" system/controlDict | awk '{print $2}' | tr -d ';')
    startTime=$(grep -E "^\s*startTime\s+" system/controlDict | awk '{print $2}' | tr -d ';')
    endTime=$(grep -E "^\s*endTime\s+" system/controlDict | awk '{print $2}' | tr -d ';')

    info "Time step settings:"
    info "  deltaT: $deltaT"
    info "  maxCo: $maxCo"
    info "  maxAlphaCo: $maxAlphaCo"
    info "  startTime: $startTime"
    info "  endTime: $endTime"

    # Validate Courant numbers
    if [ -n "$maxCo" ]; then
        if (( $(echo "$maxCo > 1" | bc -l) )); then
            warning "maxCo = $maxCo > 1.0 may cause instability"
        else
            success "maxCo = $maxCo is within stable range"
        fi
    fi

    if [ -n "$maxAlphaCo" ]; then
        if (( $(echo "$maxAlphaCo > 0.5" | bc -l) )); then
            warning "maxAlphaCo = $maxAlphaCo > 0.5 may cause interface smearing"
        else
            success "maxAlphaCo = $maxAlphaCo is conservative for interface capturing"
        fi
    fi
fi

# Validate laser properties
section "7. LASER PROPERTIES VALIDATION"

check
if [ -f "constant/laserProperties" ]; then
    pulseEnergy=$(grep -E "^\s*pulseEnergy\s+" constant/laserProperties | awk '{print $2}' | tr -d ';')
    pulseWidth=$(grep -E "^\s*pulseWidth\s+" constant/laserProperties | awk '{print $2}' | tr -d ';')
    wavelength=$(grep -E "^\s*wavelength\s+" constant/laserProperties | awk '{print $2}' | tr -d ';')
    spotSize=$(grep -E "^\s*spotSize\s+" constant/laserProperties | awk '{print $2}' | tr -d ';')
    absorptionCoeff=$(grep -E "^\s*absorptionCoeff\s+" constant/laserProperties | awk '{print $2}' | tr -d ';')

    info "Laser parameters:"
    info "  Pulse energy: $pulseEnergy J"
    info "  Pulse width: $pulseWidth s"
    info "  Wavelength: $wavelength m"
    info "  Spot size: $spotSize m"
    info "  Absorption coefficient: $absorptionCoeff 1/m"

    # Calculate fluence
    if [ -n "$pulseEnergy" ] && [ -n "$spotSize" ]; then
        fluence=$(echo "scale=6; $pulseEnergy / (3.14159265359 * ($spotSize/2)^2)" | bc -l)
        info "  Calculated fluence: $fluence J/m² ($(echo "scale=2; $fluence / 10000" | bc -l) J/cm²)"
    fi

    # Validate pulse width range (10fs to 10ps)
    if [ -n "$pulseWidth" ]; then
        if (( $(echo "$pulseWidth < 1e-14" | bc -l) )); then
            warning "Pulse width $pulseWidth s < 10 fs may be unrealistic"
        elif (( $(echo "$pulseWidth > 1e-11" | bc -l) )); then
            warning "Pulse width $pulseWidth s > 10 ps exceeds femtosecond regime"
        else
            success "Pulse width within femtosecond laser range"
        fi
    fi
fi

# Validate phase change coefficients
section "8. PHASE CHANGE VALIDATION"

check
if [ -f "system/controlDict" ]; then
    # Check phase change activation time
    activationTime=$(grep -A1 "activationTime" system/controlDict | grep -E "^\s*\(\(" | sed 's/[() ;]//g')
    info "Phase change activation time: $activationTime"

    # Check if activation covers laser pulse
    laserStartTime=$(grep -E "^\s*laserStartTime\s+" system/controlDict | awk '{print $2}' | tr -d ';')
    laserEndTime=$(grep -E "^\s*laserEndTime\s+" system/controlDict | awk '{print $2}' | tr -d ';')

    info "Laser timing:"
    info "  Start: $laserStartTime s"
    info "  End: $laserEndTime s"

    if [ -n "$endTime" ] && [ -n "$laserEndTime" ]; then
        if (( $(echo "$endTime < $laserEndTime" | bc -l) )); then
            warning "Simulation endTime ($endTime) < laserEndTime ($laserEndTime)"
        else
            success "Simulation duration covers laser pulse period"
        fi
    fi
fi

# Validate solver settings
section "9. SOLVER SETTINGS VALIDATION"

check
if [ -f "system/fvSolution" ]; then
    # Check PIMPLE settings
    nOuterCorrectors=$(grep -A20 "^PIMPLE" system/fvSolution | grep -E "^\s*nOuterCorrectors\s+" | awk '{print $2}' | tr -d ';')
    nCorrectors=$(grep -A20 "^PIMPLE" system/fvSolution | grep -E "^\s*nCorrectors\s+" | awk '{print $2}' | tr -d ';')
    nNonOrthogonalCorrectors=$(grep -A20 "^PIMPLE" system/fvSolution | grep -E "^\s*nNonOrthogonalCorrectors\s+" | awk '{print $2}' | tr -d ';')

    info "PIMPLE settings:"
    info "  nOuterCorrectors: $nOuterCorrectors"
    info "  nCorrectors: $nCorrectors"
    info "  nNonOrthogonalCorrectors: $nNonOrthogonalCorrectors"

    if [ -n "$nOuterCorrectors" ]; then
        if [ "$nOuterCorrectors" -lt 2 ]; then
            warning "nOuterCorrectors = $nOuterCorrectors may be too low for coupled solver"
        else
            success "nOuterCorrectors adequate for pressure-velocity coupling"
        fi
    fi
fi

# Check discretization schemes
section "10. DISCRETIZATION SCHEMES VALIDATION"

check
if [ -f "system/fvSchemes" ]; then
    # Check time scheme
    ddtScheme=$(grep -A5 "^ddtSchemes" system/fvSchemes | grep "default" | awk '{print $2}' | tr -d ';')
    info "Time discretization: $ddtScheme"

    if [ "$ddtScheme" == "Euler" ]; then
        success "Using Euler (1st order) time scheme - stable for small time steps"
    fi

    # Check alpha schemes
    alphaDiv=$(grep "div(phi,alpha)" system/fvSchemes | awk '{print $3}' | tr -d ';')
    info "Alpha convection scheme: $alphaDiv"

    if [ "$alphaDiv" == "vanLeer" ]; then
        success "Using vanLeer scheme for alpha - good for interface capturing"
    fi
fi

# Validate boundary conditions consistency
section "11. BOUNDARY CONDITION CONSISTENCY"

check
if [ -d "0.orig" ]; then
    # Get list of boundaries from a field file
    if [ -f "0.orig/U" ]; then
        info "Checking boundary consistency across fields..."

        # Extract boundary names from U file
        boundaries=$(sed -n '/^boundaryField/,/^}/p' 0.orig/U | grep -E "^\s+[a-zA-Z]" | awk '{print $1}')

        for field in T Te Tl p p_rgh alpha.metal alpha.air; do
            if [ -f "0.orig/$field" ]; then
                for bc in $boundaries; do
                    if ! grep -q "$bc" "0.orig/$field"; then
                        warning "Boundary '$bc' missing in 0.orig/$field"
                    fi
                done
            fi
        done
        success "Boundary condition consistency check completed"
    fi
fi

# Check for common issues
section "12. COMMON ISSUES CHECK"

# Check if 0 directory exists (should be created from 0.orig)
check
if [ ! -d "0" ]; then
    warning "Directory '0' not found - it will be created from 0.orig on first run"
    info "Run: cp -r 0.orig 0"
fi

# Check for old time directories that might interfere
check
old_dirs=$(find . -maxdepth 1 -type d -name "[0-9]*" -o -name "[0-9]*.[0-9]*" | grep -v "0.orig" | wc -l)
if [ "$old_dirs" -gt 0 ]; then
    warning "Found $old_dirs old time directories - consider running ./Allclean"
fi

# Check write permissions
check
if [ -w "." ]; then
    success "Write permissions OK in case directory"
else
    error "No write permission in case directory"
fi

# Memory estimation
section "13. RESOURCE ESTIMATION"

if [ -f "system/blockMeshDict" ]; then
    # Extract cell counts
    nx=$(grep -oP 'hex.*?\(\s*\d+\s+\d+\s+\d+\s*\)' system/blockMeshDict | sed 's/[^0-9]/ /g' | awk '{print $1}' | head -1)
    ny=$(grep -oP 'hex.*?\(\s*\d+\s+\d+\s+\d+\s*\)' system/blockMeshDict | sed 's/[^0-9]/ /g' | awk '{print $2}' | head -1)
    nz=$(grep -oP 'hex.*?\(\s*\d+\s+\d+\s+\d+\s*\)' system/blockMeshDict | sed 's/[^0-9]/ /g' | awk '{print $3}' | head -1)

    if [ -n "$nx" ] && [ -n "$ny" ] && [ -n "$nz" ]; then
        # Sum all blocks
        total_cells=0
        while read -r cells; do
            x=$(echo "$cells" | awk '{print $1}')
            y=$(echo "$cells" | awk '{print $2}')
            z=$(echo "$cells" | awk '{print $3}')
            block_cells=$((x * y * z))
            total_cells=$((total_cells + block_cells))
            info "Block: ${x}x${y}x${z} = $block_cells cells"
        done < <(grep -oP 'hex.*?\(\s*\d+\s+\d+\s+\d+\s*\)' system/blockMeshDict | sed 's/[^0-9]/ /g')

        info "Estimated total cells: $total_cells"

        # Rough memory estimate (assuming ~10-15 fields, ~500 bytes per cell per field)
        mem_mb=$((total_cells * 15 * 500 / 1024 / 1024))
        info "Estimated memory usage: ~${mem_mb} MB"
    fi
fi

# Time step estimation
if [ -n "$deltaT" ] && [ -n "$endTime" ] && [ -n "$startTime" ]; then
    nSteps=$(echo "scale=0; ($endTime - $startTime) / $deltaT" | bc -l)
    info "Estimated number of time steps: $nSteps"

    if [ -n "$nSteps" ] && [ "$nSteps" -gt 1000000 ]; then
        warning "Large number of time steps ($nSteps) - simulation may take considerable time"
    fi
fi

# Final summary
section "PRE-SIMULATION CHECK SUMMARY"

echo ""
echo "Total checks performed: $CHECKS"
echo -e "${GREEN}Passed: $((CHECKS - ERRORS - WARNINGS))${NC}"
echo -e "${YELLOW}Warnings: $WARNINGS${NC}"
echo -e "${RED}Errors: $ERRORS${NC}"
echo ""

if [ $ERRORS -eq 0 ] && [ $WARNINGS -eq 0 ]; then
    echo -e "${GREEN}✓ All checks passed! Case is ready for simulation.${NC}"
    exit 0
elif [ $ERRORS -eq 0 ]; then
    echo -e "${YELLOW}⚠ Case has warnings but should run. Review warnings above.${NC}"
    exit 0
else
    echo -e "${RED}✗ Case has errors that must be fixed before running.${NC}"
    exit 1
fi
