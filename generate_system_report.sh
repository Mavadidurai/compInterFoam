#!/bin/bash
#------------------------------------------------------------------------------
# System Information Report Generator for OpenFOAM
#
# This script generates a comprehensive system information report including:
# - Hardware specifications (CPU, RAM, GPU)
# - Operating system details
# - OpenFOAM environment information
# - Installed dependencies and libraries
#
# Usage: ./generate_system_report.sh [output_file]
#------------------------------------------------------------------------------

# Set output file (default: system_info_report.txt)
OUTPUT_FILE="${1:-system_info_report.txt}"

# Function to print section header
print_header() {
    echo ""
    echo "================================================================================"
    echo " $1"
    echo "================================================================================"
    echo ""
}

# Start report
{
    print_header "System Information Report - $(date)"

    # -------------------------------------------------------------------------
    # System Overview
    # -------------------------------------------------------------------------
    print_header "SYSTEM OVERVIEW"

    echo "Hostname: $(hostname)"
    echo "User: $(whoami)"
    echo "Uptime: $(uptime -p 2>/dev/null || uptime)"
    echo ""

    # -------------------------------------------------------------------------
    # Operating System
    # -------------------------------------------------------------------------
    print_header "OPERATING SYSTEM"

    if [ -f /etc/os-release ]; then
        cat /etc/os-release
    fi
    echo ""
    echo "Kernel: $(uname -r)"
    echo "Architecture: $(uname -m)"
    echo ""

    # -------------------------------------------------------------------------
    # CPU Information
    # -------------------------------------------------------------------------
    print_header "CPU INFORMATION"

    if [ -f /proc/cpuinfo ]; then
        echo "Model: $(grep 'model name' /proc/cpuinfo | head -1 | cut -d':' -f2 | xargs)"
        echo "Physical Cores: $(grep 'physical id' /proc/cpuinfo | sort -u | wc -l)"
        echo "Logical Cores: $(grep 'processor' /proc/cpuinfo | wc -l)"
        echo "CPU MHz: $(grep 'cpu MHz' /proc/cpuinfo | head -1 | cut -d':' -f2 | xargs)"
        echo ""
        echo "CPU Flags (first processor):"
        grep 'flags' /proc/cpuinfo | head -1 | cut -d':' -f2
    fi
    echo ""

    # -------------------------------------------------------------------------
    # Memory Information
    # -------------------------------------------------------------------------
    print_header "MEMORY INFORMATION"

    if [ -f /proc/meminfo ]; then
        echo "Total Memory: $(grep MemTotal /proc/meminfo | awk '{print $2/1024/1024 " GB"}')"
        echo "Free Memory: $(grep MemFree /proc/meminfo | awk '{print $2/1024/1024 " GB"}')"
        echo "Available Memory: $(grep MemAvailable /proc/meminfo | awk '{print $2/1024/1024 " GB"}')"
        echo "Swap Total: $(grep SwapTotal /proc/meminfo | awk '{print $2/1024/1024 " GB"}')"
        echo "Swap Free: $(grep SwapFree /proc/meminfo | awk '{print $2/1024/1024 " GB"}')"
    fi
    echo ""

    # -------------------------------------------------------------------------
    # GPU Information
    # -------------------------------------------------------------------------
    print_header "GPU INFORMATION"

    if command -v nvidia-smi &> /dev/null; then
        nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader
    else
        lspci | grep -i vga
    fi
    echo ""

    # -------------------------------------------------------------------------
    # Disk Information
    # -------------------------------------------------------------------------
    print_header "DISK INFORMATION"

    df -h | grep -E '^(/dev/|Filesystem)'
    echo ""

    # -------------------------------------------------------------------------
    # OpenFOAM Environment
    # -------------------------------------------------------------------------
    print_header "OPENFOAM ENVIRONMENT"

    if [ -n "$WM_PROJECT" ]; then
        echo "OpenFOAM is loaded"
        echo "  Project: $WM_PROJECT"
        echo "  Version: $WM_PROJECT_VERSION"
        echo "  Architecture: $WM_ARCH"
        echo "  Compiler: $WM_COMPILER"
        echo "  Compile Option: $WM_COMPILE_OPTION"
        echo "  Precision: $WM_PRECISION_OPTION"
        echo "  Label Size: $WM_LABEL_SIZE"
        echo "  MPI: $WM_MPLIB"
        echo ""
        echo "  Installation Dir: $WM_PROJECT_DIR"
        echo "  User Dir: $WM_PROJECT_USER_DIR"
        echo ""
    else
        echo "OpenFOAM environment not loaded"
        echo ""
        echo "To load OpenFOAM, source the appropriate bashrc file:"
        echo "  source /opt/openfoam/etc/bashrc"
        echo "  or"
        echo "  source ~/OpenFOAM/OpenFOAM-v2406/etc/bashrc"
        echo ""
    fi

    # -------------------------------------------------------------------------
    # OpenFOAM Executables
    # -------------------------------------------------------------------------
    print_header "OPENFOAM EXECUTABLES"

    for cmd in blockMesh decomposePar reconstructPar snappyHexMesh checkMesh; do
        if command -v $cmd &> /dev/null; then
            echo "$cmd: $(which $cmd)"
        else
            echo "$cmd: NOT FOUND"
        fi
    done
    echo ""

    # -------------------------------------------------------------------------
    # MPI Information
    # -------------------------------------------------------------------------
    print_header "MPI INFORMATION"

    if command -v mpirun &> /dev/null; then
        echo "mpirun: $(which mpirun)"
        mpirun --version 2>&1 | head -3
    else
        echo "mpirun: NOT FOUND"
    fi
    echo ""

    # -------------------------------------------------------------------------
    # Compiler Information
    # -------------------------------------------------------------------------
    print_header "COMPILER INFORMATION"

    if command -v gcc &> /dev/null; then
        echo "GCC Version:"
        gcc --version | head -1
        echo ""
    fi

    if command -v g++ &> /dev/null; then
        echo "G++ Version:"
        g++ --version | head -1
        echo ""
    fi

    if command -v gfortran &> /dev/null; then
        echo "GFortran Version:"
        gfortran --version | head -1
        echo ""
    fi

    # -------------------------------------------------------------------------
    # Python Information
    # -------------------------------------------------------------------------
    print_header "PYTHON INFORMATION"

    if command -v python3 &> /dev/null; then
        echo "Python3: $(which python3)"
        python3 --version
        echo ""
        echo "Key Python packages:"
        python3 -c "import sys; print('  numpy:', end=' '); import numpy; print(numpy.__version__)" 2>/dev/null || echo "  numpy: NOT FOUND"
        python3 -c "import sys; print('  matplotlib:', end=' '); import matplotlib; print(matplotlib.__version__)" 2>/dev/null || echo "  matplotlib: NOT FOUND"
        python3 -c "import sys; print('  scipy:', end=' '); import scipy; print(scipy.__version__)" 2>/dev/null || echo "  scipy: NOT FOUND"
    else
        echo "Python3: NOT FOUND"
    fi
    echo ""

    # -------------------------------------------------------------------------
    # ParaView Information
    # -------------------------------------------------------------------------
    print_header "PARAVIEW INFORMATION"

    if command -v paraview &> /dev/null; then
        echo "ParaView: $(which paraview)"
        paraview --version 2>&1 | head -1
    elif command -v paraFoam &> /dev/null; then
        echo "paraFoam: $(which paraFoam)"
    else
        echo "ParaView/paraFoam: NOT FOUND"
    fi
    echo ""

    # -------------------------------------------------------------------------
    # Mesh Utilities
    # -------------------------------------------------------------------------
    print_header "THIRD-PARTY MESH UTILITIES"

    for util in gmsh salome cfMesh; do
        if command -v $util &> /dev/null; then
            echo "$util: $(which $util)"
        else
            echo "$util: NOT FOUND"
        fi
    done
    echo ""

    # -------------------------------------------------------------------------
    # Key Libraries
    # -------------------------------------------------------------------------
    print_header "KEY LIBRARIES"

    echo "Checking for important libraries:"
    ldconfig -p 2>/dev/null | grep -E 'libmpi|libscotch|libmetis|libcgal|libboost' | head -20 || echo "  Unable to query ldconfig"
    echo ""

    # -------------------------------------------------------------------------
    # Environment Variables (OpenFOAM related)
    # -------------------------------------------------------------------------
    print_header "OPENFOAM-RELATED ENVIRONMENT VARIABLES"

    env | grep -E '^(WM_|FOAM_|MPI_)' | sort
    echo ""

    # -------------------------------------------------------------------------
    # End of Report
    # -------------------------------------------------------------------------
    print_header "END OF REPORT"

    echo "Report generated: $(date)"
    echo "Output file: $OUTPUT_FILE"

} | tee "$OUTPUT_FILE"

echo ""
echo "System information report saved to: $OUTPUT_FILE"
