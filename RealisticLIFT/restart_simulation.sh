#!/bin/bash
# Clean previous run and restart with corrected laser timing

# Resume the stopped job first
fg

# Or if that doesn't work, kill it
# killall compInterFoam

# Clean the case
./Allclean

# Reinitialize and run
# (Add your initialization commands here based on your Allrun script)
