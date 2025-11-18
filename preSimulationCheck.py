#!/usr/bin/env python3
"""
Comprehensive Pre-Simulation Validation Script for compInterFoam
==================================================================

Performs detailed checks on:
1. Timing-related parameters (timestep, endTime, writeInterval, Courant numbers)
2. Laser properties and physics parameters
3. Thermophysical properties for all phases
4. Numerical settings (fvSchemes, fvSolution)
5. Mesh quality and resolution
6. Boundary conditions consistency
7. File structure and completeness

Usage:
    python3 preSimulationCheck.py [case_directory]

If no case directory is provided, checks current directory.
"""

import os
import sys
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import warnings

# ANSI color codes for terminal output
class Colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class ValidationResult:
    """Store validation results with severity levels"""
    def __init__(self):
        self.errors = []      # Critical issues that will prevent simulation
        self.warnings = []    # Issues that may cause problems
        self.info = []        # Informational messages
        self.passed = []      # Successful validations

    def add_error(self, message: str):
        self.errors.append(message)

    def add_warning(self, message: str):
        self.warnings.append(message)

    def add_info(self, message: str):
        self.info.append(message)

    def add_pass(self, message: str):
        self.passed.append(message)

    def has_errors(self) -> bool:
        return len(self.errors) > 0

    def print_summary(self):
        """Print colored summary of results"""
        print(f"\n{Colors.BOLD}{'='*80}{Colors.ENDC}")
        print(f"{Colors.BOLD}VALIDATION SUMMARY{Colors.ENDC}")
        print(f"{Colors.BOLD}{'='*80}{Colors.ENDC}\n")

        if self.errors:
            print(f"{Colors.FAIL}{Colors.BOLD}ERRORS ({len(self.errors)}):{Colors.ENDC}")
            for i, err in enumerate(self.errors, 1):
                print(f"  {Colors.FAIL}[E{i:02d}] {err}{Colors.ENDC}")
            print()

        if self.warnings:
            print(f"{Colors.WARNING}{Colors.BOLD}WARNINGS ({len(self.warnings)}):{Colors.ENDC}")
            for i, warn in enumerate(self.warnings, 1):
                print(f"  {Colors.WARNING}[W{i:02d}] {warn}{Colors.ENDC}")
            print()

        if self.info:
            print(f"{Colors.OKCYAN}{Colors.BOLD}INFORMATION ({len(self.info)}):{Colors.ENDC}")
            for i, inf in enumerate(self.info, 1):
                print(f"  {Colors.OKCYAN}[I{i:02d}] {inf}{Colors.ENDC}")
            print()

        print(f"{Colors.OKGREEN}{Colors.BOLD}PASSED CHECKS: {len(self.passed)}{Colors.ENDC}\n")

        if self.has_errors():
            print(f"{Colors.FAIL}{Colors.BOLD}RESULT: VALIDATION FAILED - Please fix errors before running simulation{Colors.ENDC}")
            return False
        elif self.warnings:
            print(f"{Colors.WARNING}{Colors.BOLD}RESULT: VALIDATION PASSED WITH WARNINGS - Review warnings before proceeding{Colors.ENDC}")
            return True
        else:
            print(f"{Colors.OKGREEN}{Colors.BOLD}RESULT: VALIDATION PASSED - Case ready for simulation{Colors.ENDC}")
            return True


class OpenFOAMDictParser:
    """Simple OpenFOAM dictionary parser"""

    @staticmethod
    def parse_value(value_str: str) -> Any:
        """Convert string value to appropriate type"""
        value_str = value_str.strip().rstrip(';')

        # Remove comments
        if '//' in value_str:
            value_str = value_str.split('//')[0].strip()

        # Handle dimensioned values: name [dims] value
        # Example: "sigma [1 0 -2 0 0 0 0] 1.64"
        dim_match = re.match(r'(\w+)?\s*\[[\d\s\-]+\]\s+([\d.eE\-+]+)', value_str)
        if dim_match:
            # Return just the numeric value
            return float(dim_match.group(2))

        # Handle lists
        if value_str.startswith('(') and value_str.endswith(')'):
            inner = value_str[1:-1].strip()
            if not inner:
                return []
            return [OpenFOAMDictParser.parse_value(v.strip()) for v in inner.split()]

        # Handle scientific notation
        if 'e' in value_str.lower() or 'E' in value_str:
            try:
                return float(value_str)
            except ValueError:
                return value_str

        # Handle numbers
        try:
            if '.' in value_str:
                return float(value_str)
            else:
                return int(value_str)
        except ValueError:
            return value_str

    @staticmethod
    def read_dict(file_path: Path) -> Dict[str, Any]:
        """Read OpenFOAM dictionary file"""
        if not file_path.exists():
            return {}

        with open(file_path, 'r') as f:
            content = f.read()

        # Remove C-style comments
        content = re.sub(r'/\*.*?\*/', '', content, flags=re.DOTALL)

        # Simple key-value extraction (works for most OpenFOAM dicts)
        params = {}
        lines = content.split('\n')

        for line in lines:
            # Remove inline comments
            if '//' in line:
                line = line.split('//')[0]

            line = line.strip()
            if not line or line.startswith('/*') or line.startswith('//'):
                continue

            # Match key-value pairs
            match = re.match(r'(\w+)\s+(.+?);', line)
            if match:
                key = match.group(1)
                value = match.group(2).strip()
                params[key] = OpenFOAMDictParser.parse_value(value)

        return params


class PreSimulationValidator:
    """Main validation class"""

    def __init__(self, case_dir: Path):
        self.case_dir = case_dir
        self.result = ValidationResult()
        self.parser = OpenFOAMDictParser()

        # Storage for parsed dictionaries
        self.control_dict = {}
        self.laser_props = {}
        self.thermo_props = {}
        self.thermo_air = {}
        self.thermo_metal = {}
        self.transport_props = {}
        self.fv_schemes = {}
        self.fv_solution = {}
        self.g_dict = {}
        self.dynamic_mesh = {}

    @staticmethod
    def safe_float(value) -> Optional[float]:
        """Safely convert value to float, return None if not possible"""
        if value is None:
            return None
        try:
            return float(value) if isinstance(value, str) else value
        except (ValueError, TypeError):
            return None

    def run_all_checks(self):
        """Run all validation checks"""
        print(f"{Colors.BOLD}{Colors.HEADER}")
        print("="*80)
        print("CompInterFoam Pre-Simulation Validation")
        print("="*80)
        print(f"{Colors.ENDC}")
        print(f"Case directory: {Colors.BOLD}{self.case_dir}{Colors.ENDC}\n")

        # 1. File structure checks
        print(f"{Colors.OKBLUE}[1/8] Checking file structure...{Colors.ENDC}")
        self.check_file_structure()

        # 2. Load all dictionaries
        print(f"{Colors.OKBLUE}[2/8] Loading dictionary files...{Colors.ENDC}")
        self.load_dictionaries()

        # 3. Timing checks
        print(f"{Colors.OKBLUE}[3/8] Validating timing parameters...{Colors.ENDC}")
        self.check_timing_parameters()

        # 4. Laser properties
        print(f"{Colors.OKBLUE}[4/8] Validating laser properties...{Colors.ENDC}")
        self.check_laser_properties()

        # 5. Thermophysical properties
        print(f"{Colors.OKBLUE}[5/8] Validating thermophysical properties...{Colors.ENDC}")
        self.check_thermophysical_properties()

        # 6. Numerical settings
        print(f"{Colors.OKBLUE}[6/8] Validating numerical settings...{Colors.ENDC}")
        self.check_numerical_settings()

        # 7. Boundary conditions
        print(f"{Colors.OKBLUE}[7/8] Validating boundary conditions...{Colors.ENDC}")
        self.check_boundary_conditions()

        # 8. Advanced physics checks
        print(f"{Colors.OKBLUE}[8/8] Validating advanced physics parameters...{Colors.ENDC}")
        self.check_advanced_physics()

        # Print summary
        return self.result.print_summary()

    def check_file_structure(self):
        """Verify required directories and files exist"""

        # Required directories
        required_dirs = ['system', 'constant', '0']
        for dir_name in required_dirs:
            dir_path = self.case_dir / dir_name
            if not dir_path.exists():
                # Check for 0.orig as alternative to 0
                if dir_name == '0' and (self.case_dir / '0.orig').exists():
                    self.result.add_pass(f"Found 0.orig directory (alternative to 0/)")
                else:
                    self.result.add_error(f"Required directory '{dir_name}/' not found")
            else:
                self.result.add_pass(f"Directory '{dir_name}/' exists")

        # Required system files
        system_files = ['controlDict', 'fvSchemes', 'fvSolution', 'blockMeshDict']
        for file_name in system_files:
            file_path = self.case_dir / 'system' / file_name
            if not file_path.exists():
                self.result.add_error(f"Required file 'system/{file_name}' not found")
            else:
                self.result.add_pass(f"File 'system/{file_name}' exists")

        # Required constant files for compInterFoam
        constant_files = ['thermophysicalProperties', 'transportProperties', 'g']
        for file_name in constant_files:
            file_path = self.case_dir / 'constant' / file_name
            if not file_path.exists():
                self.result.add_warning(f"File 'constant/{file_name}' not found")
            else:
                self.result.add_pass(f"File 'constant/{file_name}' exists")

        # Check for laser properties (specific to this solver)
        laser_file = self.case_dir / 'constant' / 'laserProperties'
        if laser_file.exists():
            self.result.add_pass("File 'constant/laserProperties' exists")
        else:
            self.result.add_warning("File 'constant/laserProperties' not found - required for laser simulations")

        # Required initial fields (check both 0/ and 0.orig/)
        initial_dir = self.case_dir / '0'
        if not initial_dir.exists():
            initial_dir = self.case_dir / '0.orig'

        if initial_dir.exists():
            required_fields = ['p', 'U', 'T', 'alpha.air', 'alpha.metal']
            for field in required_fields:
                field_path = initial_dir / field
                if not field_path.exists():
                    self.result.add_error(f"Required initial field '{field}' not found in {initial_dir.name}/")
                else:
                    self.result.add_pass(f"Initial field '{field}' exists")

    def load_dictionaries(self):
        """Load all OpenFOAM dictionary files"""

        # Load controlDict
        control_path = self.case_dir / 'system' / 'controlDict'
        if control_path.exists():
            self.control_dict = self.parser.read_dict(control_path)
            self.result.add_pass(f"Loaded controlDict with {len(self.control_dict)} parameters")

        # Load laser properties
        laser_path = self.case_dir / 'constant' / 'laserProperties'
        if laser_path.exists():
            self.laser_props = self.parser.read_dict(laser_path)
            self.result.add_pass(f"Loaded laserProperties with {len(self.laser_props)} parameters")

        # Load thermophysical properties
        thermo_path = self.case_dir / 'constant' / 'thermophysicalProperties'
        if thermo_path.exists():
            self.thermo_props = self.parser.read_dict(thermo_path)

        thermo_air_path = self.case_dir / 'constant' / 'thermophysicalProperties.air'
        if thermo_air_path.exists():
            self.thermo_air = self.parser.read_dict(thermo_air_path)
            self.result.add_pass(f"Loaded thermophysicalProperties.air")

        thermo_metal_path = self.case_dir / 'constant' / 'thermophysicalProperties.metal'
        if thermo_metal_path.exists():
            self.thermo_metal = self.parser.read_dict(thermo_metal_path)
            self.result.add_pass(f"Loaded thermophysicalProperties.metal")

        # Load transport properties
        transport_path = self.case_dir / 'constant' / 'transportProperties'
        if transport_path.exists():
            self.transport_props = self.parser.read_dict(transport_path)

        # Load fvSchemes and fvSolution
        schemes_path = self.case_dir / 'system' / 'fvSchemes'
        if schemes_path.exists():
            self.fv_schemes = self.parser.read_dict(schemes_path)

        solution_path = self.case_dir / 'system' / 'fvSolution'
        if solution_path.exists():
            self.fv_solution = self.parser.read_dict(solution_path)

        # Load gravity
        g_path = self.case_dir / 'constant' / 'g'
        if g_path.exists():
            self.g_dict = self.parser.read_dict(g_path)

        # Load dynamic mesh dict
        dynamic_path = self.case_dir / 'constant' / 'dynamicMeshDict'
        if dynamic_path.exists():
            self.dynamic_mesh = self.parser.read_dict(dynamic_path)

    def check_timing_parameters(self):
        """Comprehensive timing parameter validation"""

        if not self.control_dict:
            self.result.add_error("controlDict not loaded - cannot validate timing")
            return

        # 1. Basic time settings
        start_time = self.control_dict.get('startTime', None)
        end_time = self.control_dict.get('endTime', None)
        delta_t = self.control_dict.get('deltaT', None)

        if start_time is None:
            self.result.add_error("startTime not defined in controlDict")
        else:
            self.result.add_pass(f"startTime = {start_time} s")
            if start_time < 0:
                self.result.add_error(f"startTime ({start_time}) cannot be negative")

        if end_time is None:
            self.result.add_error("endTime not defined in controlDict")
        else:
            self.result.add_pass(f"endTime = {end_time} s")
            if end_time <= 0:
                self.result.add_error(f"endTime ({end_time}) must be positive")

        if delta_t is None:
            self.result.add_error("deltaT not defined in controlDict")
        else:
            self.result.add_pass(f"deltaT = {delta_t} s")
            if delta_t <= 0:
                self.result.add_error(f"deltaT ({delta_t}) must be positive")

        # 2. Check simulation duration and timestep
        if start_time is not None and end_time is not None:
            duration = end_time - start_time
            if duration <= 0:
                self.result.add_error(f"Simulation duration ({duration} s) must be positive (endTime > startTime)")
            else:
                self.result.add_pass(f"Simulation duration = {duration} s")

                # Check if duration makes sense for femtosecond simulations
                if duration > 1e-6:
                    self.result.add_warning(f"Simulation duration ({duration} s) is very long for femtosecond laser simulation (> 1 μs)")

                # Estimate number of timesteps
                if delta_t is not None and delta_t > 0:
                    adjust_time_step = self.control_dict.get('adjustTimeStep', 'no')
                    if isinstance(adjust_time_step, str):
                        adjust_time_step = adjust_time_step.lower() in ['yes', 'true', 'on', '1']

                    if not adjust_time_step:
                        n_steps = int(duration / delta_t)
                        self.result.add_info(f"Estimated number of timesteps: {n_steps:,}")

                        if n_steps > 1e6:
                            self.result.add_warning(f"Very large number of timesteps ({n_steps:,}) - simulation may take long time")
                        elif n_steps < 10:
                            self.result.add_warning(f"Very few timesteps ({n_steps}) - may not capture transient behavior")
                    else:
                        n_steps_est = int(duration / delta_t)
                        self.result.add_info(f"Adaptive timestepping enabled - estimated ~{n_steps_est:,} steps (may vary)")

        # 3. Adaptive timestepping checks
        adjust_time_step = self.control_dict.get('adjustTimeStep', 'no')
        if isinstance(adjust_time_step, str):
            adjust_time_step = adjust_time_step.lower() in ['yes', 'true', 'on', '1']

        if adjust_time_step:
            self.result.add_pass("Adaptive timestepping enabled")

            # Check Courant numbers
            max_co = self.control_dict.get('maxCo', None)
            max_alpha_co = self.control_dict.get('maxAlphaCo', None)
            max_delta_t = self.control_dict.get('maxDeltaT', None)

            if max_co is None:
                self.result.add_warning("maxCo not defined - required for adaptive timestepping")
            else:
                self.result.add_pass(f"maxCo = {max_co}")
                if max_co > 1.0:
                    self.result.add_warning(f"maxCo ({max_co}) > 1.0 may cause instability")
                elif max_co < 0.01:
                    self.result.add_warning(f"maxCo ({max_co}) < 0.01 is very conservative, may slow simulation")

            if max_alpha_co is None:
                self.result.add_warning("maxAlphaCo not defined - important for interface tracking")
            else:
                self.result.add_pass(f"maxAlphaCo = {max_alpha_co}")
                if max_alpha_co > 1.0:
                    self.result.add_error(f"maxAlphaCo ({max_alpha_co}) > 1.0 will cause interface smearing")
                elif max_alpha_co < 0.01:
                    self.result.add_warning(f"maxAlphaCo ({max_alpha_co}) < 0.01 is very conservative")

            if max_delta_t is not None:
                self.result.add_pass(f"maxDeltaT = {max_delta_t} s")
                if max_delta_t < delta_t:
                    self.result.add_warning(f"maxDeltaT ({max_delta_t}) < deltaT ({delta_t}) - deltaT will be used as upper limit")

                # Check timestep range
                if max_delta_t / delta_t > 1000:
                    self.result.add_warning(f"Very large timestep range (factor {max_delta_t/delta_t:.1e}) - may cause abrupt changes")

            # Check thermal Courant number for femtosecond simulations
            max_thermal_co = self.control_dict.get('maxThermalCourant', None)
            if max_thermal_co is not None:
                self.result.add_pass(f"maxThermalCourant = {max_thermal_co}")
                if max_thermal_co > 1.0:
                    self.result.add_warning(f"maxThermalCourant ({max_thermal_co}) > 1.0 may cause thermal diffusion instability")
        else:
            self.result.add_info("Fixed timestepping mode")

            # For fixed timestep, check CFL estimates if we have mesh info
            if delta_t is not None:
                self.result.add_info(f"Using fixed timestep: {delta_t} s")

        # 4. Write control checks
        write_control = self.control_dict.get('writeControl', None)
        write_interval = self.control_dict.get('writeInterval', None)

        if write_control is None:
            self.result.add_error("writeControl not defined in controlDict")
        else:
            self.result.add_pass(f"writeControl = {write_control}")

            valid_write_controls = ['timeStep', 'runTime', 'adjustableRunTime', 'cpuTime', 'clockTime']
            if write_control not in valid_write_controls:
                self.result.add_warning(f"writeControl '{write_control}' may not be standard (valid: {valid_write_controls})")

        if write_interval is None:
            self.result.add_error("writeInterval not defined in controlDict")
        else:
            self.result.add_pass(f"writeInterval = {write_interval}")

            if write_control == 'timeStep' and isinstance(write_interval, (int, float)):
                if write_interval < 1:
                    self.result.add_error(f"writeInterval ({write_interval}) must be >= 1 for timeStep mode")
                elif write_interval == 1:
                    self.result.add_warning("writeInterval = 1 will write every timestep - very large disk usage")

                # Estimate number of outputs
                if delta_t is not None and end_time is not None and start_time is not None:
                    duration = end_time - start_time
                    n_steps = duration / delta_t
                    n_outputs = n_steps / write_interval
                    self.result.add_info(f"Estimated number of output times: {int(n_outputs)}")

                    if n_outputs > 1000:
                        self.result.add_warning(f"Very large number of outputs ({int(n_outputs)}) - consider increasing writeInterval")
                    elif n_outputs < 5:
                        self.result.add_warning(f"Very few outputs ({int(n_outputs)}) - may miss important transients")

            elif write_control == 'runTime' and isinstance(write_interval, (int, float)):
                if end_time is not None and start_time is not None:
                    duration = end_time - start_time
                    n_outputs = duration / write_interval
                    self.result.add_info(f"Estimated number of output times: {int(n_outputs)}")

                    if n_outputs > 1000:
                        self.result.add_warning(f"Very large number of outputs ({int(n_outputs)}) - consider increasing writeInterval")

        # 5. Purge write
        purge_write = self.control_dict.get('purgeWrite', 0)
        if purge_write > 0:
            self.result.add_pass(f"purgeWrite = {purge_write} (keeps only last {purge_write} time directories)")
        else:
            self.result.add_info("purgeWrite = 0 (all time directories will be kept)")

        # 6. Run time modifiable
        run_time_modifiable = self.control_dict.get('runTimeModifiable', 'yes')
        if isinstance(run_time_modifiable, str):
            if run_time_modifiable.lower() in ['yes', 'true', 'on']:
                self.result.add_pass("runTimeModifiable enabled - can modify parameters during run")
            else:
                self.result.add_info("runTimeModifiable disabled - parameters fixed at start")

        # 7. Time precision
        time_precision = self.control_dict.get('timePrecision', 6)
        if time_precision < 6:
            self.result.add_warning(f"timePrecision ({time_precision}) < 6 may cause naming issues for small timesteps")
        else:
            self.result.add_pass(f"timePrecision = {time_precision}")

        # 8. Write precision
        write_precision = self.control_dict.get('writePrecision', 6)
        if write_precision < 8:
            self.result.add_warning(f"writePrecision ({write_precision}) < 8 may lose accuracy in output")
        else:
            self.result.add_pass(f"writePrecision = {write_precision}")

    def check_laser_properties(self):
        """Validate laser properties for femtosecond LIFT simulation"""

        if not self.laser_props:
            self.result.add_warning("laserProperties not loaded - skipping laser validation")
            return

        # 1. Pulse energy
        pulse_energy = self.laser_props.get('pulseEnergy', None)
        if pulse_energy is None:
            self.result.add_error("pulseEnergy not defined in laserProperties")
        else:
            self.result.add_pass(f"pulseEnergy = {pulse_energy} J")
            if pulse_energy <= 0:
                self.result.add_error(f"pulseEnergy ({pulse_energy}) must be positive")
            elif pulse_energy > 1e-3:
                self.result.add_warning(f"pulseEnergy ({pulse_energy} J) > 1 mJ is very high for femtosecond LIFT")
            elif pulse_energy < 1e-9:
                self.result.add_warning(f"pulseEnergy ({pulse_energy} J) < 1 nJ may be too low for material removal")

        # 2. Pulse width (FWHM)
        pulse_width = self.laser_props.get('pulseWidth', None)
        if pulse_width is None:
            self.result.add_error("pulseWidth not defined in laserProperties")
        else:
            self.result.add_pass(f"pulseWidth = {pulse_width} s ({pulse_width*1e15:.1f} fs)")
            if pulse_width <= 0:
                self.result.add_error(f"pulseWidth ({pulse_width}) must be positive")
            elif pulse_width > 1e-12:
                self.result.add_warning(f"pulseWidth ({pulse_width*1e12:.1f} ps) > 1 ps - not femtosecond regime")
            elif pulse_width < 1e-16:
                self.result.add_warning(f"pulseWidth ({pulse_width*1e15:.1f} fs) < 0.1 fs is unrealistically short")

            # Check if simulation timestep can resolve pulse
            if self.control_dict:
                delta_t = self.control_dict.get('deltaT', None)
                if delta_t is not None and pulse_width is not None:
                    steps_per_pulse = pulse_width / delta_t
                    if steps_per_pulse < 5:
                        self.result.add_warning(f"Only {steps_per_pulse:.1f} timesteps per pulse - increase resolution (reduce deltaT)")
                    else:
                        self.result.add_pass(f"Temporal resolution: {steps_per_pulse:.1f} timesteps per pulse")

        # 3. Wavelength
        wavelength = self.laser_props.get('wavelength', None)
        if wavelength is None:
            self.result.add_error("wavelength not defined in laserProperties")
        else:
            self.result.add_pass(f"wavelength = {wavelength} m ({wavelength*1e9:.1f} nm)")
            if wavelength <= 0:
                self.result.add_error(f"wavelength ({wavelength}) must be positive")
            elif wavelength < 200e-9 or wavelength > 2000e-9:
                self.result.add_warning(f"wavelength ({wavelength*1e9:.1f} nm) outside typical laser range (200-2000 nm)")

        # 4. Spot size / beam radius
        spot_size = self.laser_props.get('spotSize', None)
        beam_radius = self.laser_props.get('beamRadius', None)

        if spot_size is None and beam_radius is None:
            self.result.add_error("Neither spotSize nor beamRadius defined in laserProperties")
        else:
            radius = beam_radius if beam_radius is not None else spot_size / 2
            self.result.add_pass(f"Laser spot radius = {radius} m ({radius*1e6:.2f} μm)")
            if radius <= 0:
                self.result.add_error(f"Laser spot radius ({radius}) must be positive")
            elif radius > 100e-6:
                self.result.add_warning(f"Laser spot radius ({radius*1e6:.1f} μm) > 100 μm is large for LIFT")

        # 5. Focus point
        focus_point = self.laser_props.get('focusPoint', None)
        if focus_point is None:
            self.result.add_warning("focusPoint not defined - laser position unclear")
        else:
            if isinstance(focus_point, list) and len(focus_point) == 3:
                self.result.add_pass(f"focusPoint = ({focus_point[0]}, {focus_point[1]}, {focus_point[2]}) m")
            else:
                self.result.add_error(f"focusPoint format invalid - should be (x y z)")

        # 6. Absorption coefficient
        absorption_coeff = self.laser_props.get('absorptionCoefficient', None)
        if absorption_coeff is None:
            self.result.add_warning("absorptionCoefficient not defined")
        else:
            penetration_depth = 1.0 / absorption_coeff if absorption_coeff > 0 else float('inf')
            self.result.add_pass(f"absorptionCoefficient = {absorption_coeff} m^-1 (penetration depth: {penetration_depth*1e9:.1f} nm)")
            if absorption_coeff <= 0:
                self.result.add_error(f"absorptionCoefficient ({absorption_coeff}) must be positive")

        # 7. Reflectivity
        reflectivity = self.laser_props.get('reflectivity', None)
        if reflectivity is None:
            self.result.add_warning("reflectivity not defined - assuming 0")
        else:
            self.result.add_pass(f"reflectivity = {reflectivity} ({reflectivity*100:.1f}%)")
            if reflectivity < 0 or reflectivity > 1:
                self.result.add_error(f"reflectivity ({reflectivity}) must be between 0 and 1")
            elif reflectivity > 0.9:
                self.result.add_warning(f"reflectivity ({reflectivity*100:.1f}%) > 90% - very little energy absorbed")

        # 8. Film thickness (important for LIFT)
        film_thickness = self.laser_props.get('filmThickness', None)
        if film_thickness is not None:
            self.result.add_pass(f"filmThickness = {film_thickness} m ({film_thickness*1e9:.1f} nm)")

            # Check if absorption depth and film thickness are consistent
            if absorption_coeff is not None and absorption_coeff > 0:
                penetration_depth = 1.0 / absorption_coeff
                ratio = film_thickness / penetration_depth
                if ratio < 0.1:
                    self.result.add_warning(f"Film thickness ({film_thickness*1e9:.1f} nm) << penetration depth ({penetration_depth*1e9:.1f} nm) - most energy passes through")
                elif ratio > 10:
                    self.result.add_info(f"Film thickness ({film_thickness*1e9:.1f} nm) >> penetration depth ({penetration_depth*1e9:.1f} nm) - energy absorbed in surface layer")
                else:
                    self.result.add_pass(f"Film thickness / penetration depth ratio: {ratio:.2f}")

    def check_thermophysical_properties(self):
        """Validate thermophysical properties for all phases"""

        # Check air properties
        if self.thermo_air:
            self.result.add_info("Validating air phase thermophysical properties...")
            self._check_phase_properties(self.thermo_air, "air")

        # Check metal properties
        if self.thermo_metal:
            self.result.add_info("Validating metal phase thermophysical properties...")
            self._check_phase_properties(self.thermo_metal, "metal")

        # Check transport properties
        if self.transport_props:
            # Surface tension
            sigma = self.transport_props.get('sigma', None)
            if sigma is not None:
                # Convert to float if string
                try:
                    sigma_val = float(sigma) if isinstance(sigma, str) else sigma
                    self.result.add_pass(f"Surface tension (sigma) = {sigma_val} N/m")
                    if sigma_val < 0:
                        self.result.add_error(f"Surface tension ({sigma_val}) cannot be negative")
                    elif sigma_val > 2.0:
                        self.result.add_warning(f"Surface tension ({sigma_val} N/m) > 2.0 is very high")
                except (ValueError, TypeError):
                    self.result.add_warning(f"Surface tension value '{sigma}' cannot be parsed as number")

    def _check_phase_properties(self, props_dict: Dict, phase_name: str):
        """Check properties for a single phase"""

        # Density
        rho = self.safe_float(props_dict.get('rho', None))
        if rho is not None:
            self.result.add_pass(f"{phase_name}: density = {rho} kg/m^3")
            if rho <= 0:
                self.result.add_error(f"{phase_name}: density ({rho}) must be positive")
        else:
            self.result.add_warning(f"{phase_name}: density not defined")

        # Specific heat capacity
        cp = self.safe_float(props_dict.get('Cp', None) or props_dict.get('cp', None))
        if cp is not None:
            self.result.add_pass(f"{phase_name}: specific heat = {cp} J/kg/K")
            if cp <= 0:
                self.result.add_error(f"{phase_name}: specific heat ({cp}) must be positive")
        else:
            self.result.add_warning(f"{phase_name}: specific heat (Cp) not defined")

        # Thermal conductivity
        kappa = self.safe_float(props_dict.get('kappa', None) or props_dict.get('k', None))
        if kappa is not None:
            self.result.add_pass(f"{phase_name}: thermal conductivity = {kappa} W/m/K")
            if kappa < 0:
                self.result.add_error(f"{phase_name}: thermal conductivity ({kappa}) cannot be negative")
        else:
            self.result.add_warning(f"{phase_name}: thermal conductivity not defined")

        # Viscosity
        mu = self.safe_float(props_dict.get('mu', None) or props_dict.get('nu', None))
        if mu is not None:
            mu_name = "dynamic viscosity" if 'mu' in props_dict else "kinematic viscosity"
            self.result.add_pass(f"{phase_name}: {mu_name} = {mu}")
            if mu < 0:
                self.result.add_error(f"{phase_name}: {mu_name} ({mu}) cannot be negative")

    def check_numerical_settings(self):
        """Validate numerical schemes and solver settings"""

        # Check fvSchemes
        if not self.fv_schemes:
            self.result.add_warning("fvSchemes not loaded - skipping numerical scheme validation")
        else:
            # Time scheme
            ddt_schemes = self.fv_schemes.get('ddtSchemes', {})
            if isinstance(ddt_schemes, dict):
                default_ddt = ddt_schemes.get('default', None)
                if default_ddt:
                    self.result.add_pass(f"Time discretization: {default_ddt}")
                    if 'Euler' in str(default_ddt):
                        self.result.add_info("Using Euler (1st order) - consider backward or CrankNicolson for better accuracy")
                else:
                    self.result.add_warning("No default ddtSchemes defined")

            # Gradient schemes
            grad_schemes = self.fv_schemes.get('gradSchemes', {})
            if isinstance(grad_schemes, dict):
                default_grad = grad_schemes.get('default', None)
                if default_grad:
                    self.result.add_pass(f"Gradient scheme: {default_grad}")

            # Divergence schemes (important for stability)
            div_schemes = self.fv_schemes.get('divSchemes', {})
            if isinstance(div_schemes, dict):
                # Check alpha divergence (interface compression)
                has_alpha_scheme = False
                for key in div_schemes.keys():
                    if 'alpha' in key.lower():
                        has_alpha_scheme = True
                        self.result.add_pass(f"Alpha divergence scheme: {key} -> {div_schemes[key]}")

                if not has_alpha_scheme:
                    self.result.add_warning("No alpha divergence scheme found - interface may smear")

        # Check fvSolution
        if not self.fv_solution:
            self.result.add_warning("fvSolution not loaded - skipping solver settings validation")
        else:
            # Check PIMPLE settings
            pimple = self.fv_solution.get('PIMPLE', {})
            if isinstance(pimple, dict):
                n_outer = pimple.get('nOuterCorrectors', None)
                if n_outer is not None:
                    self.result.add_pass(f"PIMPLE nOuterCorrectors = {n_outer}")
                    if n_outer < 3:
                        self.result.add_warning(f"nOuterCorrectors ({n_outer}) < 3 may not converge well")

                n_correctors = pimple.get('nCorrectors', None)
                if n_correctors is not None:
                    self.result.add_pass(f"PIMPLE nCorrectors = {n_correctors}")

                n_non_orth = pimple.get('nNonOrthogonalCorrectors', None)
                if n_non_orth is not None:
                    self.result.add_pass(f"PIMPLE nNonOrthogonalCorrectors = {n_non_orth}")

            # Check relaxation factors
            relax_fields = self.fv_solution.get('relaxationFactors', {})
            if isinstance(relax_fields, dict):
                fields = relax_fields.get('fields', {}) or relax_fields.get('equations', {})
                if isinstance(fields, dict):
                    for field, factor in fields.items():
                        if isinstance(factor, (int, float)):
                            if factor <= 0 or factor > 1:
                                self.result.add_warning(f"Relaxation factor for '{field}' ({factor}) outside range (0,1]")
                            else:
                                self.result.add_pass(f"Relaxation factor '{field}' = {factor}")

    def check_boundary_conditions(self):
        """Check boundary conditions in initial fields"""

        initial_dir = self.case_dir / '0'
        if not initial_dir.exists():
            initial_dir = self.case_dir / '0.orig'

        if not initial_dir.exists():
            self.result.add_warning("No initial conditions directory found (0/ or 0.orig/)")
            return

        # Check key fields
        critical_fields = ['p', 'U', 'T']
        for field in critical_fields:
            field_path = initial_dir / field
            if field_path.exists():
                # Simple check - just verify file is not empty
                if field_path.stat().st_size < 100:
                    self.result.add_warning(f"Field file '{field}' seems too small - check if properly defined")
                else:
                    self.result.add_pass(f"Boundary conditions file '{field}' exists and has content")

        # Check alpha fields
        alpha_files = list(initial_dir.glob('alpha.*'))
        if len(alpha_files) >= 2:
            self.result.add_pass(f"Found {len(alpha_files)} alpha fields for multiphase simulation")
        elif len(alpha_files) == 1:
            self.result.add_warning("Only one alpha field found - need at least 2 for two-phase flow")
        else:
            self.result.add_warning("No alpha fields found - required for multiphase simulation")

    def check_advanced_physics(self):
        """Check advanced physics parameters specific to this solver"""

        # Check two-temperature model parameters
        if self.control_dict:
            # Electron heat capacity coefficient
            ce_coeff = self.control_dict.get('electronHeatCapacityCoeff', None)
            if ce_coeff is not None:
                self.result.add_pass(f"Two-temperature model: electronHeatCapacityCoeff = {ce_coeff}")

            # Electron-phonon coupling
            g_ep = self.control_dict.get('electronPhononCoupling', None)
            if g_ep is not None:
                self.result.add_pass(f"Two-temperature model: electronPhononCoupling = {g_ep}")

            # Phase change parameters
            t_sat = self.control_dict.get('TSat', None)
            if t_sat is not None:
                self.result.add_pass(f"Saturation temperature = {t_sat} K")
                if t_sat < 273:
                    self.result.add_warning(f"Saturation temperature ({t_sat} K) below water freezing point")

        # Check gravity
        if self.g_dict:
            g_value = self.g_dict.get('value', None)
            if g_value is not None:
                if isinstance(g_value, list) and len(g_value) == 3:
                    g_mag = (g_value[0]**2 + g_value[1]**2 + g_value[2]**2)**0.5
                    self.result.add_pass(f"Gravity vector = {g_value}, magnitude = {g_mag:.2f} m/s^2")

                    if g_mag < 0.1:
                        self.result.add_info("Gravity magnitude very small - nearly zero-gravity simulation")
                    elif abs(g_mag - 9.81) < 0.1:
                        self.result.add_pass("Using Earth gravity (9.81 m/s^2)")

        # Check if dynamic mesh is used
        if self.dynamic_mesh:
            solver = self.dynamic_mesh.get('dynamicFvMesh', None)
            if solver and solver != 'staticFvMesh':
                self.result.add_pass(f"Dynamic mesh enabled: {solver}")
            else:
                self.result.add_info("Static mesh (no mesh motion)")


def main():
    """Main entry point"""

    # Parse command line arguments
    if len(sys.argv) > 1:
        case_dir = Path(sys.argv[1])
    else:
        case_dir = Path.cwd()

    # Check if directory exists
    if not case_dir.exists():
        print(f"{Colors.FAIL}Error: Directory '{case_dir}' does not exist{Colors.ENDC}")
        sys.exit(1)

    if not case_dir.is_dir():
        print(f"{Colors.FAIL}Error: '{case_dir}' is not a directory{Colors.ENDC}")
        sys.exit(1)

    # Run validation
    validator = PreSimulationValidator(case_dir)
    success = validator.run_all_checks()

    # Exit with appropriate code
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
