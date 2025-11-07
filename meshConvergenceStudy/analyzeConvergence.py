#!/usr/bin/env python3
"""
Mesh Convergence Study - Analysis Script
========================================
Analyzes simulation results and computes Grid Convergence Index (GCI)
using Richardson extrapolation method.

Author: Generated for compInterFoam mesh convergence study
Date: 2025-11-07

References:
- Roache, P.J. (1994). Perspective: A Method for Uniform Reporting of
  Grid Refinement Studies. Journal of Fluids Engineering, 116(3), 405-413.
- ASME V&V 20-2009: Standard for Verification and Validation in
  Computational Fluid Dynamics and Heat Transfer
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from pathlib import Path
import re


class ConvergenceAnalyzer:
    """Analyzes mesh convergence using GCI methodology"""

    def __init__(self, study_dir="meshStudy"):
        self.study_dir = Path(study_dir)
        self.results = {}

        # Refinement factors between meshes
        self.refinement_ratios = {
            ('coarse', 'medium'): 2.0,
            ('medium', 'fine'): 1.5,
            ('fine', 'very_fine'): 1.33,
        }

    def extract_mesh_metrics(self, case_dir):
        """Extract mesh statistics from checkMesh log"""

        log_file = case_dir / "log.checkMesh"

        if not log_file.exists():
            print(f"  WARNING: {log_file} not found")
            return None

        metrics = {}

        with open(log_file, 'r') as f:
            content = f.read()

        # Extract total cells
        match = re.search(r'cells:\s*(\d+)', content)
        if match:
            metrics['total_cells'] = int(match.group(1))

        # Extract points
        match = re.search(r'points:\s*(\d+)', content)
        if match:
            metrics['total_points'] = int(match.group(1))

        # Extract mesh quality metrics
        match = re.search(r'max.*non-orthogonality.*?(\d+\.?\d*)', content, re.IGNORECASE)
        if match:
            metrics['max_non_orthogonality'] = float(match.group(1))

        match = re.search(r'max.*skewness.*?(\d+\.?\d*)', content, re.IGNORECASE)
        if match:
            metrics['max_skewness'] = float(match.group(1))

        # Calculate characteristic cell size (cubic root of domain volume / cell count)
        if 'total_cells' in metrics:
            # Domain: 50×32.0714×10 µm³ = 1.60357e-11 m³
            domain_volume = 50e-6 * 32.0714e-6 * 10e-6
            metrics['h'] = (domain_volume / metrics['total_cells']) ** (1/3)

        return metrics

    def extract_solution_metrics(self, case_dir, level_name):
        """Extract key solution metrics from simulation results"""

        metrics = {}

        # Find latest time directory
        time_dirs = []
        for item in case_dir.iterdir():
            if item.is_dir() and item.name.replace('.', '').replace('e-', '').replace('+', '').isdigit():
                try:
                    time_dirs.append((float(item.name), item))
                except ValueError:
                    pass

        if not time_dirs:
            print(f"  WARNING: No time directories found in {case_dir}")
            return None

        time_dirs.sort()
        latest_time = time_dirs[-1][1]
        metrics['latest_time'] = time_dirs[-1][0]

        print(f"  Analyzing results at t = {metrics['latest_time']:.3e} s")

        # Read field data
        fields_to_check = ['T', 'Te', 'Tl', 'U', 'p', 'alpha.metal']

        for field in fields_to_check:
            field_file = latest_time / field

            if field_file.exists():
                data = self._read_field_data(field_file)
                if data is not None:
                    metrics[f'{field}_max'] = np.max(data)
                    metrics[f'{field}_min'] = np.min(data)
                    metrics[f'{field}_mean'] = np.mean(data)
                    metrics[f'{field}_std'] = np.std(data)

                    # Special handling for temperature fields
                    if field in ['T', 'Te', 'Tl']:
                        # Peak temperature location (indicator of laser heating)
                        metrics[f'{field}_peak'] = np.max(data)

                    # Velocity magnitude
                    if field == 'U':
                        # U is a vector field
                        if len(data.shape) > 1 and data.shape[1] == 3:
                            u_mag = np.sqrt(np.sum(data**2, axis=1))
                            metrics['U_mag_max'] = np.max(u_mag)
                            metrics['U_mag_mean'] = np.mean(u_mag)

        return metrics

    def _read_field_data(self, field_file):
        """Read OpenFOAM field data (simplified parser)"""

        try:
            with open(field_file, 'r') as f:
                content = f.read()

            # Find internalField section
            match = re.search(r'internalField\s+(?:uniform|nonuniform)\s+(.*?);',
                             content, re.DOTALL)

            if not match:
                return None

            data_str = match.group(1).strip()

            # Handle uniform fields
            if 'uniform' in content.split('internalField')[1].split(';')[0]:
                # Uniform field - single value
                try:
                    # Try scalar
                    value = float(data_str.strip('()'))
                    return np.array([value])
                except:
                    # Try vector
                    values = re.findall(r'-?\d+\.?\d*(?:[eE][+-]?\d+)?', data_str)
                    if values:
                        return np.array([float(v) for v in values])

            # Handle nonuniform fields
            if data_str.startswith('List'):
                # Remove List wrapper
                data_str = re.sub(r'List\s*<\w+>\s*', '', data_str)

            # Extract numeric data
            # Try scalar field first
            scalars = re.findall(r'-?\d+\.?\d*(?:[eE][+-]?\d+)?', data_str)
            if scalars:
                data = np.array([float(v) for v in scalars])

                # Check if it's actually a vector field (groups of 3)
                if len(data) % 3 == 0 and '(' in data_str:
                    # Reshape to vectors
                    data = data.reshape(-1, 3)

                return data

        except Exception as e:
            print(f"  WARNING: Error reading {field_file}: {e}")

        return None

    def compute_gci(self, phi_coarse, phi_medium, phi_fine, r_cm, r_mf, safety_factor=1.25):
        """
        Compute Grid Convergence Index (GCI)

        Parameters:
        -----------
        phi_coarse, phi_medium, phi_fine : float
            Solution values on coarse, medium, and fine meshes
        r_cm : float
            Refinement ratio between coarse and medium meshes
        r_mf : float
            Refinement ratio between medium and fine meshes
        safety_factor : float
            Safety factor (1.25 for 3+ grids, 3.0 for 2 grids)

        Returns:
        --------
        dict : GCI metrics and convergence analysis
        """

        # Richardson extrapolation
        epsilon_cm = phi_medium - phi_coarse
        epsilon_mf = phi_fine - phi_medium

        # Observed order of accuracy
        if epsilon_mf != 0 and epsilon_cm != 0:
            p = np.abs(np.log(np.abs(epsilon_cm / epsilon_mf)) / np.log(r_mf))
        else:
            p = np.nan

        # GCI for medium-fine
        if not np.isnan(p) and phi_fine != 0:
            gci_mf = (safety_factor * np.abs(epsilon_mf / phi_fine)) / (r_mf**p - 1)
        else:
            gci_mf = np.nan

        # GCI for coarse-medium
        if not np.isnan(p) and phi_medium != 0:
            gci_cm = (safety_factor * np.abs(epsilon_cm / phi_medium)) / (r_cm**p - 1)
        else:
            gci_cm = np.nan

        # Richardson extrapolation to zero grid spacing
        if not np.isnan(p):
            phi_extrapolated = phi_fine + epsilon_mf / (r_mf**p - 1)
        else:
            phi_extrapolated = np.nan

        # Check for monotonic convergence
        if epsilon_mf * epsilon_cm > 0:
            convergence_type = "monotonic"
        elif epsilon_mf * epsilon_cm < 0:
            convergence_type = "oscillatory"
        else:
            convergence_type = "uncertain"

        return {
            'order_of_accuracy': p,
            'gci_fine': gci_mf * 100,  # Convert to percentage
            'gci_coarse': gci_cm * 100,
            'extrapolated_value': phi_extrapolated,
            'convergence_type': convergence_type,
            'relative_error_mf': np.abs(epsilon_mf / phi_fine) * 100 if phi_fine != 0 else np.nan,
            'relative_error_cm': np.abs(epsilon_cm / phi_medium) * 100 if phi_medium != 0 else np.nan,
        }

    def analyze_all_cases(self, levels=None):
        """Analyze all mesh levels and compute convergence metrics"""

        if levels is None:
            levels = ['coarse', 'medium', 'fine', 'very_fine']

        print("="*70)
        print("MESH CONVERGENCE ANALYSIS")
        print("="*70)

        # Extract data from all cases
        case_data = {}

        for level in levels:
            case_dir = self.study_dir / level

            if not case_dir.exists():
                print(f"\nWARNING: Case directory not found: {case_dir}")
                continue

            print(f"\nAnalyzing: {level}")
            print("-" * 40)

            # Extract mesh metrics
            mesh_metrics = self.extract_mesh_metrics(case_dir)
            if mesh_metrics:
                print(f"  Mesh cells: {mesh_metrics.get('total_cells', 'N/A'):,}")
                print(f"  Char. cell size (h): {mesh_metrics.get('h', 0)*1e9:.4f} nm")

            # Extract solution metrics
            solution_metrics = self.extract_solution_metrics(case_dir, level)

            case_data[level] = {
                'mesh': mesh_metrics,
                'solution': solution_metrics,
            }

        # Perform convergence analysis
        print("\n" + "="*70)
        print("CONVERGENCE ANALYSIS")
        print("="*70)

        # Analyze key metrics
        metrics_to_analyze = [
            'T_peak', 'Te_peak', 'Tl_peak',
            'U_mag_max', 'p_max',
            'alpha.metal_max',
        ]

        convergence_results = {}

        if len(case_data) >= 3:
            # Need at least 3 meshes for GCI
            coarse = case_data.get('coarse', {}).get('solution', {})
            medium = case_data.get('medium', {}).get('solution', {})
            fine = case_data.get('fine', {}).get('solution', {})

            for metric in metrics_to_analyze:
                if metric in coarse and metric in medium and metric in fine:
                    phi_c = coarse[metric]
                    phi_m = medium[metric]
                    phi_f = fine[metric]

                    print(f"\n{metric}:")
                    print(f"  Coarse:  {phi_c:.6e}")
                    print(f"  Medium:  {phi_m:.6e}")
                    print(f"  Fine:    {phi_f:.6e}")

                    # Compute GCI
                    r_cm = 2.0  # Coarse to medium refinement ratio
                    r_mf = 1.5  # Medium to fine refinement ratio

                    gci_results = self.compute_gci(phi_c, phi_m, phi_f, r_cm, r_mf)
                    convergence_results[metric] = gci_results

                    print(f"  Order of accuracy: {gci_results['order_of_accuracy']:.3f}")
                    print(f"  GCI (fine): {gci_results['gci_fine']:.4f}%")
                    print(f"  Convergence: {gci_results['convergence_type']}")
                    print(f"  Extrapolated value: {gci_results['extrapolated_value']:.6e}")

        # Generate plots
        self._generate_plots(case_data, convergence_results, levels)

        # Generate report
        self._generate_report(case_data, convergence_results, levels)

        print("\n" + "="*70)
        print("Analysis complete!")
        print("  Results saved to: convergenceResults/")
        print("="*70)

    def _generate_plots(self, case_data, convergence_results, levels):
        """Generate convergence plots"""

        output_dir = Path("convergenceResults")
        output_dir.mkdir(exist_ok=True)

        # Plot 1: Mesh refinement vs solution
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle('Mesh Convergence Study - Key Metrics', fontsize=14, fontweight='bold')

        metrics_to_plot = [
            ('T_peak', 'Peak Temperature [K]'),
            ('Te_peak', 'Peak Electron Temperature [K]'),
            ('U_mag_max', 'Max Velocity Magnitude [m/s]'),
            ('p_max', 'Max Pressure [Pa]'),
            ('alpha.metal_max', 'Max Metal Volume Fraction [-]'),
            ('total_cells', 'Total Mesh Cells'),
        ]

        for idx, (metric, ylabel) in enumerate(metrics_to_plot):
            ax = axes[idx // 3, idx % 3]

            x_vals = []
            y_vals = []
            labels = []

            for level in levels:
                if level in case_data:
                    # Get cell size
                    h = case_data[level]['mesh'].get('h', 0) if case_data[level]['mesh'] else 0

                    # Get metric value
                    if metric == 'total_cells':
                        value = case_data[level]['mesh'].get(metric, np.nan) if case_data[level]['mesh'] else np.nan
                    else:
                        value = case_data[level]['solution'].get(metric, np.nan) if case_data[level]['solution'] else np.nan

                    if h > 0 and not np.isnan(value):
                        x_vals.append(h * 1e9)  # Convert to nm
                        y_vals.append(value)
                        labels.append(level)

            if x_vals:
                ax.plot(x_vals, y_vals, 'o-', linewidth=2, markersize=8)
                for i, label in enumerate(labels):
                    ax.annotate(label, (x_vals[i], y_vals[i]),
                               xytext=(5, 5), textcoords='offset points', fontsize=8)

                ax.set_xlabel('Characteristic Cell Size [nm]', fontsize=10)
                ax.set_ylabel(ylabel, fontsize=10)
                ax.grid(True, alpha=0.3)
                ax.set_xscale('log')

        plt.tight_layout()
        plt.savefig(output_dir / 'convergence_plots.png', dpi=300, bbox_inches='tight')
        print(f"\n  ✓ Saved convergence plots: {output_dir / 'convergence_plots.png'}")

        plt.close()

        # Plot 2: GCI bar chart
        if convergence_results:
            fig, ax = plt.subplots(figsize=(10, 6))

            metrics = list(convergence_results.keys())
            gci_values = [convergence_results[m]['gci_fine'] for m in metrics]

            bars = ax.bar(range(len(metrics)), gci_values, color='steelblue', alpha=0.7)

            ax.set_xlabel('Metric', fontsize=12)
            ax.set_ylabel('GCI [%]', fontsize=12)
            ax.set_title('Grid Convergence Index (GCI) for Fine Mesh', fontsize=14, fontweight='bold')
            ax.set_xticks(range(len(metrics)))
            ax.set_xticklabels([m.replace('_', ' ') for m in metrics], rotation=45, ha='right')
            ax.grid(True, axis='y', alpha=0.3)

            # Add value labels on bars
            for bar, value in zip(bars, gci_values):
                height = bar.get_height()
                if not np.isnan(height):
                    ax.text(bar.get_x() + bar.get_width()/2., height,
                           f'{height:.2f}%', ha='center', va='bottom', fontsize=9)

            # Add threshold line (1% is typically good)
            ax.axhline(y=1.0, color='r', linestyle='--', linewidth=2, label='1% threshold')
            ax.legend()

            plt.tight_layout()
            plt.savefig(output_dir / 'gci_values.png', dpi=300, bbox_inches='tight')
            print(f"  ✓ Saved GCI plot: {output_dir / 'gci_values.png'}")

            plt.close()

    def _generate_report(self, case_data, convergence_results, levels):
        """Generate text report of convergence study"""

        output_dir = Path("convergenceResults")
        output_dir.mkdir(exist_ok=True)

        report_file = output_dir / "convergence_report.txt"

        with open(report_file, 'w') as f:
            f.write("="*70 + "\n")
            f.write("MESH CONVERGENCE STUDY REPORT\n")
            f.write("CompInterFoam LIFT Simulation\n")
            f.write("="*70 + "\n\n")

            # Mesh summary
            f.write("MESH SUMMARY\n")
            f.write("-"*70 + "\n")
            f.write(f"{'Level':<15} {'Cells':<15} {'h [nm]':<15} {'Max Non-Orth':<15}\n")
            f.write("-"*70 + "\n")

            for level in levels:
                if level in case_data and case_data[level]['mesh']:
                    mesh = case_data[level]['mesh']
                    cells = mesh.get('total_cells', 'N/A')
                    h = mesh.get('h', 0) * 1e9
                    non_orth = mesh.get('max_non_orthogonality', 'N/A')

                    f.write(f"{level:<15} {cells:<15} {h:<15.4f} {non_orth:<15}\n")

            # Solution summary
            f.write("\n\nSOLUTION SUMMARY\n")
            f.write("-"*70 + "\n")

            for level in levels:
                if level in case_data and case_data[level]['solution']:
                    f.write(f"\n{level.upper()}:\n")
                    sol = case_data[level]['solution']

                    for key, value in sorted(sol.items()):
                        if isinstance(value, (int, float)):
                            f.write(f"  {key:<25} {value:.6e}\n")

            # Convergence analysis
            if convergence_results:
                f.write("\n\nCONVERGENCE ANALYSIS\n")
                f.write("="*70 + "\n")

                for metric, results in convergence_results.items():
                    f.write(f"\n{metric}:\n")
                    f.write("-"*70 + "\n")
                    f.write(f"  Order of accuracy:      {results['order_of_accuracy']:.4f}\n")
                    f.write(f"  GCI (fine mesh):        {results['gci_fine']:.4f}%\n")
                    f.write(f"  GCI (coarse mesh):      {results['gci_coarse']:.4f}%\n")
                    f.write(f"  Convergence type:       {results['convergence_type']}\n")
                    f.write(f"  Extrapolated value:     {results['extrapolated_value']:.6e}\n")
                    f.write(f"  Relative error (m-f):   {results['relative_error_mf']:.4f}%\n")

                f.write("\n" + "="*70 + "\n")
                f.write("CONVERGENCE ASSESSMENT\n")
                f.write("="*70 + "\n")

                # Check if mesh is converged
                gci_threshold = 1.0  # 1% is typically considered good
                converged_metrics = sum(1 for r in convergence_results.values()
                                      if not np.isnan(r['gci_fine']) and r['gci_fine'] < gci_threshold)

                f.write(f"\nMetrics with GCI < {gci_threshold}%: {converged_metrics}/{len(convergence_results)}\n")

                if converged_metrics == len(convergence_results):
                    f.write("\n✓ MESH IS CONVERGED - All metrics meet GCI threshold\n")
                elif converged_metrics > len(convergence_results) / 2:
                    f.write("\n⚠ PARTIAL CONVERGENCE - Most metrics meet threshold\n")
                    f.write("  Consider further refinement for critical metrics\n")
                else:
                    f.write("\n✗ MESH NOT CONVERGED - Further refinement recommended\n")
                    f.write("  Increase mesh resolution and re-run study\n")

            f.write("\n" + "="*70 + "\n")
            f.write("END OF REPORT\n")
            f.write("="*70 + "\n")

        print(f"  ✓ Saved report: {report_file}")


def main():
    """Main execution"""

    analyzer = ConvergenceAnalyzer()

    levels = ['coarse', 'medium', 'fine', 'very_fine']

    if len(sys.argv) > 1:
        levels = sys.argv[1:]

    try:
        analyzer.analyze_all_cases(levels=levels)

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
