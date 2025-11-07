#!/usr/bin/env python3
"""
Mesh Convergence Study - Mesh Generation Script
================================================
Generates multiple blockMeshDict files with different refinement levels
for the LIFT simulation case.

Author: Generated for compInterFoam mesh convergence study
Date: 2025-11-07
"""

import os
import sys
from pathlib import Path

class MeshGenerator:
    """Generates blockMeshDict files with varying refinement levels"""

    def __init__(self, base_case="LiftTest1"):
        self.base_case = base_case

        # Base mesh dimensions (from LiftTest1)
        # Domain: X=[0,50], Y=[0,32.0714], Z=[0,10] micrometers
        self.base_mesh = {
            'substrate':   {'nx': 80, 'ny': 400, 'nz': 40},   # y: 0-8 µm
            'air_gap':     {'nx': 80, 'ny': 400, 'nz': 60},   # y: 8-20 µm
            'ti_film':     {'nx': 80, 'ny': 400, 'nz': 16},   # y: 20-20.0714 µm (71.4 nm, CRITICAL)
            'donor_glass': {'nx': 80, 'ny': 400, 'nz': 60},   # y: 20.0714-32.0714 µm
        }

        # Geometric boundaries (micrometers)
        self.geometry = {
            'x_range': (0, 50e-6),
            'y_substrate': (0, 8e-6),
            'y_air_gap': (8e-6, 20e-6),
            'y_ti_film': (20e-6, 20.0714e-6),  # 71.4 nm thick
            'y_donor': (20.0714e-6, 32.0714e-6),
            'z_range': (0, 10e-6),
        }

        # Refinement levels (relative to base)
        self.refinement_levels = {
            'coarse':     0.5,   # 50% of base resolution
            'medium':     1.0,   # Base resolution (LiftTest1)
            'fine':       1.5,   # 150% of base
            'very_fine':  2.0,   # 200% of base
            'ultra_fine': 2.5,   # 250% of base (if needed)
        }

    def generate_mesh_dict(self, level_name, refinement_factor):
        """Generate blockMeshDict content for a given refinement level"""

        # Calculate cell counts
        mesh = {}
        for region, cells in self.base_mesh.items():
            mesh[region] = {
                'nx': int(cells['nx'] * refinement_factor),
                'ny': int(cells['ny'] * refinement_factor),
                'nz': int(cells['nz'] * refinement_factor),
            }

        # Ensure minimum cells in Ti film (critical region)
        if mesh['ti_film']['nz'] < 8:
            print(f"WARNING: Ti film has only {mesh['ti_film']['nz']} cells in thickness!")
            print(f"         Forcing minimum of 8 cells for {level_name}")
            mesh['ti_film']['nz'] = 8

        total_cells = sum(
            m['nx'] * m['ny'] * m['nz'] for m in mesh.values()
        )

        print(f"\n{level_name.upper()} Mesh (factor={refinement_factor}):")
        print(f"  Substrate:   {mesh['substrate']['nx']} × {mesh['substrate']['ny']} × {mesh['substrate']['nz']} = {mesh['substrate']['nx']*mesh['substrate']['ny']*mesh['substrate']['nz']:,}")
        print(f"  Air gap:     {mesh['air_gap']['nx']} × {mesh['air_gap']['ny']} × {mesh['air_gap']['nz']} = {mesh['air_gap']['nx']*mesh['air_gap']['ny']*mesh['air_gap']['nz']:,}")
        print(f"  Ti film:     {mesh['ti_film']['nx']} × {mesh['ti_film']['ny']} × {mesh['ti_film']['nz']} = {mesh['ti_film']['nx']*mesh['ti_film']['ny']*mesh['ti_film']['nz']:,}")
        print(f"  Donor glass: {mesh['donor_glass']['nx']} × {mesh['donor_glass']['ny']} × {mesh['donor_glass']['nz']} = {mesh['donor_glass']['nx']*mesh['donor_glass']['ny']*mesh['donor_glass']['nz']:,}")
        print(f"  TOTAL CELLS: {total_cells:,}")

        # Ti film cell size
        ti_thickness = 71.4e-9  # meters
        ti_cell_size = ti_thickness / mesh['ti_film']['nz']
        print(f"  Ti film cell height: {ti_cell_size*1e9:.3f} nm")

        return self._create_blockMeshDict_content(mesh, level_name, total_cells)

    def _create_blockMeshDict_content(self, mesh, level_name, total_cells):
        """Create the actual blockMeshDict file content"""

        content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  v2406
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Mesh Convergence Study - {level_name.upper()} mesh
// Total cells: {total_cells:,}
// Generated automatically for LIFT simulation

scale   1e-6;  // Convert from micrometers to meters

// Geometry parameters
x0      0;
x1      50;     // Domain width: 50 µm

y0      0;      // Bottom of substrate
y1      8;      // Top of substrate / bottom of air gap
y2      20;     // Top of air gap / bottom of Ti film
y3      20.0714;// Top of Ti film / bottom of donor glass (71.4 nm thick)
y4      32.0714;// Top of donor glass

z0      0;
z1      10;     // Domain depth: 10 µm

vertices
(
    // Substrate layer (region 0)
    ($x0 $y0 $z0)    // 0
    ($x1 $y0 $z0)    // 1
    ($x1 $y1 $z0)    // 2
    ($x0 $y1 $z0)    // 3
    ($x0 $y0 $z1)    // 4
    ($x1 $y0 $z1)    // 5
    ($x1 $y1 $z1)    // 6
    ($x0 $y1 $z1)    // 7

    // Air gap layer (region 1)
    ($x0 $y1 $z0)    // 8
    ($x1 $y1 $z0)    // 9
    ($x1 $y2 $z0)    // 10
    ($x0 $y2 $z0)    // 11
    ($x0 $y1 $z1)    // 12
    ($x1 $y1 $z1)    // 13
    ($x1 $y2 $z1)    // 14
    ($x0 $y2 $z1)    // 15

    // Ti film layer (region 2) - CRITICAL: 71.4 nm thick
    ($x0 $y2 $z0)    // 16
    ($x1 $y2 $z0)    // 17
    ($x1 $y3 $z0)    // 18
    ($x0 $y3 $z0)    // 19
    ($x0 $y2 $z1)    // 20
    ($x1 $y2 $z1)    // 21
    ($x1 $y3 $z1)    // 22
    ($x0 $y3 $z1)    // 23

    // Donor glass layer (region 3)
    ($x0 $y3 $z0)    // 24
    ($x1 $y3 $z0)    // 25
    ($x1 $y4 $z0)    // 26
    ($x0 $y4 $z0)    // 27
    ($x0 $y3 $z1)    // 28
    ($x1 $y3 $z1)    // 29
    ($x1 $y4 $z1)    // 30
    ($x0 $y4 $z1)    // 31
);

blocks
(
    // Substrate block
    hex (0 1 2 3 4 5 6 7)
    ({mesh['substrate']['nx']} {mesh['substrate']['ny']} {mesh['substrate']['nz']})
    simpleGrading (1 1 1)

    // Air gap block
    hex (8 9 10 11 12 13 14 15)
    ({mesh['air_gap']['nx']} {mesh['air_gap']['ny']} {mesh['air_gap']['nz']})
    simpleGrading (1 1 1)

    // Ti film block (ultra-fine in y-direction)
    hex (16 17 18 19 20 21 22 23)
    ({mesh['ti_film']['nx']} {mesh['ti_film']['ny']} {mesh['ti_film']['nz']})
    simpleGrading (1 1 1)

    // Donor glass block
    hex (24 25 26 27 28 29 30 31)
    ({mesh['donor_glass']['nx']} {mesh['donor_glass']['ny']} {mesh['donor_glass']['nz']})
    simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    substrate
    {{
        type wall;
        faces
        (
            (0 4 5 1)  // Bottom substrate wall
        );
    }}

    donor
    {{
        type wall;
        faces
        (
            (27 26 30 31)  // Top donor wall
        );
    }}

    left
    {{
        type symmetryPlane;
        faces
        (
            (0 3 7 4)    // Substrate left
            (8 11 15 12)  // Air gap left
            (16 19 23 20) // Ti film left
            (24 27 31 28) // Donor left
        );
    }}

    right
    {{
        type symmetryPlane;
        faces
        (
            (1 5 6 2)    // Substrate right
            (9 13 14 10)  // Air gap right
            (17 21 22 18) // Ti film right
            (25 29 30 26) // Donor right
        );
    }}

    front
    {{
        type symmetryPlane;
        faces
        (
            (0 1 2 3)    // Substrate front
            (8 9 10 11)   // Air gap front
            (16 17 18 19) // Ti film front
            (24 25 26 27) // Donor front
        );
    }}

    back
    {{
        type symmetryPlane;
        faces
        (
            (4 7 6 5)    // Substrate back
            (12 15 14 13) // Air gap back
            (20 23 22 21) // Ti film back
            (28 31 30 29) // Donor back
        );
    }}
);

mergePatchPairs
(
);

// ************************************************************************* //
"""
        return content

    def generate_all_meshes(self, output_dir="meshes", levels=None):
        """Generate blockMeshDict files for all refinement levels"""

        if levels is None:
            levels = ['coarse', 'medium', 'fine', 'very_fine']

        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)

        print("="*70)
        print("MESH CONVERGENCE STUDY - Generating Mesh Definitions")
        print("="*70)

        mesh_info = {}

        for level in levels:
            if level not in self.refinement_levels:
                print(f"WARNING: Unknown refinement level '{level}', skipping...")
                continue

            factor = self.refinement_levels[level]
            content = self.generate_mesh_dict(level, factor)

            # Write to file
            output_file = output_path / f"blockMeshDict.{level}"
            with open(output_file, 'w') as f:
                f.write(content)

            print(f"  ✓ Written to: {output_file}")

        print("\n" + "="*70)
        print("Mesh generation complete!")
        print("="*70)

        return mesh_info


def main():
    """Main execution"""

    generator = MeshGenerator()

    # Generate standard mesh levels
    levels = ['coarse', 'medium', 'fine', 'very_fine']

    if len(sys.argv) > 1:
        # Custom levels from command line
        levels = sys.argv[1:]

    generator.generate_all_meshes(levels=levels)

    print("\nNext steps:")
    print("  1. Run setupCases.py to create case directories")
    print("  2. Use runStudy.sh to execute all simulations")
    print("  3. Analyze results with analyzeConvergence.py")


if __name__ == "__main__":
    main()
