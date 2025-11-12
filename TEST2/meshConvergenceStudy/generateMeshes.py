#!/usr/bin/env python3
"""
Generate blockMeshDict files for mesh convergence study - TEST2 (Graded Mesh)

This script creates multiple refinement levels based on TEST2's graded mesh configuration.
TEST2 uses non-uniform grading to refine critical regions:
- Substrate: grading (1 1 0.5) - finer at top
- Air gap:   grading (1 1 2.0) - finer at bottom
- Ti film:   grading (1 1 0.67) - finer at top

Refinement levels:
- coarse:    0.5x (40×200×20 substrate, 40×200×40 air, 40×200×18 Ti)
- medium:    1.0x (80×400×40 substrate, 80×400×80 air, 80×400×36 Ti) [baseline]
- fine:      1.5x (120×600×60 substrate, 120×600×120 air, 120×600×54 Ti)
- very_fine: 2.0x (160×800×80 substrate, 160×800×160 air, 160×800×72 Ti)
"""

import os
import sys

# Base mesh configuration from TEST2 (graded mesh)
BASE_CONFIG = {
    'substrate': {'Nx': 80, 'Ny': 400, 'Nz': 40, 'grading': '(1 1 0.5)'},
    'air_gap': {'Nx': 80, 'Ny': 400, 'Nz': 80, 'grading': '(1 1 2.0)'},
    'ti_film': {'Nx': 80, 'Ny': 400, 'Nz': 36, 'grading': '(1 1 0.67)'}
}

# Refinement factors
REFINEMENT_LEVELS = {
    'coarse': 0.5,
    'medium': 1.0,
    'fine': 1.5,
    'very_fine': 2.0
}

BLOCKMESH_TEMPLATE = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2406                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Generated for mesh convergence study - TEST2 GRADED MESH
// Refinement level: {level}
// Refinement factor: {factor}x baseline
// Total cells: ~{total_cells:.2f}M

scale 1e-6; // coordinates in micrometers

// OPTIMIZED GRADED MESH FOR fs-LIFT EJECTION
// ==========================================
// Improvements over uniform mesh:
//   1. Ti film: Finer at top (peak laser absorption, laser enters from top)
//   2. Air gap: Finer near Ti (jet formation)
//   3. Substrate: Finer at top (interface with air)

vertices
(
    // y0 = 0.0000 (substrate/receiver bottom)
    (0   0.0000  0)  (50  0.0000  0)  (50  0.0000 10)  (0   0.0000 10)

    // y1 = 8.0000 (substrate top / air gap bottom)
    (0   8.0000  0)  (50  8.0000  0)  (50  8.0000 10)  (0   8.0000 10)

    // y2 = 20.0000 (air gap top / Ti film bottom)
    (0  20.0000  0)  (50 20.0000  0)  (50 20.0000 10)  (0  20.0000 10)

    // y3 = 20.0714 (Ti film top) - 71.4 nm thick film
    (0  20.0714  0)  (50 20.0714  0)  (50 20.0714 10)  (0  20.0714 10)
);

blocks
(
    // SUBSTRATE (receiver, 8 μm thick)
    // Grading: (1 1 0.5) = finer cells at TOP (near air gap interface)
    // Benefits: Better thermal coupling, smooth transition to air gap
    hex ( 0  3  2  1   4  7  6  5)
        ({sub_nx} {sub_ny} {sub_nz})
        simpleGrading {sub_grading}

    // AIR GAP (12 μm, jet trajectory zone) - CRITICAL FOR EJECTION
    // Grading: (1 1 2.0) = finer cells at BOTTOM (near Ti film)
    // Benefits:
    //   - Captures jet formation with fine resolution near Ti
    //   - Resolves recoil pressure gradients
    //   - Bottom cells finer than top cells
    hex ( 4  7  6  5   8 11 10  9)
        ({air_nx} {air_ny} {air_nz})
        simpleGrading {air_grading}

    // Ti FILM (71.4 nm, ablation zone) - MOST CRITICAL
    // Grading: (1 1 0.67) = finer cells at TOP (laser entry surface)
    //
    // IMPORTANT: Laser enters from TOP (donor side, y3) going DOWN
    // Beer-Lambert absorption: I(y) = I₀ exp(-α(y3-y))
    // Peak absorption at top surface → needs finest cells there
    //
    // Result: Finer cells at top (y3) → coarser at bottom (y2)
    // Gives better resolution in first penetration depth
    hex ( 8 11 10  9  12 15 14 13)
        ({ti_nx} {ti_ny} {ti_nz})
        simpleGrading {ti_grading}
);

edges ();

boundary
(
    // Symmetry planes
    left          {{ type symmetryPlane; faces ((0 3 7 4)(4 7 11 8)(8 11 15 12)); }}
    right         {{ type symmetryPlane; faces ((1 5 6 2)(5 9 10 6)(9 13 14 10)); }}
    front         {{ type symmetryPlane; faces ((0 4 5 1)(4 8 9 5)(8 12 13 9)); }}
    back          {{ type symmetryPlane; faces ((3 2 6 7)(7 6 10 11)(11 10 14 15)); }}

    // Physical walls
    substrate     {{ type wall;  faces ((0 3 2 1)); }}      // BOTTOM
    donor         {{ type wall;  faces ((12 15 14 13)); }}   // TOP
);

mergePatchPairs();

// Cell count information:
// Substrate: {sub_nx} × {sub_ny} × {sub_nz} (graded 0.5 - finer at top)
// Air gap:   {air_nx} × {air_ny} × {air_nz} (graded 2.0 - finer at bottom)
// Ti film:   {ti_nx} × {ti_ny} × {ti_nz} (graded 0.67 - finer at top)

// Nominal cell sizes (actual varies due to grading):
// Substrate: {sub_dx:.3f} × {sub_dy:.3f} × {sub_dz:.3f} μm (avg)
// Air gap:   {air_dx:.3f} × {air_dy:.3f} × {air_dz:.3f} μm (avg)
// Ti film:   {ti_dx:.3f} × {ti_dy:.3f} × {ti_dz:.4f} μm ({ti_dz_nm:.2f} nm avg in z)

// ************************************************************************* //
"""

def calculate_cells(base_count, factor):
    """Calculate cell count with refinement factor, ensuring even numbers"""
    refined = int(base_count * factor)
    # Ensure even number for better mesh quality
    if refined % 2 == 1:
        refined += 1
    return refined

def generate_mesh(level, factor, output_dir):
    """Generate blockMeshDict for a given refinement level"""

    # Calculate cell counts for each block
    sub_nx = calculate_cells(BASE_CONFIG['substrate']['Nx'], factor)
    sub_ny = calculate_cells(BASE_CONFIG['substrate']['Ny'], factor)
    sub_nz = calculate_cells(BASE_CONFIG['substrate']['Nz'], factor)

    air_nx = calculate_cells(BASE_CONFIG['air_gap']['Nx'], factor)
    air_ny = calculate_cells(BASE_CONFIG['air_gap']['Ny'], factor)
    air_nz = calculate_cells(BASE_CONFIG['air_gap']['Nz'], factor)

    ti_nx = calculate_cells(BASE_CONFIG['ti_film']['Nx'], factor)
    ti_ny = calculate_cells(BASE_CONFIG['ti_film']['Ny'], factor)
    ti_nz = calculate_cells(BASE_CONFIG['ti_film']['Nz'], factor)

    # Calculate total cells (in millions)
    total_cells = (sub_nx * sub_ny * sub_nz +
                   air_nx * air_ny * air_nz +
                   ti_nx * ti_ny * ti_nz) / 1e6

    # Calculate nominal cell sizes (domain dimensions / cell count)
    # Note: Actual cell sizes vary due to grading
    # X: 0 to 50 μm, Y: substrate 8 μm, air 12 μm, Ti 0.0714 μm, Z: 0 to 10 μm
    sub_dx = 50.0 / sub_nx
    sub_dy = 8.0 / sub_ny
    sub_dz = 10.0 / sub_nz

    air_dx = 50.0 / air_nx
    air_dy = 12.0 / air_ny
    air_dz = 10.0 / air_nz

    ti_dx = 50.0 / ti_nx
    ti_dy = 0.0714 / ti_ny  # 71.4 nm
    ti_dz = 10.0 / ti_nz
    ti_dz_nm = ti_dy * 1000  # Convert to nm

    # Format the blockMeshDict
    content = BLOCKMESH_TEMPLATE.format(
        level=level,
        factor=factor,
        total_cells=total_cells,
        sub_nx=sub_nx, sub_ny=sub_ny, sub_nz=sub_nz,
        sub_grading=BASE_CONFIG['substrate']['grading'],
        air_nx=air_nx, air_ny=air_ny, air_nz=air_nz,
        air_grading=BASE_CONFIG['air_gap']['grading'],
        ti_nx=ti_nx, ti_ny=ti_ny, ti_nz=ti_nz,
        ti_grading=BASE_CONFIG['ti_film']['grading'],
        sub_dx=sub_dx, sub_dy=sub_dy, sub_dz=sub_dz,
        air_dx=air_dx, air_dy=air_dy, air_dz=air_dz,
        ti_dx=ti_dx, ti_dy=ti_dy, ti_dz=ti_dz, ti_dz_nm=ti_dz_nm
    )

    # Write to file
    filename = os.path.join(output_dir, f'blockMeshDict.{level}')
    with open(filename, 'w') as f:
        f.write(content)

    print(f"Generated {level:10s} ({factor:.1f}x): {total_cells:6.2f}M cells - {filename}")
    print(f"  Substrate: {sub_nx}×{sub_ny}×{sub_nz} (graded 0.5)")
    print(f"  Air gap:   {air_nx}×{air_ny}×{air_nz} (graded 2.0)")
    print(f"  Ti film:   {ti_nx}×{ti_ny}×{ti_nz} (graded 0.67, avg z-cell: {ti_dz_nm:.2f} nm)")
    print()

def main():
    """Main function to generate all mesh levels"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    meshes_dir = os.path.join(script_dir, 'meshes')

    # Create output directory
    os.makedirs(meshes_dir, exist_ok=True)

    print("="*70)
    print("Mesh Convergence Study - TEST2 (Graded Mesh)")
    print("="*70)
    print()

    # Generate meshes for all levels
    for level, factor in REFINEMENT_LEVELS.items():
        generate_mesh(level, factor, meshes_dir)

    print("="*70)
    print("All mesh files generated successfully!")
    print(f"Output directory: {meshes_dir}")
    print("="*70)

if __name__ == "__main__":
    main()
