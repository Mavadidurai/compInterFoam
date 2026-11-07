#!/usr/bin/env python3
"""Write system/blockMeshDict for one level of the mesh family of cases/fsLIFT_332mJ_Ti_meshStudy.

    python3 tools/make_mesh_level.py <case dir> <level>      level: 1 (coarse), 2 (medium), 3 (fine)

Every cell count is that of level 2 times f = 2**((level - 2) / 2), so cell sizes change by sqrt(2)
in both directions between levels and the gradings (expansion ratios) are the same on every level.
Level 2 has the axial distribution of cases/fsLIFT_332mJ_axisymmetric (film 16 cells of 4.46 nm,
6 nm cells over the 3 um of gap below the film, cells graded to 4.8 nm above the film on the glass
side) and a radial distribution refined near the beam axis: 35 nm cells for r < 3.5 um, where the
film is launched and the droplet forms, then graded by 4x to 8 um (the reference case has 50 nm
everywhere).

Also prints the largest time steps for stages 2 and 3, scaled by f**-1.5 from the reference values
(1 ps, 2.5 ps) because the capillary time-step limit scales with (cell size)**1.5.
"""
import math, os, sys

case, level = sys.argv[1], int(sys.argv[2])
f = 2.0 ** ((level - 2) / 2.0)
n = lambda base: max(1, int(round(base * f)))

R_CORE, R_OUT = 3.5, 4.5                     # um
nr_core, nr_out = n(100), n(50)
nr = nr_core + nr_out
xg = "((%g %d 1) (%g %d 4))" % (R_CORE, nr_core, R_OUT, nr_out)

n_rec = n(24)
g1, g2, g3 = n(170), n(24), n(500)           # gap: 8.5 um at 50 nm, 0.5 um graded to 6 nm, 3 um at 6 nm
n_film = n(16)
s1, s2 = n(26), n(68)                        # glass side: 0.5 um graded from ~4.8 nm, 11.5 um graded

c, s = math.cos(math.radians(2.5)), math.sin(math.radians(2.5))
R = 8.0
X, Z = R * c, R * s
ys = [0.0, 8.0, 20.0, 20.0714, 32.0714]
verts = "\n".join("    (0 %g 0) (%.6f %g %.6f) (%.6f %g %.6f)" % (y, X, y, Z, X, y, -Z) for y in ys)

ny = n_rec + g1 + g2 + g3 + n_film + s1 + s2
text = """/*--------------------------------*- C++ -*----------------------------------*\\
  Mesh level %(level)d of the mesh study (tools/make_mesh_level.py, f = %(f).4f).
  Axisymmetric case, 5 degree wedge. Beam axis = y axis, radial extent 8 um.
  Layers (y in um): receiver 0-8, gap 8-20, titanium film 20-20.0714,
  gas standing in for the donor glass 20.0714-32.0714.
  Radial: %(nr_core)d cells of %(dr_core).1f nm for r < 3.5 um, %(nr_out)d cells graded by 4 to 8 um.
  Film: %(n_film)d cells of %(dy_film).2f nm.  Cells: %(nr)d x %(ny)d = %(ncell)d.
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }
scale 1e-6;
vertices
(
%(verts)s
);
blocks
(
    hex (0 1 2 0 3 4 5 3)     (%(nr)d 1 %(n_rec)d)  simpleGrading (%(xg)s 1 0.05)
    hex (3 4 5 3 6 7 8 6)     (%(nr)d 1 %(ngap)d) simpleGrading (%(xg)s 1 ((8.5 %(g1)d 1) (0.5 %(g2)d 0.12) (3.0 %(g3)d 1)))
    hex (6 7 8 6 9 10 11 9)   (%(nr)d 1 %(n_film)d)  simpleGrading (%(xg)s 1 1)
    hex (9 10 11 9 12 13 14 12) (%(nr)d 1 %(nglass)d) simpleGrading (%(xg)s 1 ((0.5 %(s1)d 10) (11.5 %(s2)d 8)))
);
edges ();
boundary
(
    front     { type wedge;         faces ((0 3 4 1) (3 6 7 4) (6 9 10 7) (9 12 13 10)); }
    back      { type wedge;         faces ((0 2 5 3) (3 5 8 6) (6 8 11 9) (9 11 14 12)); }
    right     { type patch;         faces ((1 4 5 2) (4 7 8 5) (7 10 11 8) (10 13 14 11)); }
    axis      { type empty;         faces ((0 3 3 0) (3 6 6 3) (6 9 9 6) (9 12 12 9)); }
    substrate { type wall;          faces ((0 1 2 0)); }
    donor     { type wall;          faces ((12 14 13 12)); }
);
""" % dict(level=level, f=f, nr_core=nr_core, nr_out=nr_out, dr_core=R_CORE * 1e3 / nr_core,
           n_film=n_film, dy_film=71.4 / n_film, nr=nr, ny=ny, ncell=nr * ny, verts=verts, n_rec=n_rec,
           ngap=g1 + g2 + g3, g1=g1, g2=g2, g3=g3, nglass=s1 + s2, s1=s1, s2=s2, xg=xg)
open(os.path.join(case, "system", "blockMeshDict"), "w").write(text)
print("level %d  f %.4f  cells %d x %d = %d  core dr %.1f nm  film dy %.2f nm  gap dy %.2f nm"
      % (level, f, nr, ny, nr * ny, R_CORE * 1e3 / nr_core, 71.4 / n_film, 3000.0 / g3))
print("MAXDT2=%.3g MAXDT3=%.3g" % (1e-12 * f ** -1.5, 2.5e-12 * f ** -1.5))
