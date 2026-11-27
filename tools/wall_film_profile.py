"""Metal film on the rigid receiver, per animation frame: equivalent thickness against radius.

    xvfb-run -a pvpython tools/wall_film_profile.py <decomposed case> <times.txt> <out.csv>

On the coarse mesh the metal that lands spreads into a film thinner than one cell (71 nm), so the
alpha.metal = 0.5 surface drawn in the animation no longer shows it. For each radial cell column, the
equivalent film thickness is h(r) = sum of alpha.metal * dy over the first two cell layers above the
wall (y < 8.142 um). Metal colder than 400 K is left out: it is receiver metal that lay above y = 8 um
in run A when the receiver was cut off at 21.5 ns (about 0.009 um^3, full circle), not deposited metal.
Output: one row per frame, t_ns followed by h in nm at the radii of the header (cell centres, um).
"""
import os, sys
import numpy as np
from paraview.simple import OpenFOAMReader, MergeBlocks, CellCenters, CellSize, servermanager
from vtk.util.numpy_support import vtk_to_numpy as v2n
case, tf, out = sys.argv[1:4]
times = [float(l.split()[1]) for l in open(tf)]
r = OpenFOAMReader(FileName=os.path.join(case, 'case.foam'))
r.MeshRegions = ['internalMesh']; r.CaseType = 'Decomposed Case'; r.CellArrays = ['alpha.metal', 'Tl']
r.UpdatePipelineInformation()
cc = CellCenters(Input=CellSize(Input=MergeBlocks(Input=r)))
rows = []; R = None
for t in times:
    cc.UpdatePipeline(t); d = servermanager.Fetch(cc)
    X = v2n(d.GetPoints().GetData()); P = d.GetPointData()
    a = v2n(P.GetArray('alpha.metal')); V = v2n(P.GetArray('Volume')); Tl = v2n(P.GetArray('Tl'))
    rr = np.hypot(X[:, 0], X[:, 2]); y = X[:, 1]
    m = (y < 8.142e-6) & (rr < 3e-6)
    if R is None:
        R = np.unique(np.round(rr[m & (y < 8.071e-6)] * 1e9) / 1e3)  # radial cell centres, um
    ri = np.searchsorted(R, np.round(rr[m] * 1e9) / 1e3).clip(0, len(R) - 1)
    hy = 71e-9  # cell height near the wall on L1
    h = np.bincount(ri, weights=(a[m] * (Tl[m] > 400.0)) * hy, minlength=len(R)) * 1e9
    rows.append([t * 1e9] + list(h))
    print('%.3f ns  max h %.1f nm  r(h>2 nm) %.2f um' % (t * 1e9, h.max(), R[h > 2].max() if (h > 2).any() else 0), flush=True)
np.savetxt(out, np.array(rows), delimiter=',', comments='', fmt='%.5g', header='t_ns,' + ','.join('%.4f' % v for v in R))
