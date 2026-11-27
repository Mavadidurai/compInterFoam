"""Per animation frame: gap to the receiver, speed and temperature of the metal near the axis.

    xvfb-run -a pvpython tools/full_run_curves.py <decomposed case> <times.txt> <out.csv>

<times.txt> is the index/time list written by tools/render_full_run.py. Columns:
  t_ns              time [ns]
  gap_um            lowest cell centre with alpha.metal > 0.5 and y > 8.05 um, minus 8 um (receiver top).
                    The first cell layer above the receiver (centre 8.035 um on L1) is excluded, so
                    0.106 um means the metal is in the second layer or lower: contact.
  v_axis_ms         alpha-weighted speed towards the receiver of the metal with r < 0.5 um, 8.05 < y < 20.2 um
  Tl_axis_K         alpha-weighted lattice temperature of the same metal
  Tl_min_flying_K   lowest Tl of cells with alpha > 0.5 for 8.05 < y < 18.5 um (after impact this includes
                    receiver metal thrown up at the rim)
"""
import os, sys
import numpy as np
from paraview.simple import OpenFOAMReader, MergeBlocks, CellCenters, CellSize, servermanager
from vtk.util.numpy_support import vtk_to_numpy as v2n
case, tf, out = sys.argv[1:4]
times = [float(l.split()[1]) for l in open(tf)]
r = OpenFOAMReader(FileName=os.path.join(case, 'case.foam'))
r.MeshRegions = ['internalMesh']; r.CaseType = 'Decomposed Case'; r.CellArrays = ['alpha.metal', 'U', 'Tl']
r.UpdatePipelineInformation()
cc = CellCenters(Input=CellSize(Input=MergeBlocks(Input=r)))
rows = []
for t in times:
    cc.UpdatePipeline(t); d = servermanager.Fetch(cc)
    X = v2n(d.GetPoints().GetData()); P = d.GetPointData()
    a = v2n(P.GetArray('alpha.metal')); U = v2n(P.GetArray('U')); V = v2n(P.GetArray('Volume')); Tl = v2n(P.GetArray('Tl'))
    rr = np.hypot(X[:, 0], X[:, 2]); y = X[:, 1]
    lo = (y > 8.05e-6) & (a > 0.5)
    ax = (rr < 0.5e-6) & (y > 8.05e-6) & (y < 20.2e-6)
    w = a * V * ax
    fl = (y > 8.05e-6) & (y < 18.5e-6) & (a > 0.5)
    rows.append((t * 1e9, y[lo].min() * 1e6 - 8.0 if lo.any() else np.nan,
                 -(w * U[:, 1]).sum() / w.sum() if w.sum() > 0 else np.nan,
                 np.average(Tl[ax], weights=w[ax]) if w.sum() > 0 else np.nan,
                 Tl[fl].min() if fl.any() else np.nan))
    print('%.4f ns gap %.3f um  v_axis %.1f  Tl_axis %.0f  Tl_min_flying %s' % rows[-1], flush=True)
np.savetxt(out, np.array(rows), delimiter=',', comments='', fmt='%.6g',
           header='t_ns,gap_um,v_axis_ms,Tl_axis_K,Tl_min_flying_K')
