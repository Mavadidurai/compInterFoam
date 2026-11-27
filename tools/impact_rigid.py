"""Impact of the jet on the rigid receiver (y = 8 um wall) of the rigid-receiver run, per written time.

    xvfb-run -a pvpython tools/impact_rigid.py <decomposed case> <out.csv>

Cells with alpha.metal > 0.5 count as metal. Columns (lengths in um, full-circle volumes in um^3):
  t_ns            time [ns]
  gap_um          lowest metal cell centre minus 8 um; the first cell layer (centre 8.035 um on L1)
                  gives 0.035, i.e. contact
  r_wall_um       largest radius of a metal cell in the first cell layer (y < 8.071 um): the radius of
                  the metal footprint on the wall, to within one radial cell (49 nm)
  V_wall_um3      metal volume (sum of alpha V, x 72 for the 5 deg wedge) in the first cell layer
  V_below9_um3    metal volume within 1 um of the wall
  v_down_max_ms   largest speed towards the wall of metal within 1 um of it
  v_up_max_ms     largest speed away from the wall of metal within 1 um of it (rebound / splash)
  v_rad_max_ms    largest radial speed of metal within 1 um of the wall (spreading)
  p_wall_max_MPa  largest p in the first cell layer
  Tl_wall_min_K   lowest lattice temperature of metal in the first cell layer
  Tl_wall_mean_K  alpha V weighted lattice temperature of metal in the first cell layer
"""
import os, sys
import numpy as np
from paraview.simple import OpenFOAMReader, MergeBlocks, CellCenters, CellSize, servermanager
from vtk.util.numpy_support import vtk_to_numpy as v2n
case, out = sys.argv[1:3]
r = OpenFOAMReader(FileName=os.path.join(case, 'case.foam'))
r.MeshRegions = ['internalMesh']; r.CaseType = 'Decomposed Case'; r.CellArrays = ['alpha.metal', 'U', 'Tl', 'p']
r.UpdatePipelineInformation()
times = [t for t in r.TimestepValues]
cc = CellCenters(Input=CellSize(Input=MergeBlocks(Input=r)))
nan = float('nan')
rows = []
for t in times:
    cc.UpdatePipeline(t); d = servermanager.Fetch(cc)
    X = v2n(d.GetPoints().GetData()); P = d.GetPointData()
    a = v2n(P.GetArray('alpha.metal')); U = v2n(P.GetArray('U')); V = v2n(P.GetArray('Volume'))
    Tl = v2n(P.GetArray('Tl')); p = v2n(P.GetArray('p'))
    rr = np.hypot(X[:, 0], X[:, 2]); y = X[:, 1]
    ur = (U[:, 0]*X[:, 0] + U[:, 2]*X[:, 2]) / np.maximum(rr, 1e-12)
    m = a > 0.5
    L1 = y < 8.071e-6
    nb = m & (y < 9.0e-6)
    w1 = m & L1
    rows.append((t*1e9,
                 y[m].min()*1e6 - 8.0 if m.any() else nan,
                 rr[w1].max()*1e6 if w1.any() else 0.0,
                 (a*V)[L1].sum()*72*1e18,
                 (a*V)[y < 9.0e-6].sum()*72*1e18,
                 -U[nb, 1].min() if nb.any() else nan,
                 U[nb, 1].max() if nb.any() else nan,
                 ur[nb].max() if nb.any() else nan,
                 p[L1].max()*1e-6,
                 Tl[w1].min() if w1.any() else nan,
                 np.average(Tl[w1], weights=(a*V)[w1]) if w1.any() else nan))
    print(('%.3f ns gap %.3f r_wall %.3f V_wall %.4f V_below9 %.4f v_down %.0f v_up %.0f v_rad %.0f '
           'p_wall %.0f MPa Tl_wall min %.0f mean %.0f') % rows[-1], flush=True)
np.savetxt(out, np.array(rows), delimiter=',', comments='', fmt='%.6g',
           header='t_ns,gap_um,r_wall_um,V_wall_um3,V_below9_um3,v_down_max_ms,v_up_max_ms,v_rad_max_ms,'
                  'p_wall_max_MPa,Tl_wall_min_K,Tl_wall_mean_K')
