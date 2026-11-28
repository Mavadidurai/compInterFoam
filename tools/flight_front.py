"""Position, speed and resolution of the metal in flight below the donor film.

    xvfb-run -a pvpython tools/flight_front.py <case> <t_min [s]> <stride> <out.csv>

Reads every <stride>-th written time after t_min (and always the last one) of a decomposed or
reconstructed axisymmetric 5-degree wedge case, and writes one row per time:

  t_ns                     time [ns]
  y_low_alpha05_um         lowest cell centre with alpha.metal > 0.5 above the receiver (y > 8.05 um) [um]
  y_low_alpha005_um        the same with alpha.metal > 0.05 (includes smeared metal) [um]
  V_below18p5_um3          metal volume below y = 18.5 um, i.e. the dome and the free cap (full circle) [um3]
  v_below18p5_ms           volume-weighted vertical velocity of that metal [m/s, negative = towards the receiver]
  alpha_max_axis_below18p5 largest alpha.metal for r < 1 um and y < 18.5 um: below about 0.95 the cap is
                           thinner than one cell and no longer resolved
  V_metal_total_um3        all metal above the receiver (full circle) [um3]

The receiver top is at y = 8 um and the film at y = 20-20.07 um, as in cases/fsLIFT_332mJ_Ti_meshStudy.
"""
import os
import sys

import numpy as np
from paraview.simple import OpenFOAMReader, MergeBlocks, CellCenters, CellSize, servermanager
from vtk.util.numpy_support import vtk_to_numpy as v2n

case, tmin, stride, out = os.path.abspath(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
f = os.path.join(case, 'case.foam')
open(f, 'a').close()
r = OpenFOAMReader(FileName=f)
r.MeshRegions = ['internalMesh']
r.CaseType = 'Decomposed Case' if os.path.isdir(os.path.join(case, 'processor0')) else 'Reconstructed Case'
r.CellArrays = ['alpha.metal', 'U']
r.UpdatePipelineInformation()
cc = CellCenters(Input=CellSize(Input=MergeBlocks(Input=r)))
W = 360.0 / 5.0
ts = [t for t in r.TimestepValues if t > tmin]
ts = ts[::stride] + ([ts[-1]] if ts and (len(ts) - 1) % stride else [])
rows = []
for t in ts:
    cc.UpdatePipeline(t)
    d = servermanager.Fetch(cc)
    X = v2n(d.GetPoints().GetData())
    P = d.GetPointData()
    a = v2n(P.GetArray('alpha.metal'))
    U = v2n(P.GetArray('U'))
    V = v2n(P.GetArray('Volume'))
    rr = np.hypot(X[:, 0], X[:, 2])
    y = X[:, 1]
    m = y > 8.05e-6
    am = a * V * m
    lo5, lo05 = m & (a > 0.5), m & (a > 0.05)
    det = am * (y < 18.5e-6)
    ax = m & (rr < 1e-6) & (y < 18.5e-6)
    rows.append((t * 1e9,
                 y[lo5].min() * 1e6 if lo5.any() else np.nan,
                 y[lo05].min() * 1e6 if lo05.any() else np.nan,
                 det.sum() * W * 1e18,
                 (det * U[:, 1]).sum() / det.sum() if det.sum() > 0 else np.nan,
                 a[ax].max() if ax.any() else 0.0,
                 am.sum() * W * 1e18))
    print('t = %.3f ns  lowest metal %.3f um  V(y<18.5) %.4f um3  v %.1f m/s  alpha_max %.2f'
          % (rows[-1][:2] + rows[-1][3:6]),
          flush=True)
np.savetxt(out, np.array(rows), delimiter=',', comments='', fmt='%.6g',
           header='t_ns,y_low_alpha05_um,y_low_alpha005_um,V_below18p5_um3,v_below18p5_ms,'
                  'alpha_max_axis_below18p5,V_metal_total_um3')
