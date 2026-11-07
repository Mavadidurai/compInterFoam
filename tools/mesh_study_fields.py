#!/usr/bin/env pvpython
"""Mesh-independent field quantities of a run at every written time, for the mesh study.

    xvfb-run -a pvpython tools/mesh_study_fields.py <case dir> [rho_metal, default 4506]

Writes <case dir>/mesh_study_fields.csv with, at each written time t:
  V_metal      metal volume above the receiver (y > 8.05 um), full circle [um3]
  Py           axial momentum of that metal, full circle [kg m/s]
  KE           kinetic energy of that metal, full circle [J]
  v_r1, v_r25  metal-volume-weighted mean axial velocity for r < 1 um and r < 2.5 um [m/s]
  y_bot_axis   lowest y at which alpha.metal crosses 0.5 in the cell column on the axis [um]
  h_axis       metal thickness on the axis, sum of alpha*dy in that column [nm]
  Tl_max       largest lattice temperature in cells with alpha.metal > 0.5 [K]
The wedge is 5 degrees, so full-circle values are the wedge sums times 72.
All definitions use physical positions only, so they can be compared between meshes.
"""
import os, sys
import numpy as np
from paraview.simple import OpenFOAMReader, MergeBlocks, CellCenters, CellSize, servermanager
from vtk.util.numpy_support import vtk_to_numpy as v2n

case = os.path.abspath(sys.argv[1]); RHO = float(sys.argv[2]) if len(sys.argv) > 2 else 4506.0
f = os.path.join(case, 'case.foam'); open(f, 'a').close()
r = OpenFOAMReader(FileName=f); r.MeshRegions = ['internalMesh']
r.CaseType = 'Decomposed Case' if os.path.isdir(os.path.join(case, 'processor0')) else 'Reconstructed Case'
r.CellArrays = ['alpha.metal', 'U', 'Tl']; r.UpdatePipelineInformation()
cs = CellSize(Input=MergeBlocks(Input=r)); cs.ComputeVolume = 1; cc = CellCenters(Input=cs)
W = 360.0 / 5.0
rows = []
for t in r.TimestepValues:
    if t <= 0:
        continue
    cc.UpdatePipeline(t); d = servermanager.Fetch(cc)
    X = v2n(d.GetPoints().GetData()); P = d.GetPointData()
    a = v2n(P.GetArray('alpha.metal')); U = v2n(P.GetArray('U')); V = v2n(P.GetArray('Volume'))
    Tl = v2n(P.GetArray('Tl'))
    rr = np.hypot(X[:, 0], X[:, 2]); y = X[:, 1]
    m = y > 8.05e-6
    am = a[m] * V[m]
    Vm = am.sum() * W
    Py = (RHO * am * U[m, 1]).sum() * W
    KE = 0.5 * RHO * (am * (U[m] ** 2).sum(axis=1)).sum() * W
    vavg = []
    for R in (1e-6, 2.5e-6):
        k = m & (rr < R)
        w = a[k] * V[k]
        vavg.append((w * U[k, 1]).sum() / w.sum() if w.sum() > 0 else np.nan)
    # cell column on the axis: the cells with the smallest radius
    col = m & (rr < rr[m].min() * 1.5 + 1e-12)
    o = np.argsort(y[col]); yc = y[col][o]; ac = a[col][o]
    dy = np.gradient(yc) if yc.size > 1 else np.array([0.0])
    h_axis = (ac * dy).sum()
    above = np.flatnonzero(ac >= 0.5)
    if above.size and above[0] > 0:
        i = above[0]
        y_bot = yc[i - 1] + (0.5 - ac[i - 1]) / (ac[i] - ac[i - 1]) * (yc[i] - yc[i - 1])
    else:
        y_bot = np.nan
    hot = a > 0.5
    rows.append((t, Vm * 1e18, Py, KE, vavg[0], vavg[1], y_bot * 1e6, h_axis * 1e9,
                 Tl[hot].max() if hot.any() else np.nan))
    print('t = %.4e s  V = %.4f um3  Py = %.4e  v(r<1um) = %.2f m/s  y_bot = %.4f um  h = %.1f nm'
          % (t, Vm * 1e18, Py, vavg[0], y_bot * 1e6, h_axis * 1e9), flush=True)
np.savetxt(os.path.join(case, 'mesh_study_fields.csv'), np.array(rows), delimiter=',',
           header='t_s,V_metal_um3,Py_kgms,KE_J,v_r1_ms,v_r25_ms,y_bot_axis_um,h_axis_nm,Tl_max_K',
           comments='')
