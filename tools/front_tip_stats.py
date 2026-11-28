"""Leading metal of the flying jet at chosen times (receiver top at y = 8 um).

    xvfb-run -a pvpython tools/front_tip_stats.py <decomposed case> <out.csv> <t1_ns> <t2_ns> ...

Flying metal: cell centres with 8.05 < y < 19 um. Columns (full-circle volumes, sum alpha V x 72):
  gap_a05_um / gap_a01_um   lowest flying cell with alpha > 0.5 / > 0.1, minus 8 um
  V_8to9/10/12_um3          flying metal volume below y = 9 / 10 / 12 um
  V_tip_um3                 metal (alpha > 0.5) within 0.5 um above the lowest alpha > 0.5 cell
  v_tip_max_ms / _mean_ms   largest / alpha-V weighted speed towards the receiver of that tip metal
  r_tip_um                  largest radius of a tip cell centre
"""
import os, sys
import numpy as np
from paraview.simple import OpenFOAMReader, MergeBlocks, CellCenters, CellSize, servermanager
from vtk.util.numpy_support import vtk_to_numpy as v2n
case, out = sys.argv[1], sys.argv[2]; want = [float(x)*1e-9 for x in sys.argv[3:]]
r = OpenFOAMReader(FileName=os.path.join(case, 'case.foam'))
r.MeshRegions = ['internalMesh']; r.CaseType = 'Decomposed Case'; r.CellArrays = ['alpha.metal', 'U']
r.UpdatePipelineInformation(); tv = np.array(r.TimestepValues)
cc = CellCenters(Input=CellSize(Input=MergeBlocks(Input=r)))
rows = []
for w in want:
    t = tv[np.argmin(abs(tv - w))]
    cc.UpdatePipeline(t); d = servermanager.Fetch(cc)
    X = v2n(d.GetPoints().GetData()); P = d.GetPointData()
    a = v2n(P.GetArray('alpha.metal')); U = v2n(P.GetArray('U')); V = v2n(P.GetArray('Volume'))
    y = X[:, 1]*1e6; rr = np.hypot(X[:, 0], X[:, 2])*1e6
    fly = (y > 8.05) & (y < 19.0)
    m5 = fly & (a > 0.5); m1 = fly & (a > 0.1)
    def vol(mask): return (a*V)[mask].sum()*72*1e18
    yl5 = y[m5].min() if m5.any() else np.nan; yl1 = y[m1].min() if m1.any() else np.nan
    tip = m5 & (y < yl5 + 0.5)            # metal within 0.5 um of the lowest alpha>0.5 cell
    rows.append((t*1e9, yl5-8, yl1-8,
                 vol(fly & (y < 9)), vol(fly & (y < 10)), vol(fly & (y < 12)),
                 vol(tip), -U[tip, 1].min() if tip.any() else np.nan,
                 np.average(-U[tip, 1], weights=(a*V)[tip]) if tip.any() else np.nan,
                 rr[tip].max() if tip.any() else np.nan))
    print('%.3f gap5 %.3f gap1 %.3f V<9 %.4f V<10 %.4f V<12 %.4f Vtip %.5f vtip_max %.0f vtip_mean %.0f rtip %.3f' % rows[-1], flush=True)
np.savetxt(out, np.array(rows), delimiter=',', comments='', fmt='%.6g',
           header='t_ns,gap_a05_um,gap_a01_um,V_8to9_um3,V_8to10_um3,V_8to12_um3,V_tip_um3,v_tip_max_ms,v_tip_mean_ms,r_tip_um')
