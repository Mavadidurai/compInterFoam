#!/usr/bin/env pvpython
"""Motion of the detached droplet and of all metal below the donor film.

    xvfb-run -a pvpython tools/droplet_motion.py cases/<case> t1,t2,...   (times in s) > droplet_motion.csv

For each time: the main detached body (cells with alpha > 0.5, connected, touching the beam axis,
between the receiver at 8.05 um and y = 19 um), its mass, the y position and vertical velocity of its
centre of mass, and its lowest point on the axis; and the mass, vertical momentum and mean vertical
velocity of all metal between 8.05 and 19 um (smeared cells included, weighted by alpha).
Values are for the 5-degree wedge as computed (multiply masses by 72 for the full circle).

The centre-of-mass velocity is the speed of the droplet as a body. The velocity sampled on the axis
(tools/axis_timeseries.py) is the local liquid velocity there, which includes the internal flow of a
droplet that is still changing shape.
"""
import sys, os, numpy as np
from paraview.simple import OpenFOAMReader, ForceTime, MergeBlocks, CleantoGrid, CellSize, Threshold, Connectivity, UpdatePipeline
from paraview import servermanager as sm
from vtk.numpy_interface import dataset_adapter as dsa
from vtk import vtkCellCenters
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from case_props import metal_density
case = sys.argv[1]; RHO = metal_density(case)
f = os.path.join(case, 'case.foam'); open(f, 'a').close()
rd = OpenFOAMReader(FileName=f); rd.CaseType = 'Decomposed Case'; rd.MeshRegions = ['internalMesh']; rd.CellArrays = ['alpha.metal', 'U']
ft = ForceTime(Input=rd); ft.IgnorePipelineTime = 1
cs = CellSize(Input=CleantoGrid(Input=MergeBlocks(Input=ft))); cs.ComputeVolume = 1
th = Threshold(Input=cs); th.Scalars = ['CELLS', 'alpha.metal']
try: th.LowerThreshold = 0.5; th.UpperThreshold = 1.0
except AttributeError: th.ThresholdRange = [0.5, 1.0]
cn = Connectivity(Input=th)
def centres(obj):
    c = vtkCellCenters(); c.SetInputData(obj); c.Update(); return np.asarray(dsa.WrapDataObject(c.GetOutput()).Points)
print("t_ns, body_mass_kg, body_ycom_um, body_Vy_com, body_yaxis_min_um, lower_all_metal_mass_kg, lower_all_metal_Py, lower_all_metal_Vy_mean")
for t in [float(x) for x in sys.argv[2].split(',')]:
    ft.ForcedTime = t
    UpdatePipeline(proxy=cs, time=t); full = sm.Fetch(cs); D = dsa.WrapDataObject(full)
    a = np.asarray(D.CellData['alpha.metal']); U = np.asarray(D.CellData['U']); V = np.asarray(D.CellData['Volume']); X = centres(full)
    low = (X[:, 1] > 8.05e-6) & (X[:, 1] < 19.0e-6)          # everything below the rims (rims sit above ~19.2 um at late times)
    mL = (RHO * a[low] * V[low]).sum(); PL = (RHO * a[low] * V[low] * U[low, 1]).sum()
    UpdatePipeline(proxy=cn, time=t); o = sm.Fetch(cn); W = dsa.WrapDataObject(o)
    rid = np.asarray(W.CellData['RegionId']); aa = np.asarray(W.CellData['alpha.metal']); UU = np.asarray(W.CellData['U']); VV = np.asarray(W.CellData['Volume']); XX = centres(o)
    best = None
    for r in np.unique(rid):
        m = rid == r
        if XX[m, 1].max() > 19.0e-6 or XX[m, 1].min() < 8.05e-6: continue   # skip the film on the glass and the receiver
        if np.abs(XX[m, 0]).min() > 0.1e-6: continue          # main body touches the axis
        mass = (RHO * aa[m] * VV[m]).sum()
        if best is None or mass > best[0]: best = (mass, m)
    if best is None:
        print("%.2f, no detached body on the axis" % (t * 1e9)); continue
    mass, m = best; w = RHO * aa[m] * VV[m]
    yc = (w * XX[m, 1]).sum() / w.sum(); vc = (w * UU[m, 1]).sum() / w.sum()
    ax = np.abs(XX[m, 0]) < 0.06e-6
    ymin = XX[m, 1][ax].min() if ax.any() else float('nan')
    print("%.2f, %.4e, %.4f, %.2f, %.4f, %.4e, %.4e, %.2f" % (t * 1e9, mass, yc * 1e6, vc, ymin * 1e6, mL, PL, PL / mL))
