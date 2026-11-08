#!/usr/bin/env pvpython
"""Momentum of the donor metal along the beam axis at written times.

    xvfb-run -a pvpython tools/metal_momentum.py cases/<case> t1,t2,...   (times in s)

Prints, for each time, P_y = sum(rho_metal * alpha * U_y * V) over all cells above the
receiver (y > 8.05 um), for the 5-degree wedge as computed. Comparing the change of P_y
between two times with the time average of the 'total' force in the solver's
momentumBudget log lines (same wedge) checks Newton's second law for the film.
"""
import os, sys
import numpy as np
from paraview.simple import OpenFOAMReader, MergeBlocks, CellCenters, CellSize, CleantoGrid, servermanager
from vtk.util.numpy_support import vtk_to_numpy as v2n
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from case_props import metal_density
case = os.path.abspath(sys.argv[1]); RHO = metal_density(case)
f = os.path.join(case, 'case.foam'); open(f, 'a').close()
r = OpenFOAMReader(FileName=f); r.MeshRegions = ['internalMesh']
r.CaseType = 'Decomposed Case' if os.path.isdir(os.path.join(case, 'processor0')) else 'Reconstructed Case'
r.CellArrays = ['alpha.metal', 'U']; r.UpdatePipelineInformation()
cs = CellSize(Input=MergeBlocks(Input=r)); cs.ComputeVolume = 1; cc = CellCenters(Input=cs)
for w in [float(x) for x in sys.argv[2].split(',')]:
    t = min(r.TimestepValues, key=lambda a: abs(a - w)); cc.UpdatePipeline(t); d = servermanager.Fetch(cc)
    X = v2n(d.GetPoints().GetData()); Pd = d.GetPointData()
    a = v2n(Pd.GetArray('alpha.metal')); U = v2n(Pd.GetArray('U')); V = v2n(Pd.GetArray('Volume'))
    m = X[:, 1] > 8.05e-6
    print('t = %.6e s   P_y = %.6e kg m/s   metal volume = %.6e m3' % (t, (RHO * a[m] * U[m, 1] * V[m]).sum(), (a[m] * V[m]).sum()))
