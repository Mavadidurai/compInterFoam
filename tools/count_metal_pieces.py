#!/usr/bin/env pvpython
"""Count separate metal pieces (alpha > 0.5) in a written time of a case.

    xvfb-run -a pvpython tools/count_metal_pieces.py cases/<case> <time>

CleantoGrid merges the points duplicated at processor boundaries. Without it a
decomposed case shows false splits along those boundaries (an intact dome can
appear as two pieces).
Region with y < 8.05 um = receiver; every other region = donor film or a detached piece.
"""
import sys, os
from paraview.simple import OpenFOAMReader, ForceTime, MergeBlocks, Threshold, Connectivity, UpdatePipeline, CleantoGrid
from paraview import servermanager as sm
from vtk.numpy_interface import dataset_adapter as dsa
import numpy as np
case, t = sys.argv[1], float(sys.argv[2])
foam = os.path.join(case, 'case.foam'); open(foam, 'a').close()
rd = OpenFOAMReader(FileName=foam); rd.CaseType = 'Decomposed Case'; rd.MeshRegions = ['internalMesh']; rd.CellArrays = ['alpha.metal', 'U']
ft = ForceTime(Input=rd); ft.ForcedTime = t; ft.IgnorePipelineTime = 1
th = Threshold(Input=CleantoGrid(Input=MergeBlocks(Input=ft))); th.Scalars = ['CELLS', 'alpha.metal']
try:
    th.LowerThreshold = 0.5; th.UpperThreshold = 1.0
except AttributeError:
    th.ThresholdRange = [0.5, 1.0]
cn = Connectivity(Input=th); UpdatePipeline(proxy=cn, time=t)
d = dsa.WrapDataObject(sm.Fetch(cn))
rid = np.asarray(d.CellData['RegionId']); U = np.asarray(d.CellData['U'])
from vtk import vtkCellCenters
cc = vtkCellCenters(); cc.SetInputData(sm.Fetch(cn)); cc.Update()
xyz = dsa.WrapDataObject(cc.GetOutput()).Points
y = np.asarray(xyz)[:, 1]
for r in np.unique(rid):
    m = rid == r
    print('region %d: cells %d, y %.3f-%.3f um, mean Uy %.1f m/s' % (r, m.sum(), y[m].min()*1e6, y[m].max()*1e6, U[m, 1].mean()))
