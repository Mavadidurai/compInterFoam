#!/usr/bin/env pvpython
"""
ParaView snapshots of the film dome: full cross-section through the beam
axis (the 5-degree wedge mirrored about r = 0), coloured by U_y and by T, with
the alpha = 0.5 interface drawn in black.

    xvfb-run -a pvpython tools/render_dome.py cases/<case> [t1,t2,...]   (times in s)

y is stretched (STRETCH, default 2) so the 71 nm film stays visible over a
+-4 um field of view; the stretch is printed on every image.
Frame: the bottom of the frame follows the dome. For each time it is set
0.25 um below the lowest axis film point in tools/axis_timeseries_<case>.csv (if that
file exists), and the stretch is reduced so the frame fits the image. Y0_UM
fixes the bottom of the frame by hand.
Output: figures/<case>/dome_<field>_<t>ps.png
"""
import os, sys
from paraview.simple import (Slice, Threshold, OpenFOAMReader, ForceTime, MergeBlocks, Clip, Reflect, Calculator,
    Transform, Show, CreateView, ColorBy, GetColorTransferFunction, GetScalarBar, Render,
    CellDatatoPointData, Contour, Text, SaveScreenshot, UpdatePipeline)

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.abspath(sys.argv[1])
name = os.path.basename(CASE.rstrip('/'))
OUT = os.path.join(HERE, '..', 'figures', name); os.makedirs(OUT, exist_ok=True)
STRETCH = float(os.environ.get('STRETCH', '2'))
Y0_DEFAULT, Y1, R = 19.35e-6, 20.15e-6, 4.0e-6
IMG_H_UM = 8.0 * 760 / 1600 - 0.25        # usable frame height in image units (8 um wide image)
def frame(t):
    """bottom of the frame [m] and y stretch for time t."""
    if os.environ.get('Y0_UM'):
        y0 = float(os.environ['Y0_UM']) * 1e-6
    else:
        y0 = Y0_DEFAULT
        csv = os.path.join(HERE, 'axis_timeseries_%s.csv' % os.path.basename(os.path.abspath(sys.argv[1]).rstrip('/')))
        if os.path.exists(csv):
            import numpy as np
            d = np.loadtxt(csv, delimiter=',', ndmin=2)
            d = d[d[:, 0] <= t * (1 + 1e-9)]
            if len(d):
                y0 = min(Y0_DEFAULT, d[-1, 1] - 0.25e-6)
    st = min(STRETCH, IMG_H_UM / ((Y1 - y0) * 1e6))
    return y0, st
foam = os.path.join(CASE, 'case.foam'); open(foam, 'a').close()
rd = OpenFOAMReader(FileName=foam)
rd.CaseType = 'Decomposed Case' if os.path.isdir(os.path.join(CASE, 'processor0')) else 'Reconstructed Case'
rd.MeshRegions = ['internalMesh']; rd.CellArrays = ['alpha.metal', 'U', 'T']
rd.UpdatePipelineInformation()
times = [t for t in rd.TimestepValues if t > 0]
want = [float(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else [times[-1]]
pick = sorted({min(times, key=lambda a: abs(a - w)) for w in want})

def render(t, field, rng, label, preset):
    v = CreateView('RenderView'); v.ViewSize = [1600, 760]; v.Background = [1, 1, 1]
    try: v.UseColorPaletteForBackground = 0
    except AttributeError: pass
    v.OrientationAxesVisibility = 0
    ft = ForceTime(Input=rd); ft.ForcedTime = t; ft.IgnorePipelineTime = 1
    src = MergeBlocks(Input=ft)
    Y0, ST = frame(t)
    for org, nrm, inv in (([0, Y0, 0], [0, 1, 0], 0), ([0, Y1, 0], [0, 1, 0], 1), ([R, 0, 0], [1, 0, 0], 1)):
        c = Clip(Input=src); c.ClipType = 'Plane'; c.ClipType.Origin = org; c.ClipType.Normal = nrm; c.Invert = inv; src = c
    # 2-D section: slice the 5-degree wedge on its mid-plane (z = 0), then mirror about r = 0
    sl = Slice(Input=src); sl.SliceType = 'Plane'; sl.SliceType.Origin = [0, 0, 0]; sl.SliceType.Normal = [0, 0, 1]
    calc = Calculator(Input=sl); calc.AttributeType = 'Cell Data'; calc.ResultArrayName = 'Uy'; calc.Function = 'U_Y'
    ref = Reflect(Input=calc); ref.Plane = 'X'; ref.Center = 0.0; ref.CopyInput = 1
    tr = Transform(Input=ref); tr.Transform = 'Transform'; tr.Transform.Scale = [1e6, 1e6 * ST, 1e6]
    UpdatePipeline(proxy=tr)
    gas = Show(tr, v); ColorBy(gas, None); gas.DiffuseColor = [0.93, 0.93, 0.91]; gas.AmbientColor = [0.93, 0.93, 0.91]
    gas.Ambient = 1.0; gas.Diffuse = 0.0
    met = Threshold(Input=tr); met.Scalars = ['CELLS', 'alpha.metal']
    try:
        met.LowerThreshold = 0.5; met.UpperThreshold = 1.0
    except AttributeError:
        met.ThresholdRange = [0.5, 1.0]
    d = Show(met, v); ColorBy(d, ('CELLS', field))
    lut = GetColorTransferFunction(field); lut.ApplyPreset(preset, True)
    lut.RescaleTransferFunction(*rng); lut.AutomaticRescaleRangeMode = 'Never'
    d.Ambient = 1.0; d.Diffuse = 0.0; d.SetScalarBarVisibility(v, True)
    ct = Contour(Input=CellDatatoPointData(Input=tr)); ct.ContourBy = ['POINTS', 'alpha.metal']; ct.Isosurfaces = [0.5]
    cd = Show(ct, v); ColorBy(cd, None); cd.DiffuseColor = [0, 0, 0]; cd.AmbientColor = [0, 0, 0]; cd.LineWidth = 2.0
    Render(v); v.ResetCamera()
    yc = 0.5 * (Y0 + Y1) * 1e6 * ST
    v.CameraPosition = [0, yc, 100]; v.CameraFocalPoint = [0, yc, 0]; v.CameraViewUp = [0, 1, 0]
    v.CameraParallelProjection = 1; v.CameraParallelScale = max(0.5 * (Y1 - Y0) * 1e6 * ST * 1.18, 1.06 * R * 1e6 * 760 / 1600)
    bar = GetScalarBar(lut, v); bar.Title = label; bar.ComponentTitle = ''
    bar.TitleColor = [0, 0, 0]; bar.LabelColor = [0, 0, 0]; bar.TitleFontSize = 16; bar.LabelFontSize = 14
    for txt, pos, fs in (
        ('%s   t = %.0f ps   %s' % (name, t * 1e12, label), [0.02, 0.93], 22),
        ('section through the beam axis (r = 0 at centre), width 8 um; y stretched x%g   |   coloured: metal (metal fraction > 0.5), grey: gas, black line: metal fraction = 0.5'
         '   |   glass side up, receiver below' % round(ST, 2), [0.02, 0.02], 14),
    ):
        tx = Text(Text=txt); td = Show(tx, v); td.WindowLocation = 'Any Location'; td.Position = pos; td.FontSize = fs; td.Color = [0, 0, 0]
    Render(v)
    out = os.path.join(OUT, 'dome_%s_%05.0fps.png' % (field.replace('.', '_'), t * 1e12))
    SaveScreenshot(out, v, ImageResolution=[1600, 760]); print('saved', out, flush=True)

for t in pick:
    render(t, 'Uy', (-200, 200), 'vertical velocity (m/s, negative = towards the receiver)', 'Cool to Warm')
    render(t, 'T', (300, 6000), 'T [K]', 'Inferno (matplotlib)')
