#!/usr/bin/env pvpython
"""
ParaView gallery of a pocket run: sections through the beam axis of
alpha (full 0-1 field), |U| with direction arrows, radial velocity U_r and
pressure, at several times; laser-heating close-ups (Te, Tl, Q_laser) at
early times; and 3-D views of the alpha = 0.5 surface revolved about the axis.

    xvfb-run -a pvpython tools/render_gallery.py cases/<case>

Output: figures/<case>/gallery/*.png
Section frame: r = -4..4 um (mirrored wedge, mid-plane z = 0); the bottom of
the frame follows the lowest metal (alpha > 0.5) above the receiver.
"""
import os, sys
import numpy as np
from paraview.simple import (OpenFOAMReader, ForceTime, MergeBlocks, Clip, Slice, Calculator, Reflect,
    Transform, Show, Hide, CreateView, ColorBy, GetColorTransferFunction, GetScalarBar, Render,
    CellDatatoPointData, Contour, Text, SaveScreenshot, UpdatePipeline, Glyph, Threshold,
    RotationalExtrusion, Delete)
from paraview import servermanager as sm
from vtk.numpy_interface import dataset_adapter as dsa

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.abspath(sys.argv[1])
NAME = os.path.basename(CASE.rstrip('/'))
OUT = os.path.join(HERE, '..', 'figures', NAME, 'gallery'); os.makedirs(OUT, exist_ok=True)
W, H = 1600, 760
R = 4.0e-6
Y1 = 20.15e-6
TIMES = [float(x) for x in os.environ.get('TIMES', '2e-10,1e-9,3e-9,6e-9,1e-8,1.34e-8,1.5e-8,1.98e-8').split(',')]
EARLY = [float(x) for x in os.environ.get('EARLY', '5e-13,1e-12,2e-12').split(',')]

foam = os.path.join(CASE, 'case.foam'); open(foam, 'a').close()
rd = OpenFOAMReader(FileName=foam)
rd.CaseType = 'Decomposed Case' if os.path.isdir(os.path.join(CASE, 'processor0')) else 'Reconstructed Case'
rd.MeshRegions = ['internalMesh']
rd.CellArrays = ['alpha.metal', 'U', 'p', 'T', 'Te', 'Tl', 'Q_laser']
rd.UpdatePipelineInformation()
avail = [t for t in rd.TimestepValues if t > 0]
nearest = lambda w: min(avail, key=lambda a: abs(a - w))


def base(t):
    ft = ForceTime(Input=rd); ft.ForcedTime = t; ft.IgnorePipelineTime = 1
    return MergeBlocks(Input=ft)


def lowest_metal(src, t):
    """lowest y of metal (alpha > 0.5) above the receiver, from the mid-plane slice."""
    sl = Slice(Input=src); sl.SliceType = 'Plane'; sl.SliceType.Origin = [0, 0, 0]; sl.SliceType.Normal = [0, 0, 1]
    UpdatePipeline(proxy=sl, time=t)
    from vtk import vtkCellCenters
    o = sm.Fetch(sl); cc = vtkCellCenters(); cc.SetInputData(o); cc.Update()
    P = np.asarray(dsa.WrapDataObject(cc.GetOutput()).Points); a = np.asarray(dsa.WrapDataObject(o).CellData['alpha.metal'])
    m = (a > 0.5) & (P[:, 1] > 8.1e-6) & (np.abs(P[:, 0]) < R)
    Delete(sl)
    return P[m, 1].min() if m.any() else 19.9e-6


def section(t, y0, y1, r, stretch):
    """2-D mirrored section (cell data), scaled to um with y stretched."""
    src = base(t)
    for org, nrm, inv in (([0, y0, 0], [0, 1, 0], 0), ([0, y1, 0], [0, 1, 0], 1), ([r, 0, 0], [1, 0, 0], 1)):
        c = Clip(Input=src); c.ClipType = 'Plane'; c.ClipType.Origin = org; c.ClipType.Normal = nrm; c.Invert = inv; src = c
    sl = Slice(Input=src); sl.SliceType = 'Plane'; sl.SliceType.Origin = [0, 0, 0]; sl.SliceType.Normal = [0, 0, 1]
    c1 = Calculator(Input=sl); c1.AttributeType = 'Cell Data'; c1.ResultArrayName = 'Umag'; c1.Function = 'mag(U)'
    c2 = Calculator(Input=c1); c2.AttributeType = 'Cell Data'; c2.ResultArrayName = 'Ur'; c2.Function = 'U_X'
    c3 = Calculator(Input=c2); c3.AttributeType = 'Cell Data'; c3.ResultArrayName = 'Uy'; c3.Function = 'U_Y'
    # Ur is computed before mirroring, so on both halves it stays "positive = away from the axis"
    ref = Reflect(Input=c3); ref.Plane = 'X'; ref.Center = 0.0; ref.CopyInput = 1
    tr = Transform(Input=ref); tr.Transform = 'Transform'; tr.Transform.Scale = [1e6, 1e6 * stretch, 1e6]
    UpdatePipeline(proxy=tr, time=t)
    return tr


def view():
    v = CreateView('RenderView'); v.ViewSize = [W, H]; v.Background = [1, 1, 1]
    try: v.UseColorPaletteForBackground = 0
    except AttributeError: pass
    v.OrientationAxesVisibility = 0
    return v


def camera(v, y0, y1, r, stretch):
    Render(v); v.ResetCamera()
    yc = 0.5 * (y0 + y1) * 1e6 * stretch
    v.CameraPosition = [0, yc, 100]; v.CameraFocalPoint = [0, yc, 0]; v.CameraViewUp = [0, 1, 0]
    v.CameraParallelProjection = 1
    v.CameraParallelScale = max(0.5 * (y1 - y0) * 1e6 * stretch * 1.18, 1.06 * r * 1e6 * H / W)


def captions(v, lines):
    for txt, pos, fs in lines:
        tx = Text(Text=txt); td = Show(tx, v); td.WindowLocation = 'Any Location'; td.Position = pos; td.FontSize = fs; td.Color = [0, 0, 0]


def contour_line(v, tr, color=(0, 0, 0), width=2.0):
    ct = Contour(Input=CellDatatoPointData(Input=tr)); ct.ContourBy = ['POINTS', 'alpha.metal']; ct.Isosurfaces = [0.5]
    cd = Show(ct, v); ColorBy(cd, None); cd.DiffuseColor = list(color); cd.AmbientColor = list(color); cd.LineWidth = width


def field_on(v, tr, field, rng, preset, label, log=False):
    d = Show(tr, v); ColorBy(d, ('CELLS', field))
    lut = GetColorTransferFunction(field); lut.ApplyPreset(preset, True)
    lut.UseLogScale = 1 if log else 0
    lut.RescaleTransferFunction(*rng); lut.AutomaticRescaleRangeMode = 'Never'
    d.Ambient = 1.0; d.Diffuse = 0.0; d.SetScalarBarVisibility(v, True)
    bar = GetScalarBar(lut, v); bar.Title = label; bar.ComponentTitle = ''
    bar.TitleColor = [0, 0, 0]; bar.LabelColor = [0, 0, 0]; bar.TitleFontSize = 16; bar.LabelFontSize = 14
    return d


def save(v, fname):
    Render(v); out = os.path.join(OUT, fname)
    SaveScreenshot(out, v, ImageResolution=[W, H]); print('saved', out, flush=True)
    Delete(v)


def frame_for(t):
    ylow = lowest_metal(base(t), t)
    y0 = min(19.35e-6, ylow - 0.25e-6)
    stretch = min(2.0, (8.0 * H / W - 0.25) / ((Y1 - y0) * 1e6))
    return y0, stretch


FOOT = 'section through the beam axis (r = 0 at centre), width 8 um; y stretched x%g   |   black line: metal fraction = 0.5   |   glass side up, receiver below'

# ---------- A: sections ----------
for w in ([] if os.environ.get('ONLY3D') else TIMES):
    t = nearest(w); tag = '%05.0fps' % (t * 1e12)
    y0, st = frame_for(t)
    tr = section(t, y0, Y1, R, st)
    # Image 1: alpha, full field
    v = view(); field_on(v, tr, 'alpha.metal', (0, 1), 'Viridis (matplotlib)', 'metal volume fraction')
    contour_line(v, tr, (1, 1, 1), 1.5)
    camera(v, y0, Y1, R, st)
    captions(v, [('%s   t = %.0f ps   metal volume fraction (0 = gas, 1 = metal)' % (NAME, t * 1e12), [0.02, 0.93], 22),
                 (FOOT.replace('black line', 'white line') % round(st, 2), [0.02, 0.02], 14)])
    save(v, 'alpha_%s.png' % tag)
    # Image 2: |U| with arrows (gas and metal)
    v = view(); field_on(v, tr, 'Umag', (0, 250), 'Turbo', 'speed (m/s)')
    pts = CellDatatoPointData(Input=tr)
    g = Glyph(Input=pts, GlyphType='Arrow'); g.OrientationArray = ['POINTS', 'U']; g.ScaleArray = ['POINTS', 'No scale array']
    g.ScaleFactor = 0.18; g.GlyphMode = 'Uniform Spatial Distribution (Bounds Based)'; g.MaximumNumberOfSamplePoints = 900
    gd = Show(g, v); ColorBy(gd, None); gd.DiffuseColor = [1, 1, 1]; gd.AmbientColor = [1, 1, 1]
    contour_line(v, tr)
    camera(v, y0, Y1, R, st)
    captions(v, [('%s   t = %.0f ps   speed (m/s), arrows = flow direction (gas and metal)' % (NAME, t * 1e12), [0.02, 0.93], 22),
                 (FOOT % round(st, 2), [0.02, 0.02], 14)])
    save(v, 'Umag_%s.png' % tag)
    # Image 3: radial velocity
    v = view(); field_on(v, tr, 'Ur', (-100, 100), 'Cool to Warm', 'radial velocity (m/s)')
    contour_line(v, tr)
    camera(v, y0, Y1, R, st)
    captions(v, [('%s   t = %.0f ps   radial velocity (m/s): blue = towards the axis, red = outwards' % (NAME, t * 1e12), [0.02, 0.93], 22),
                 (FOOT % round(st, 2), [0.02, 0.02], 14)])
    save(v, 'Ur_%s.png' % tag)
    # Image 4: pressure (range from the 2-98 % percentiles of the section)
    pv = np.asarray(dsa.WrapDataObject(sm.Fetch(tr)).CellData['p'])
    lo, hi = np.percentile(pv, 2), np.percentile(pv, 98)
    if hi - lo < 1e3: lo, hi = lo - 5e2, hi + 5e2
    v = view(); field_on(v, tr, 'p', (lo, hi), 'Cool to Warm', 'pressure (Pa)')
    contour_line(v, tr)
    camera(v, y0, Y1, R, st)
    captions(v, [('%s   t = %.0f ps   pressure (Pa), colour range = 2-98 %% of the section' % (NAME, t * 1e12), [0.02, 0.93], 22),
                 (FOOT % round(st, 2), [0.02, 0.02], 14)])
    save(v, 'p_%s.png' % tag)

# ---------- B: laser-heating close-ups (film only) ----------
y0e, y1e, re_, ste = 19.985e-6, 20.085e-6, 4.0e-6, 15.0
for w in ([] if os.environ.get('ONLY3D') else EARLY):
    t = nearest(w); tag = '%07.0ffs' % (t * 1e15)
    tr = section(t, y0e, y1e, re_, ste)
    for field, rng, preset, label, log in (('Te', (300, 20000), 'Inferno (matplotlib)', 'electron temperature (K)', False),
                                           ('Tl', (300, 15000), 'Inferno (matplotlib)', 'lattice (atom) temperature (K)', False),
                                           ('Q_laser', (1e18, 1e23), 'Plasma (matplotlib)', 'laser source (W/m^3)', True)):
        if field == 'Q_laser' and w != EARLY[0]:
            continue   # the pulse is over after ~0.5 ps; later frames are blank
        v = view(); field_on(v, tr, field, rng, preset, label, log)
        contour_line(v, tr, (1, 1, 1), 1.5)
        camera(v, y0e, y1e, re_, ste)
        captions(v, [('%s   t = %.1f ps   %s in the 71 nm film (laser enters from the glass side, top)' % (NAME, t * 1e12, label), [0.02, 0.93], 22),
                     ('film close-up, width 8 um; y stretched x%g   |   white line: metal fraction = 0.5' % ste, [0.02, 0.02], 14)])
        save(v, 'heat_%s_%s.png' % (field, tag))

# ---------- C: 3-D views (alpha = 0.5 surface revolved about the axis) ----------
for w in [float(x) for x in os.environ.get('TIMES3D', '1e-9,6e-9,1.34e-8,1.98e-8').split(',')]:
    t = nearest(w); tag = '%05.0fps' % (t * 1e12)
    src = base(t)
    c = Clip(Input=src); c.ClipType = 'Plane'; c.ClipType.Origin = [3.5e-6, 0, 0]; c.ClipType.Normal = [1, 0, 0]; c.Invert = 1
    c2 = Clip(Input=c); c2.ClipType = 'Plane'; c2.ClipType.Origin = [0, 8.2e-6, 0]; c2.ClipType.Normal = [0, 1, 0]; c2.Invert = 0
    sl = Slice(Input=c2); sl.SliceType = 'Plane'; sl.SliceType.Origin = [0, 0, 0]; sl.SliceType.Normal = [0, 0, 1]
    ca = Calculator(Input=sl); ca.AttributeType = 'Cell Data'; ca.ResultArrayName = 'Uy'; ca.Function = 'U_Y'
    ct = Contour(Input=CellDatatoPointData(Input=ca)); ct.ContourBy = ['POINTS', 'alpha.metal']; ct.Isosurfaces = [0.5]
    ct.ComputeScalars = 1
    # (r, y, 0) -> (r, 0, y) in um, then sweep about the z axis
    tf = Transform(Input=ct); tf.Transform = 'Transform'; tf.Transform.Rotate = [90, 0, 0]; tf.Transform.Scale = [1e6, 1e6, 1e6]
    # 270 degrees: a cut-away that shows the profile; the camera looks into the missing quarter
    rx = RotationalExtrusion(Input=tf); rx.Resolution = 135; rx.Angle = 270.0
    UpdatePipeline(proxy=rx, time=t)
    v = CreateView('RenderView'); v.ViewSize = [W, H]; v.Background = [1, 1, 1]
    try: v.UseColorPaletteForBackground = 0
    except AttributeError: pass
    v.OrientationAxesVisibility = 0
    d = Show(rx, v); ColorBy(d, ('POINTS', 'Uy'))
    lut = GetColorTransferFunction('Uy'); lut.ApplyPreset('Cool to Warm', True); lut.RescaleTransferFunction(-200, 200)
    lut.AutomaticRescaleRangeMode = 'Never'; d.SetScalarBarVisibility(v, True)
    bar = GetScalarBar(lut, v); bar.Title = 'vertical velocity (m/s)'; bar.ComponentTitle = ''; bar.TitleColor = [0, 0, 0]; bar.LabelColor = [0, 0, 0]
    d.Ambient = 0.25; d.Diffuse = 0.85; d.Specular = 0.2
    Render(v); v.ResetCamera()
    v.CameraViewUp = [0, 0, 1]; v.CameraFocalPoint = [0, 0, 18.9]
    v.CameraPosition = [11.0, -11.0, 15.5]; v.CameraViewAngle = 28
    Render(v)
    captions(v, [('%s   t = %.0f ps   2-D AXISYMMETRIC result drawn in 3-D (metal surface revolved about the axis) - NOT a 3-D simulation' % (NAME, t * 1e12), [0.02, 0.93], 20),
                 ('the solution is the same at every angle by construction (5 deg wedge); non-axisymmetric break-up (rings into droplets) cannot appear   |   270 deg cut-away, r < 3.5 um, colour = vertical velocity, true scale', [0.02, 0.02], 14)])
    save(v, '3d_%s.png' % tag)
