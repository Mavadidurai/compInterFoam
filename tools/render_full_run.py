"""Panels for the whole-run animation (laser pulse to impact) of a decomposed axisymmetric wedge case.

    LP_NUM_THREADS=1 CAM3='0,0.3,24,-24,9,21' xvfb-run -a pvpython tools/render_full_run.py \
        <case> <outdir> <x_half_um> <y0_um> <y1_um> <SY> <W> [times|all]

Writes into <outdir>, for each time i:
  T_%04d.png  section through the axis, mirrored, metal = alpha.metal >= 0.5, coloured by the
              lattice temperature Tl (log scale 300-10,000 K)
  U_%04d.png  the same coloured by speed (0-UMAX m/s, default 500)
  D_%04d.png  the alpha = 0.5 surface above y = CLIP3 um (default 8.02) revolved 270 deg about the
              axis, coloured by Tl, with the receiver drawn as a grey disk at y = 8 um
and times.txt (index, time). The full run used <x_half_um> 5, y 7.0-20.6 um, W 1000 for the main
section and 1.5, 7.0-9.4 um, W 900 for the zoom on the landing site; then tools/compose_full_run.py.
Run it with LP_NUM_THREADS=1 next to a running MPI case: the software renderer otherwise takes every
core and the solver ranks, which busy-wait on each other, slow down about five-fold.
"""
import os, sys
from paraview.simple import *
case, out = sys.argv[1], sys.argv[2]; XH, Y0, Y1 = [float(v)*1e-6 for v in sys.argv[3:6]]; SY = float(sys.argv[6]); W = int(sys.argv[7])
sel = sys.argv[8:] or ['all']
os.makedirs(out, exist_ok=True)
rd = OpenFOAMReader(FileName=os.path.join(case, 'case.foam'))
rd.CaseType = 'Decomposed Case'; rd.MeshRegions = ['internalMesh']; rd.CellArrays = ['alpha.metal', 'U', 'Tl']
rd.UpdatePipelineInformation()
allt = [t for t in rd.TimestepValues if t > 0]
times = allt if sel == ['all'] else [min(allt, key=lambda a: abs(a-float(s))) for s in sel]
mb = MergeBlocks(Input=rd)
sl = Slice(Input=mb); sl.SliceType.Normal = [0, 0, 1]; sl.SliceType.Origin = [0, 0, 0]
pd = CellDatatoPointData(Input=sl)
um = Calculator(Input=pd); um.ResultArrayName = 'Umag'; um.Function = 'mag(U)'
# section: mirrored, stretched, metal region bounded by the alpha = 0.5 iso-line
rf = Reflect(Input=um); rf.Plane = 'X'; rf.CopyInput = 1
tr = Transform(Input=rf); tr.Transform.Scale = [1, SY, 1]
cl = Clip(Input=tr); cl.ClipType = 'Scalar'; cl.Scalars = ['POINTS', 'alpha.metal']; cl.Value = 0.5; cl.Invert = 0
H = int(round(W*SY*(Y1-Y0)/(2*XH))); cy = 0.5*SY*(Y0+Y1)
def cam(v):
    v.CameraPosition = [0, cy, 1e-4]; v.CameraFocalPoint = [0, cy, 0]; v.CameraViewUp = [0, 1, 0]; v.CameraParallelProjection = 1
    v.CameraParallelScale = 0.5*SY*(Y1-Y0)
def newview(w, h, bg):
    v = CreateView('RenderView'); v.ViewSize = [w, h]; v.OrientationAxesVisibility = 0
    try: v.UseColorPaletteForBackground = 0
    except Exception: pass
    v.Background = bg; return v
def lut(field, preset, rng):
    l = GetColorTransferFunction(field); l.ApplyPreset(preset, True); l.RescaleTransferFunction(*rng); l.AutomaticRescaleRangeMode = 'Never'; return l
RNG = {'Tl': ('Inferno (matplotlib)', (300, 10000), True), 'U': ('Viridis (matplotlib)', (0, float(os.environ.get('UMAX','500'))), False)}
def fix_luts():
    for f, (p, r, lg) in RNG.items():
        l = lut(f, p, r)
        l.UseLogScale = 1 if lg else 0
        l.RescaleTransferFunction(*r)
        if f == 'U': l.VectorMode = 'Magnitude'
views = {}
for name, field in (('T', ('POINTS', 'Tl')), ('U', ('POINTS', 'U', 'Magnitude'))):
    v = newview(W, H, [0.86, 0.86, 0.87]); d = Show(cl, v); ColorBy(d, field); d.SetScalarBarVisibility(v, False); fix_luts()
    v.ViewTime = times[0]; Render(v); cam(v); views[name] = v
# 3-D: alpha = 0.5 line of the unmirrored half-section, (r, y) -> (r, 0, y) in um, swept about z
ct = Contour(Input=um); ct.ContourBy = ['POINTS', 'alpha.metal']; ct.Isosurfaces = [0.5]; ct.ComputeScalars = 1
c3 = Clip(Input=ct); c3.ClipType = 'Plane'; c3.ClipType.Origin = [0, float(os.environ.get('CLIP3','8.02'))*1e-6, 0]; c3.ClipType.Normal = [0, 1, 0]; c3.Invert = 0
c4 = Clip(Input=c3); c4.ClipType = 'Plane'; c4.ClipType.Origin = [XH, 0, 0]; c4.ClipType.Normal = [1, 0, 0]; c4.Invert = 1
tf = Transform(Input=c4); tf.Transform.Rotate = [90, 0, 0]; tf.Transform.Scale = [1e6, 1e6, 1e6]
es = ExtractSurface(Input=tf)
rx = RotationalExtrusion(Input=es); rx.Resolution = 180; rx.Angle = 270.0
W3, H3 = int(W*0.5), int(W*0.4)
v3 = newview(W3, H3, [1, 1, 1]); d3 = Show(rx, v3); ColorBy(d3, ('POINTS', 'Tl')); d3.SetScalarBarVisibility(v3, False); fix_luts()
dk = Disk(); dk.InnerRadius = 0.0; dk.OuterRadius = XH*1e6; dk.CircumferentialResolution = 180
dkt = Transform(Input=dk); dkt.Transform.Translate = [0, 0, 8.0]
dd = Show(dkt, v3); dd.ColorArrayName = [None, ""]; dd.DiffuseColor = [0.72, 0.74, 0.78]; dd.AmbientColor = [0.72, 0.74, 0.78]; dd.Opacity = 1.0
d3.Ambient = 0.2; d3.Diffuse = 0.85; d3.Specular = 0.35; d3.SpecularPower = 40
zc = 0.5*(Y0+Y1)*1e6
def cam3(v):
    v.CameraViewUp = [0, 0, 1]; C=[float(c) for c in os.environ.get('CAM3','0,-0.3,7.2,-7.2,5.0,30').split(',')]; v.CameraFocalPoint = [0, 0, zc + C[1]]; v.CameraPosition = [C[2], C[3], zc + C[4]]; v.CameraViewAngle = C[5]
v3.ViewTime = times[0]; UpdatePipeline(proxy=rx, time=times[0]); print('3-D bounds', rx.GetDataInformation().GetBounds(), flush=True)
Render(v3); v3.ResetCamera(); cam3(v3); Render(v3)
try: v3.EnableRenderOnInteraction = 0
except Exception: pass
with open(os.path.join(out, 'times.txt'), 'w') as fh:
    for i, t in enumerate(times):
        for n, v in views.items():
            v.ViewTime = t; fix_luts(); cam(v); Render(v); SaveScreenshot(os.path.join(out, '%s_%04d.png' % (n, i)), v, ImageResolution=[W, H])
        v3.ViewTime = t; fix_luts(); cam3(v3); Render(v3); SaveScreenshot(os.path.join(out, 'D_%04d.png' % i), v3, ImageResolution=[W3, H3])
        fh.write('%d %.9g\n' % (i, t)); print(i, '%.4g ps' % (t*1e12), flush=True)
