"""16:9 frames of the whole-run animation (laser pulse to impact).

    python3 tools/compose_full_run.py <paneldir> <zoomdir> <curves.csv> <outdir> <dpi>

<paneldir>, <zoomdir>: output of tools/render_full_run.py for the main section and the zoom, same times.
<curves.csv>: output of tools/full_run_curves.py for the same times. dpi 100 gives 1920x1080, 200 gives 3840x2160.
RIGID=1: the receiver (y < 8 um) is drawn as a solid wall on every frame and the caption says that the
impact is on a rigid wall (frames after 21.5 ns from the rigid-receiver run, docs/INVESTIGATION_LOG.md section 16).
FILM=<csv>: output of tools/wall_film_profile.py; the metal film on the wall that is thinner than half a cell
(so missing from the alpha = 0.5 rendering) is drawn in cyan at its equivalent thickness, true scale, in
the section and the zoom, behind the rendered metal.
Encode with ffmpeg -framerate 20 -i <outdir>/c_%04d.png -c:v libx264 -pix_fmt yuv420p -crf 18 out.mp4
"""
import sys, os, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib import cm, colors
P, Z, F, O, DPI = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]); os.makedirs(O, exist_ok=True)
WUM, HUM = 10.0, 13.6           # section window (um)
ZW, ZH = 3.0, 2.4               # zoom window (um)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 15})
T = np.atleast_2d(np.loadtxt(os.path.join(P, 'times.txt')))
d = np.genfromtxt(F, delimiter=',', names=True)
RIGID = os.environ.get('RIGID', '0') == '1'
tc = d['t_ns']; gap = np.where(d['gap_um'] < 0.12, 0.0, d['gap_um']); va = d['v_axis_ms']
def stage(tns):
    if tns < 0.002: return 'laser heats the film (electrons, then lattice)'
    if tns < 0.3: return 'trapped vapour pushes the hot centre forward'
    if tns < 4.1: return 'film stretches into a dome'
    if tns < TC: return 'dome wall breaks; a thin jet flies to the receiver'
    return 'jet hits the rigid receiver and spreads' if RIGID else 'jet hits the receiver'
FILM = None
if os.environ.get('FILM'):
    _f = np.loadtxt(os.environ['FILM'], delimiter=',', skiprows=1, ndmin=2)
    FILM = (_f[:, 0], np.array([float(v) for v in open(os.environ['FILM']).readline().strip().split(',')[1:]]), _f[:, 1:])
BG = np.array([0.86, 0.86, 0.87])
def see_through(img):
    # background pixels of a rendered panel transparent, so that the film drawn behind shows
    a = np.ones(img.shape[:2]); a[np.abs(img[..., :3] - BG).sum(-1) < 0.03] = 0.0
    return np.dstack([img[..., :3], a])
def film(ax, tns, w, h, span, height):
    # draw the wall film (um) behind the panel image; span, height: panel window in um
    if FILM is None: return
    k = np.argmin(np.abs(FILM[0] - tns))
    if abs(FILM[0][k] - tns) > 0.02: return
    R, hh = FILM[1], FILM[2][k] * 1e-3
    if hh.max() < 0.002: return
    rr = np.concatenate([-R[::-1], R]); hv = np.concatenate([hh[::-1], hh])
    x = w/2 + rr/span*w; yw = h*(1 - 1.0/height); yt = yw - hv/height*h
    ax.fill_between(x, yw, yt, where=hv > 0.002, color='#06b6d4', lw=0, zorder=1)
    ax.plot(np.where(hv > 0.002, x, np.nan), yt, color='#06b6d4', lw=1.2, zorder=1)
TC = tc[np.argmax(gap == 0.0)] if (gap == 0.0).any() else 1e9  # first frame in contact
def cbar(fig, rect, cmap, norm, lab, ext, ticks=None, tl=None):
    cax = fig.add_axes(rect)
    cb = fig.colorbar(cm.ScalarMappable(norm, cmap), cax=cax, orientation='horizontal', extend=ext)
    if ticks is not None: cb.set_ticks(ticks); cb.set_ticklabels(tl); cb.minorticks_off()
    cb.ax.tick_params(labelsize=12.5); cax.set_title(lab, fontsize=13.5, loc='left', pad=4)
def halves(dirn, i):
    iT = plt.imread(os.path.join(dirn, 'T_%04d.png' % i)); iU = plt.imread(os.path.join(dirn, 'U_%04d.png' % i))
    w_ = iT.shape[1]; return np.concatenate([iT[:, :w_//2], iU[:, w_//2:]], axis=1)
n = 0
for i, t in T:
    i = int(i); tns = t*1e9
    fig = plt.figure(figsize=(19.2, 10.8), dpi=DPI); fig.patch.set_facecolor('white')
    fig.text(0.025, 0.94, 'Femtosecond LIFT of titanium: from the laser pulse to impact', fontsize=29, weight='bold', color='#1a1a1a')
    fig.text(0.025, 0.902, '71 nm Ti on glass  ·  343 nm, 214 fs, 332 mJ/cm²  ·  12 µm gap  ·  OpenFOAM VOF solver with a two-temperature model, evaporation and a trapped vapour pocket', fontsize=15, color='0.35')
    ts = ('t = %5.2f ns' % tns) if tns >= 1 else (('t = %4.0f ps' if tns >= 0.01 else 't = %4.1f ps') % (tns*1e3))
    fig.text(0.975, 0.94, ts, fontsize=29, ha='right', family='DejaVu Sans Mono', color='#1a1a1a')
    # 1. section, true scale
    img = halves(P, i); h_, w_ = img.shape[:2]
    H = 0.70; Wd = H*w_/h_*10.8/19.2; L, y0 = 0.025, 0.16
    ax = fig.add_axes([L, y0, Wd, H]); ax.axis('off')
    if FILM is not None:
        ax.add_patch(plt.Rectangle((-0.5, -0.5), w_, h_, fc=BG, ec='none', zorder=0)); film(ax, tns, w_, h_, WUM, HUM)
        img = see_through(img)
    ax.imshow(img, interpolation='lanczos', zorder=2)
    ax.axvline(w_//2, color='white', lw=1.5, alpha=0.9)
    ax.text(0.015, 0.985, 'temperature', transform=ax.transAxes, fontsize=13.5, color='0.25', va='top')
    ax.text(0.985, 0.985, 'speed', transform=ax.transAxes, fontsize=13.5, color='0.25', va='top', ha='right')
    if RIGID: ax.add_patch(plt.Rectangle((-0.5, h_*(1-1.0/HUM)), w_, h_/HUM, fc='#4b5563', ec='#9ca3af', lw=0, hatch='///', zorder=4))
    ax.text(0.2*w_, 0.975*h_, 'rigid receiver (300 K)' if RIGID else 'receiver (Ti, 300 K)', fontsize=12.5, color='white', ha='center', va='center', zorder=5)
    sb = 2.0/WUM*w_; ax.plot([0.05*w_, 0.05*w_+sb], [0.80*h_, 0.80*h_], color='0.2', lw=4, solid_capstyle='butt')
    ax.text(0.05*w_+sb/2, 0.79*h_, '2 µm', ha='center', va='bottom', fontsize=13, color='0.2')
    zy0 = (7.0 + ZH - 7.0)/HUM  # zoom frame outline (window y 7.0-9.4 um, |r| < 1.5 um)
    ax.add_patch(plt.Rectangle(((WUM/2-ZW/2)/WUM*w_, (1-ZH/HUM)*h_), ZW/WUM*w_, ZH/HUM*h_, fill=False, ec='#2563eb', lw=1.6, ls='--'))
    ax.set_xlim(-0.5, w_-0.5); ax.set_ylim(h_-0.5, -0.5)  # patches must not widen the image limits
    fig.text(L, y0+H+0.008, 'glass side; laser from above ↑', fontsize=12.5, color='0.4')
    # 2. 3-D view and zoom
    i3 = plt.imread(os.path.join(P, 'D_%04d.png' % i)); x2 = L + Wd + 0.025; W2 = 0.30; h3 = 0.33; w3 = h3*i3.shape[1]/i3.shape[0]*10.8/19.2
    ax = fig.add_axes([x2+(W2-w3)/2, 0.865-h3, w3, h3]); ax.imshow(i3, interpolation='lanczos'); ax.axis('off')
    fig.text(x2+W2/2, 0.865-h3+0.004, 'metal surface in 3-D (the axisymmetric\nresult swept round the axis); receiver in grey', fontsize=11.5, color='0.35', ha='center', va='top')
    iz = halves(Z, i); hz, wz = iz.shape[:2]; hzf = 0.30; wzf = hzf*wz/hz*10.8/19.2
    ax = fig.add_axes([x2 + (W2-wzf)/2, y0, wzf, hzf])
    if FILM is not None:
        ax.add_patch(plt.Rectangle((-0.5, -0.5), wz, hz, fc=BG, ec='none', zorder=0)); film(ax, tns, wz, hz, ZW, ZH)
        iz = see_through(iz)
        if tns > 21.5: ax.text(0.97*wz, 0.03*hz, 'cyan: metal film thinner than\none mesh cell, drawn at its\nmean thickness', fontsize=10.5, color='#0e7490', ha='right', va='top', zorder=5)
    ax.imshow(iz, interpolation='lanczos', zorder=2)
    if RIGID: ax.add_patch(plt.Rectangle((-0.5, hz*(1-1.0/ZH)), wz, hz/ZH, fc='#4b5563', ec='#9ca3af', lw=0, hatch='///', zorder=4))
    ax.set_xticks([]); ax.set_yticks([]); [s.set_color('#2563eb') for s in ax.spines.values()]; [s.set_linestyle('--') for s in ax.spines.values()]; [s.set_linewidth(1.6) for s in ax.spines.values()]
    sb = 0.5/ZW*wz; ax.plot([0.05*wz, 0.05*wz+sb], [0.12*hz, 0.12*hz], color='0.2', lw=4, solid_capstyle='butt')
    ax.text(0.05*wz+sb/2, 0.10*hz, '0.5 µm', ha='center', va='bottom', fontsize=12, color='0.2')
    ax.set_xlim(-0.5, wz-0.5); ax.set_ylim(hz-0.5, -0.5)
    fig.text(x2+W2/2, y0+hzf+0.006, 'zoom on the landing site (blue box)', fontsize=12.5, color='#2563eb', ha='center')
    # 3. stage, curves, colour bars, caption
    x3 = x2 + W2 + 0.03; w4 = 0.975 - x3
    fig.text(x3, 0.83, stage(tns), fontsize=17, color='#9a3412', weight='bold', wrap=True)
    for k, (y, lab, col, top) in enumerate(((gap, 'gap between the lowest metal and the receiver (µm)', '#d94801', 12.5),
                                           (va, 'speed of the metal on the beam axis (m/s)', '#238b45', 520))):
        ax = fig.add_axes([x3+0.01, 0.62-0.17*k, w4-0.01, 0.10])
        m = tc <= tns*(1+1e-9)
        ax.plot(tc, y, color='0.85', lw=2.5); ax.plot(tc[m], y[m], color=col, lw=3.2)
        ax.plot([tns], [np.interp(tns, tc, y)], 'o', color=col, ms=9)
        ax.set_xscale('log'); ax.set_xlim(1e-3, 33); ax.set_ylim(0, top)
        ax.set_xticks([1e-3, 1e-2, 1e-1, 1, 10]); ax.set_xticklabels(['1 ps', '10 ps', '0.1 ns', '1 ns', '10 ns'] if k else [])
        ax.set_title(lab, fontsize=13, loc='left'); ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=11.5)
    cbar(fig, [x3+0.01, 0.335, w4-0.02, 0.014], 'inferno', colors.LogNorm(300, 10000), 'lattice temperature (K, log scale)', 'both',
         [300, 1000, 1941, 3000, 10000], ['300', '1000', 'Ti melts\n1941', '3000', '10 000'])
    cbar(fig, [x3+0.01, 0.245, w4-0.02, 0.014], 'viridis', colors.Normalize(0, 500), 'speed (m/s)', 'max')
    fig.text(x3, 0.175, ('Cross-section through the beam axis, mirrored, true scale.\n'
             'Simulation: 2-D axisymmetric, coarse mesh (62k cells); not compared\nwith experiment yet. '
             'After about 7 ns the jet is thinner than a\nmesh cell, so flight and impact are qualitative. ' +
             ('The receiver is a\nrigid 300 K wall (no slip, 90° contact angle); the jet does not\nfreeze in the model.' if RIGID else
              'The model treats\nall metal as liquid, including the 300 K receiver, so the crater\nat impact is not physical.')),
             fontsize=11.5, color='0.35', va='top', linespacing=1.4)
    fig.savefig(os.path.join(O, 'c_%04d.png' % n), dpi=DPI); plt.close(fig); n += 1
print('composed', n)
