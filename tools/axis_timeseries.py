#!/usr/bin/env python3
"""Film position and velocity on the beam axis vs time, from postProcessing/axisLine.

Usage: python3 tools/axis_timeseries.py cases/<case>
Writes tools/axis_timeseries_<case>.csv and figures/<case>/axis_position_speed_temperature.png.
Film = the largest contiguous metal segment (alpha>0.5) on the axis between
Y_MIN (default 8.05 um, just above the receiver) and 21 um, so the film is followed
however far it travels across the gap.
The last CSV column gives the number of metal segments in the window (>1: detached
pieces or a smeared film, check in ParaView).
"""
import glob, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
Y_MIN = float(os.environ.get("Y_MIN", "8.05e-6"))
Y_MAX = 21e-6
CASE = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "case")
rows = []
for d in sorted(glob.glob(os.path.join(CASE, "postProcessing", "axisLine", "*")), key=lambda x: float(os.path.basename(x))):
    f = glob.glob(os.path.join(d, "axis_*.xy"))
    if not f:
        continue
    names = os.path.basename(f[0])[5:-3].split("_")
    try:
        A = np.loadtxt(f[0])
    except Exception:
        continue
    if A.ndim != 2:
        continue
    # column of each field in the raw file (y first; vectors take 3 columns). OpenFOAM
    # writes scalars before vectors, so the file name order, not the dictionary order,
    # gives the layout: e.g. axis_Tl_alpha.metal_p_U.xy -> y Tl alpha p Ux Uy Uz
    col, c = {}, 1
    for n in names:
        col[n] = c; c += 3 if n == "U" else 1
    ia, it, iu = col["alpha.metal"], col["Tl"], col["U"]   # U_y is column iu + 1
    y, a, T, uy = A[:, 0], A[:, ia], A[:, it], A[:, iu + 1]
    inw = (y > Y_MIN) & (y < Y_MAX)
    met = (a > 0.5) & inw
    if not met.any():
        continue
    # contiguous metal segments (in sample order) inside the window
    idx = np.flatnonzero(met)
    segs = np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1)
    big = max(segs, key=len)
    m = np.zeros_like(met); m[big] = True
    # alpha-weighted centre over the samples around the film (sub-sample accuracy)
    lo, hi = y[big].min() - 0.1e-6, y[big].max() + 0.1e-6
    w = inw & (y > lo) & (y < hi)
    ycm = (y[w] * a[w]).sum() / max(a[w].sum(), 1e-30)
    rows.append((float(os.path.basename(d)), y[m].min(), y[m].max(), uy[m].mean(), T[m].max(), int(m.sum()), ycm, len(segs)))
if not rows:
    sys.exit("no samples")
R = np.array(rows)
TAG = os.path.basename(os.path.normpath(CASE))
HDR = ("time [s], lowest point of the largest metal piece on the axis [m], its highest point [m], "
       "mean vertical velocity of that piece [m/s] (negative = towards the receiver), "
       "highest lattice temperature in it [K], number of axis samples in it, "
       "metal-fraction-weighted centre of the film on the axis [m], number of separate metal pieces on the axis")
np.savetxt(os.path.join(HERE, "axis_timeseries_%s.csv" % TAG), R, delimiter=",", header=HDR)
os.makedirs(os.path.join(HERE, "..", "figures", TAG, "data"), exist_ok=True)
np.savetxt(os.path.join(HERE, "..", "figures", TAG, "data", "axis_timeseries.csv"), R, delimiter=",", header=HDR)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
t = R[:, 0] * 1e9   # ns
ax[0].plot(t, -(R[:, 6] - R[0, 6]) * 1e6)
ax[0].set_ylabel("distance moved towards the receiver (um)")
ax[0].set_title("Film centre on the beam axis", fontsize=10)
ax[1].plot(t, -R[:, 3]); ax[1].axhline(71.3, color="grey", ls=":", lw=1)
ax[1].text(t[-1], 71.3, "speed at which inertia equals surface tension (71 m/s)  ", va="bottom", ha="right", fontsize=8)
ax[1].set_ylabel("speed towards the receiver (m/s)")
ax[1].set_title("Speed of the metal on the beam axis", fontsize=10)
ax[2].plot(t, R[:, 4]); ax[2].axhline(3560, color="grey", ls=":", lw=1)
ax[2].text(t[-1], 3560, "boiling temperature of titanium (3560 K)  ", va="bottom", ha="right", fontsize=8)
ax[2].set_ylabel("temperature (K)")
ax[2].set_title("Highest atom (lattice) temperature in the film on the axis", fontsize=10)
for a_ in ax:
    a_.set_xlabel("time after the laser pulse (ns)"); a_.grid(alpha=0.25)
fig.suptitle("History on the beam axis, case %s" % TAG, fontsize=11)
fig.tight_layout()
out = os.path.join(HERE, "..", "figures", TAG, "axis_position_speed_temperature.png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=130)
print("saved", out, "last t=%.4g ns, speed %.4g m/s, displacement %.4g nm" % (t[-1], -R[-1, 3], -(R[-1, 6] - R[0, 6]) * 1e9))
