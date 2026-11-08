#!/usr/bin/env python3
"""Results package for a trapped-vapour run, optionally compared with a second run.

    python3 tools/results_package.py cases/<case> [cases/<comparison case>]

Reads only public run output:
  tools/axis_timeseries_<case>.csv (run tools/axis_timeseries.py first),
  the solver logs <case>/log.stage1..3 (or log.phase1/log.phase2 of older runs), with the
  .partN logs of interrupted stages first (vapourPocket lines, metal volume),
  <case>/droplet_detections.txt and run_status.txt.
Writes:
  figures/<case>/results_timeseries.png  (4 panels)
  figures/<case>/data/trapped_vapour_pressure.csv
  figures/<case>/results_table.md
"""
import os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.basename(os.path.normpath(sys.argv[1]))
CASEDIR = os.path.abspath(sys.argv[1])
REF = os.path.basename(os.path.normpath(sys.argv[2])) if len(sys.argv) > 2 else None
OUT = os.path.join(HERE, "..", "figures", CASE)
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from case_props import metal_density, surface_tension
RHO, SIGMA, H = metal_density(CASEDIR, 4515.0), surface_tension(CASEDIR, 1.64), 71.4e-9


def axis(case):
    f = os.path.join(HERE, "axis_timeseries_%s.csv" % case)
    if not os.path.exists(f):
        return None
    d = np.loadtxt(f, delimiter=",", ndmin=2)
    return {"t": d[:, 0] * 1e9, "v": -d[:, 3], "disp": (d[:, 6] - d[0, 6]) * 1e9, "Tl": d[:, 4],
            "nseg": d[:, 7] if d.shape[1] > 7 else np.ones(len(d))}


def parse_log(case):
    """pocket pressures and metal volume vs time from the solver log (streamed)."""
    pk, vol = [], []
    t_cur = None
    rx_t = re.compile(r"^Time = ([0-9.eE+-]+)")
    rx_p = re.compile(r"vapourPocket\(local\): t = ([0-9.eE+-]+) m_v = ([0-9.eE+-]+) kg, p\(0\) = ([0-9.eE+-]+) Pa, "
                      r"p\(1um\) = ([0-9.eE+-]+) Pa, p\(2um\) = ([0-9.eE+-]+) Pa.*?p_max = ([0-9.eE+-]+) Pa.*?d\(0\) = ([0-9.eE+-]+) m")
    rx_v = re.compile(r"^Metal volume: ([0-9.eE+-]+)")
    last_pt = -1.0
    import glob as _g
    def stage_logs(stem):  # log.<stem>.part1, part2, ... then log.<stem>
        parts = sorted(_g.glob(os.path.join(CASEDIR, "log.%s.part*" % stem)),
                       key=lambda x: int(re.sub(r"\D", "", x.rsplit("part", 1)[1]) or 0))
        last = os.path.join(CASEDIR, "log." + stem)
        return parts + ([last] if os.path.exists(last) else [])
    # tools/run_case.sh writes log.stage1..3; older runs wrote log.phase1/log.phase2
    files = [f for st in ("stage1", "stage2", "stage3", "phase1", "phase2") for f in stage_logs(st)]
    def lines():
        for fn in files:
            with open(fn, errors="replace") as f:
                for line in f:
                    yield line
    if True:
        for line in lines():
            m = rx_t.match(line)
            if m:
                t_cur = float(m.group(1)); continue
            if "vapourPocket(local)" in line:
                m = rx_p.search(line)
                if m:
                    t = float(m.group(1))
                    if t - last_pt >= 5e-12:            # thin to <= one sample per 5 ps
                        pk.append([float(x) for x in m.groups()]); last_pt = t
                continue
            m = rx_v.match(line)
            if m and t_cur is not None:
                vol.append((t_cur, float(m.group(1))))
    return np.array(pk), np.array(vol)


A, G = axis(CASE), (axis(REF) if REF else None)
PK, VOL = parse_log(CASE)
os.makedirs(os.path.join(OUT, "data"), exist_ok=True)
np.savetxt(os.path.join(OUT, "data", "trapped_vapour_pressure.csv"), PK, delimiter=",",
           header="time [s], vapour mass in the pocket, all radial bins [kg], pressure on the axis [Pa], "
                  "pressure at r = 1 um [Pa], pressure at r = 2 um [Pa], highest pressure of all bins [Pa], "
                  "film displacement away from the glass on the axis [m]")

# ---------- figure ----------
fig, ax = plt.subplots(2, 2, figsize=(13, 8.2))
c_i, c_g = "#EA580C", "#D4A017"
a = ax[0, 0]
if G is not None: a.plot(G["t"], G["v"], color=c_g, lw=2, label=REF)
a.plot(A["t"], A["v"], color=c_i, lw=2.5, label=CASE)
a.axhline(np.sqrt(SIGMA / (RHO * H)), color="grey", ls=":", lw=1); a.text(0.2, np.sqrt(SIGMA / (RHO * H)) + 4, "inertia equals surface tension (71 m/s)", fontsize=9, color="grey")
a.set_xlabel("time (ns)"); a.set_ylabel("speed towards the receiver (m/s)"); a.set_ylim(0, 230); a.legend(frameon=False); a.set_title("Speed of the metal on the beam axis")
a = ax[0, 1]
if G is not None: a.plot(G["t"], -G["disp"] / 1000, color=c_g, lw=2)
a.plot(A["t"], -A["disp"] / 1000, color=c_i, lw=2.5)
a.set_xlabel("time (ns)"); a.set_ylabel("distance moved towards the receiver (µm)"); a.set_title("Distance moved by the film on the beam axis")
a = ax[1, 0]
if len(PK):
    tp = PK[:, 0] * 1e9
    for k, lab, col in ((2, "on the beam axis", "#6D28D9"), (3, "1 µm from the axis", "#0F8F8C"), (4, "2 µm from the axis", "#7C8798"), (5, "highest anywhere", "#C62828")):
        y = np.where(PK[:, k] > 0, PK[:, k], np.nan)
        a.semilogy(tp, y / 1e6, lw=1.8, color=col, label=lab)
a.set_xlabel("time (ns)"); a.set_ylabel("pressure (MPa, log scale)"); a.legend(frameon=False); a.set_title("Pressure of the vapour trapped between glass and film")
a = ax[1, 1]
if len(VOL):
    a.plot(VOL[:, 0] * 1e9, (VOL[:, 1] / VOL[0, 1] - 1) * 100, color="#0F8F8C", lw=2)
a.set_xlabel("time (ns)"); a.set_ylabel("change of total metal volume (%)"); a.set_title("Metal volume (the small decrease is metal turned to vapour)")
for a in ax.flat:
    a.grid(alpha=0.25); a.spines[["top", "right"]].set_visible(False)
fig.suptitle("%s: main results" % CASE, fontsize=13)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "results_timeseries.png"), dpi=150)

# ---------- table ----------
k = int(np.argmax(A["v"])); vpk, tpk = A["v"][k], A["t"][k]
We = RHO * vpk ** 2 * H / SIGMA; KE = 0.5 * RHO * H * vpk ** 2


def at(tq):
    if tq > A["t"][-1] + 0.1:
        return "—"
    if tq > A["t"][-1]:
        return "%.2f µm (at %.2f ns)" % (-A["disp"][-1] / 1000, A["t"][-1])
    return "%.2f µm" % (-np.interp(tq, A["t"], A["disp"]) / 1000)


det = os.path.join(CASEDIR, "droplet_detections.txt")
dets = open(det).read().strip().splitlines() if os.path.exists(det) else []
sf = os.path.join(CASEDIR, "run_status.txt"); status = open(sf).read() if os.path.exists(sf) else ""
finished = ("stage3 exit 0" in status) or ("phase2 exit 0" in status)
vol_change = (VOL[-1, 1] / VOL[0, 1] - 1) * 100 if len(VOL) else float("nan")
rows = [
    ("Simulated time", "%.2f ns%s" % (A["t"][-1], " (finished)" if finished else " (running)")),
    ("Highest speed of the film on the beam axis", "%.1f m/s at %.2f ns" % (vpk, tpk)),
    ("Speed of the metal on the beam axis at the end (after about 10 ns the gas behind the film pushes the detached piece, so this late value is too high; see docs/RESULTS.md)", "%.1f m/s" % A["v"][-1]),
    ("Weber number at the highest speed (film inertia divided by surface tension, ρv²h/σ)", "%.1f  (= %d · %.1f² · 71.4e-9 / 1.64)" % (We, RHO, vpk)),
    ("Kinetic energy per area divided by the energy of the two new surfaces (½ρhv² / 2σ)", "%.2f / %.2f J/m² = %.2f" % (KE, 2 * SIGMA, KE / (2 * SIGMA))),
    ("Distance moved by the film on the beam axis at 1 / 5 / 20 ns", "%s / %s / %s" % (at(1.0), at(5.0), at(20.0))),
    ("Separate metal pieces crossed by the beam axis at the end (the axis detector only sees a gas gap along the axis; a piece that separates by tearing off the axis is counted in docs/RESULTS.md)", "%d" % int(A["nseg"][-1])),
    ("Detections by the axis droplet detector", "%d%s" % (len(dets), (": first " + dets[0]) if dets else "")),
    ("Change of total metal volume (includes metal turned to vapour)", "%.4f %%" % vol_change),
]
if len(PK):
    rows.append(("Pressure of the vapour trapped on the glass side: peak / at the end", "%.1f MPa / %.3f MPa" % (PK[:, 5].max() / 1e6, PK[-1, 5] / 1e6)))
if G is not None:
    kg = int(np.argmax(G["v"]))
    rows.append(("Comparison run %s: peak axis speed" % REF, "%.1f m/s at %.2f ns" % (G["v"][kg], G["t"][kg])))
with open(os.path.join(OUT, "results_table.md"), "w") as f:
    f.write("# %s: main results\n\n| Quantity | Value |\n|---|---|\n" % CASE)
    for r in rows:
        f.write("| %s | %s |\n" % r)
    f.write("\nThe volumes and speeds of the separated pieces come from tools/count_metal_pieces.py; see docs/RESULTS.md.\n")
print(open(os.path.join(OUT, "results_table.md")).read())
