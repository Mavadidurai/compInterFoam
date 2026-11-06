#!/usr/bin/env python3
"""Compare the 1-D column run (cases/column1D_laser_TTM) with the independent
finite-difference two-temperature code tools/ttm_fd_1d.py, which uses the same inputs.

Usage, after the run:   python3 tools/compare_ttm_1d.py <run directory of column1D_laser_TTM>
Prints the peak electron temperature on the glass side of the film and the time at
which the electron and lattice temperatures meet (within 5 %), for both codes.
Figure: figures/verification/ttm_column_vs_fd.png
"""
import math, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ttm_fd_1d as fd  # noqa: E402

CASE = os.path.abspath(sys.argv[1])
OUT = os.path.join(HERE, "..", "figures", "verification"); os.makedirs(OUT, exist_ok=True)
LABELS = ["film, glass side (laser entry)", "film centre", "film, gap side"]


def read_probe(field):
    t, rows = [], []
    with open(os.path.join(CASE, "postProcessing", "filmProbes", "0", field)) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            t.append(float(parts[0])); rows.append([float(x) for x in parts[1:]])
    return np.array(t), np.array(rows)


def meet_time(t, te, tl):
    """first time after the electron peak at which |Te - Tl| / Tl < 5 %"""
    d = np.abs(te - tl) / tl
    ok = np.flatnonzero((d < 0.05) & (np.arange(len(t)) > np.argmax(te)))
    return t[ok[0]] if ok.size else float("nan")


t, Te = read_probe("Te")
_, Tl = read_probe("Tl")
# the solver centres the pulse at min(0.5*(laserEnd - laserStart), 3 sigma) when pulseCenterTime is absent
sigma = fd.P["tau"] / (2 * math.sqrt(2 * math.log(2)))
tc = min(0.5 * 2e-10, 3 * sigma)
shift = 3.0 * fd.P["tau"] - tc                     # the FD code puts the pulse peak at 3 tau
r = fd.run(N=32, dt=1e-15, t_end=t[-1] + shift, ke_model="repo", record_every=10)  # ke = Ce*De, as in the case
h = r["hist"]; tfd = h["t"] - shift

pk_s, pk_f = Te[:, 0].max(), h["Te_front"].max()
eq_s, eq_f = meet_time(t, Te[:, 0], Tl[:, 0]), meet_time(tfd, h["Te_front"], h["Tl_front"])
print("peak electron temperature, glass side: solver %.0f K, finite-difference code %.0f K (difference %+.1f %%)"
      % (pk_s, pk_f, 100 * (pk_s / pk_f - 1)))
print("electron and lattice temperatures meet (5 %%): solver %.3f ps, finite-difference code %.3f ps (difference %+.1f %%)"
      % (eq_s * 1e12, eq_f * 1e12, 100 * (eq_s / eq_f - 1)))

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8.5, 4.8))
for k, c in zip(range(3), ["C0", "C1", "C2"]):
    ax.plot(t * 1e12, Te[:, k], c, label="electrons, " + LABELS[k])
    ax.plot(t * 1e12, Tl[:, k], c, ls="--", label="lattice, " + LABELS[k])
ax.plot(tfd * 1e12, h["Te_front"], "k:", lw=1.2, label="independent finite-difference code: electrons, glass side")
ax.plot(tfd * 1e12, h["Tl_front"], color="grey", ls=":", lw=1.2, label="independent finite-difference code: lattice, glass side")
ax.set_xscale("log"); ax.set_xlim(0.05, t[-1] * 1e12)
ax.set_xlabel("time [ps]"); ax.set_ylabel("temperature [K]")
ax.set_title("1-D column: solver (solid: electrons, dashed: lattice) against the independent finite-difference code (dotted)", fontsize=9)
ax.legend(fontsize=7); fig.tight_layout()
fig.savefig(os.path.join(OUT, "ttm_column_vs_fd.png"), dpi=140)
print("saved", os.path.abspath(os.path.join(OUT, "ttm_column_vs_fd.png")))
