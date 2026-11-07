#!/usr/bin/env python3
"""Figure: speed of the detached droplet as a body against the speed sampled on the beam axis.

    python3 tools/plot_droplet_motion.py [case name, default fsLIFT_332mJ_axisymmetric]

Reads figures/<case>/data/droplet_motion.csv (tools/droplet_motion.py) and
figures/<case>/data/axis_timeseries.csv (tools/axis_timeseries.py).
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
case = sys.argv[1] if len(sys.argv) > 1 else "fsLIFT_332mJ_axisymmetric"
D = os.path.join(HERE, "..", "figures", case, "data")
M = np.loadtxt(os.path.join(D, "droplet_motion.csv"), delimiter=",", skiprows=1)
A = np.loadtxt(os.path.join(D, "axis_timeseries.csv"), delimiter=",")
w = (A[:, 0] > 9e-9)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
ax[0].plot(A[w, 0] * 1e9, -A[w, 3], color="#9E9E9E", lw=1.5, label="local liquid velocity sampled on the beam axis")
ax[0].plot(M[:, 0], -M[:, 3], "o-", color="#C62828", lw=2, label="detached droplet: centre-of-mass velocity")
ax[0].plot(M[:, 0], -M[:, 7], "s--", color="#1565C0", lw=1.2, ms=4, label="all metal below the film: mean velocity")
ax[0].axvline(13.4, color="k", ls=":", lw=0.8); ax[0].text(13.5, 40, "walls torn", fontsize=8)
ax[0].set_xlabel("time [ns]"); ax[0].set_ylabel("speed towards the receiver [m/s]"); ax[0].set_ylim(0, 160)
ax[0].set_title("The droplet flies at a constant ~100 m/s;\nthe axis sample rises because the droplet changes shape", fontsize=9.5); ax[0].legend(fontsize=7.5, loc="lower left")
ax[1].plot(M[:, 0], M[:, 2], "o-", color="#C62828", label="droplet centre of mass")
ax[1].plot(M[:, 0], M[:, 4], "^-", color="#6D28D9", label="droplet's lowest point on the axis")
fit = np.polyfit(M[:, 0] * 1e-9, M[:, 2] * 1e-6, 1)
ax[1].set_xlabel("time [ns]"); ax[1].set_ylabel("height y [µm]")
ax[1].set_title("Position: centre of mass moves at %.0f m/s (linear fit)" % (-fit[0]), fontsize=9.5); ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "..", "figures", case, "droplet_motion.png"), dpi=150)
print("centre-of-mass speed from a linear fit of position: %.1f m/s; mass (full circle) %.3e kg = %.3f um3" % (-fit[0], M[-1, 1] * 72, M[-1, 1] * 72 / 4515 * 1e18))
