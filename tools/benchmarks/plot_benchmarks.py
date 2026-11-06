#!/usr/bin/env python3
"""Figures for the drop benchmarks.

  static drop      pressure jump against the Laplace law 2*sigma/R, and spurious currents
  oscillating drop deformation of an l = 2 mode against Lamb's damped oscillation

Usage:  python3 tools/benchmarks/plot_benchmarks.py <benchmark run directory>
The run directory holds the cases built by make_drop.sh, each with the drop_metrics.csv
written by drop_analyze.py: sd_dx20, sd_dx10, sd_dx5 (static drop, cell size in nm) and
os_dx20, os_dx10 (oscillating drop). Missing cases are skipped.
Output: figures/verification/static_drop.png and oscillating_drop.png."""
import os, sys, math, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

RUNS = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else '.')
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', 'figures', 'verification'); os.makedirs(OUT, exist_ok=True)

# liquid titanium drop, R = 0.2 um
sig, rho, mu, R = 1.64, 4515.0, 2.35e-3, 0.2e-6
nu = mu / rho; dp0 = 2 * sig / R
w0 = math.sqrt(8 * sig / (rho * R**3)); beta0 = 5 * nu / R**2; wd = math.sqrt(w0**2 - beta0**2)

def load(c):
    p = os.path.join(RUNS, c, 'drop_metrics.csv')
    return np.atleast_2d(np.loadtxt(p, delimiter=',')) if os.path.exists(p) else None

def at(D, t):
    return D[np.argmin(abs(D[:, 0] - t))]

# static drop: error of the pressure jump and largest spurious velocity at t = 0.2 ns
xs, es, us = [], [], []
for dx in (20, 10, 5):
    D = load('sd_dx%d' % dx)
    if D is None: continue
    r = at(D, 2e-10); xs.append(R / (dx * 1e-9)); es.append(100 * (r[2] / dp0 - 1)); us.append(r[3])
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
ax[0].plot(xs, es, 'o-', color='tab:green')
for x, e in zip(xs, es): ax[0].annotate('%+.2f %%' % e, (x, e), textcoords='offset points', xytext=(6, 6), fontsize=9)
ax[0].axhline(0, color='k', lw=0.7); ax[0].set_xscale('log', base=2)
ax[0].set_xlabel('drop radius / cell size  R/Δx'); ax[0].set_ylabel('error of the pressure jump [%]')
ax[0].set_title('Static drop: pressure jump against the Laplace law 2σ/R = %.1f MPa' % (dp0 / 1e6), fontsize=10)
ax[1].plot(xs, us, 's-', color='tab:blue'); ax[1].set_xscale('log', base=2)
ax[1].set_xlabel('drop radius / cell size  R/Δx'); ax[1].set_ylabel('largest velocity in the field [m/s]')
ax[1].set_title('Static drop: spurious velocity (exact value 0)', fontsize=10)
for a_ in ax:
    a_.set_xticks(xs); a_.set_xticklabels(['%d' % x for x in xs]); a_.minorticks_off()
fig.suptitle('Liquid titanium drop, R = 0.2 µm, at 0.2 ns; settings of the reference case', fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'static_drop.png'), dpi=140)

# oscillating drop: deformation (b - a)/R against Lamb
fig, ax = plt.subplots(figsize=(8.5, 4.8))
t = np.linspace(0, 12e-9, 400); first = None
for c, lab, col in (('os_dx20', 'simulation, R/Δx = 10', 'tab:blue'), ('os_dx10', 'simulation, R/Δx = 20', 'navy')):
    D = load(c)
    if D is None: continue
    y = (D[:, 4] - D[:, 5]) / R
    ax.plot(D[:, 0] * 1e9, y, 'o-', ms=3, color=col, label=lab)
    if first is None: first = (D[0, 0], y[0])
if first is not None:
    A = first[1] / (math.exp(-beta0 * first[0]) * math.cos(wd * first[0]))
    ax.plot(t * 1e9, A * np.exp(-beta0 * t) * np.cos(wd * t), 'k--', lw=1.5,
            label='Lamb (l = 2): period %.2f ns, decay time %.1f ns' % (2 * math.pi / wd * 1e9, 1e9 / beta0))
ax.axhline(0, color='k', lw=0.5); ax.set_xlabel('time [ns]'); ax.set_ylabel('deformation (b − a)/R')
ax.set_title('Oscillating titanium drop (R = 0.2 µm): surface tension, inertia and viscosity', fontsize=10)
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(os.path.join(OUT, 'oscillating_drop.png'), dpi=140)
print('figures in', os.path.abspath(OUT))
