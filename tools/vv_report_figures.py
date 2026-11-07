#!/usr/bin/env python3
"""Figures for docs/VALIDATION.md.

    python3 tools/vv_report_figures.py

Every number below is either an input of cases/fsLIFT_332mJ_axisymmetric, a result stated in
docs/VERIFICATION.md, or a published value cited in docs/VALIDATION.md (source given next to it).
No figure from a publication is reproduced. Output: figures/validation/*.png
"""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures", "validation"); os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 9})
GREEN, AMBER, RED, GREY = "#2E7D32", "#F9A825", "#C62828", "#9E9E9E"

# ------------------------------------------------------------------ 1. status of the five levels
levels = [("1  Code verification", "15 checks: 12 within 5 %, force balance 6 %, drop frequency 10.5 % and converging, damping open", GREEN),
          ("2  Solution verification", "no mesh-refinement study of the reference case; separation mesh-limited", AMBER),
          ("3  Input validation", "3 titanium inputs differ from measured or first-principles values; density gives the right film mass", RED),
          ("4  Sub-model validation", "mechanisms consistent with titanium and thin-film literature; no matched numbers", AMBER),
          ("5  System validation", "no public measurement under matching conditions", RED)]
fig, ax = plt.subplots(figsize=(10, 4.6)); ax.set_xlim(0, 10); ax.set_ylim(0, 5.2); ax.axis("off")
for i, (name, note, col) in enumerate(levels):
    y0, y1 = i, i + 0.92
    hw0, hw1 = 2.6 - 0.42 * i, 2.6 - 0.42 * (i + 0.92)
    ax.add_patch(Polygon([[3 - hw0, y0], [3 + hw0, y0], [3 + hw1, y1], [3 - hw1, y1]], color=col, alpha=0.85))
    ax.text(3, y0 + 0.46, name, ha="center", va="center", fontsize=8.5, color="white", weight="bold")
    ax.text(6.0, y0 + 0.46, note, ha="left", va="center", fontsize=8.5)
ax.text(0.1, 5.05, "green: done   amber: partial   red: not done or failed", fontsize=8, color="dimgrey")
ax.set_title("Verification and validation status of compInterFoam-fsLIFT", fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "vv_status.png"), dpi=150); plt.close(fig)

# ------------------------------------------------------------------ 2. verification errors (docs/VERIFICATION.md)
tests = [  # (label, |error| in %, run with the released code?)
    ("Laser energy deposited", 0.2, False),
    ("Laser radial profile (largest deviation)", 0.02, False),
    ("Two-temperature: peak electron temperature", 0.5, True),
    ("Two-temperature: electron-lattice meeting time", 0.8, True),
    ("Hertz-Knudsen flux, cell by cell", 1e-4, False),
    ("Recoil pressure 0.54 p_sat", 0.1, False),
    ("Force on the film (worst of 3 intervals)", 6.0, True),
    ("Static drop: Laplace jump, R/dx = 10", 2.75, True),
    ("Static drop: Laplace jump, R/dx = 40", 4.09, True),
    ("Viscous decay rate", 0.12, True),
    ("Galilean test: pressure", 0.1, True),
    ("Metal volume, static drop", 1e-4, True),
    ("Serial vs parallel", 1e-6, False),
    ("Drop oscillation frequency, R/dx = 20", 10.5, True),
    ("Drop oscillation damping, R/dx = 20 (open point)", 70.0, True),
]
fig, ax = plt.subplots(figsize=(9.5, 5.6))
y = np.arange(len(tests))[::-1]
for yi, (lab, e, rel) in zip(y, tests):
    col = RED if e > 20 else AMBER if e > 5 else GREEN
    ax.barh(yi, e, color=col, alpha=0.9 if rel else 0.45, edgecolor="k" if rel else "none", lw=0.6)
    ax.text(e * 1.25, yi, ("%.2g %%" % e) + ("" if rel else "  (development run)"), va="center", fontsize=7.5)
ax.set_yticks(y); ax.set_yticklabels([t[0] for t in tests], fontsize=8)
ax.set_xscale("log"); ax.set_xlim(1e-7, 1e3); ax.set_xlabel("magnitude of the error against the exact or reference value [%]")
ax.axvline(5, color="grey", ls=":", lw=0.8)
ax.set_title("Verification tests: error against exact solutions and an independent code\n(solid bars: run with the released code; pale bars: development runs, unaffected code paths)", fontsize=9.5)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "verification_errors.png"), dpi=150); plt.close(fig)

# ------------------------------------------------------------------ 3. inputs against the literature
fig, ax = plt.subplots(2, 2, figsize=(11.5, 7.6))
T = np.linspace(1700, 4000, 200)
a = ax[0, 0]
a.plot(T, np.full_like(T, 4515), "k", lw=2, label="model: 4515 kg/m³ (constant)")
Tm = np.linspace(1640, 2090, 50); a.plot(Tm, 4193 - 0.23762 * (Tm - 1941), color=GREEN, lw=2.5, label="pure Ti, measured 1640-2090 K (Ozawa et al. 2017)")
a.plot(T, 4193 - 0.23762 * (T - 1941), color=GREEN, ls="--", lw=1, label="same fit extrapolated")
Ta = np.linspace(1661, 1997, 50); a.plot(Ta, 4123 - 0.254 * (Ta - 1943), color="#6D28D9", lw=2, label="Ti-6Al-4V, measured (Li et al. 2006)")
a.axvline(1941, color="grey", ls=":", lw=0.8); a.text(1960, 4480, "melting point", fontsize=7, color="grey")
a.set_xlabel("temperature [K]"); a.set_ylabel("density [kg/m³]"); a.set_title("(a) Density of liquid titanium"); a.legend(fontsize=7)
a = ax[0, 1]
a.plot(T, np.full_like(T, 1.64), "k", lw=2, label="model: 1.64 N/m (constant)")
a.plot(T, 1.613 - 0.2049e-3 * (T - 1941), color=GREEN, lw=2.5, label="pure Ti, levitation (Seimiya et al. 2022)")
a.axvline(1941, color="grey", ls=":", lw=0.8)
a.set_xlabel("temperature [K]"); a.set_ylabel("surface tension [N/m]"); a.set_title("(b) Surface tension of liquid titanium"); a.legend(fontsize=7)
Te = np.linspace(300, 20000, 300)
GT = np.array([(300, 1e18), (1000, 3e18), (2000, 6e18), (3000, 1e19), (5000, 2e19), (7000, 3e19), (10000, 5e19), (15000, 7e19), (20000, 1e20)])
a = ax[1, 0]
a.semilogy(Te, np.interp(Te, GT[:, 0], GT[:, 1]), "k", lw=2, label="model table")
a.axhspan(2.5e18, 3.7e18, color=GREEN, alpha=0.3, label="published peak for Ti: 2.5e18 (Hobmaier et al. 2026)\nto about 3.7e18 (Lin et al. 2008, as reported by Hobmaier)")
a.set_xlabel("electron temperature [K]"); a.set_ylabel("G [W/(m³ K)]"); a.set_title("(c) Electron-lattice coupling G(T_e): model against published peak")
a.legend(fontsize=7, loc="lower right")
a = ax[1, 1]
a.plot(Te, 630 * Te / 1e6, "k", lw=2, label="model: 630 T_e")
a.plot(Te, 318 * Te / 1e6, color=GREEN, lw=1.5, ls="--", label="Sommerfeld fit 318 T_e (Hobmaier et al. 2026)")
a.plot([12000], [2.29], "o", color=GREEN, ms=8, label="first-principles peak 2.29 MJ/(m³ K) near 12,000 K")
a.set_xlabel("electron temperature [K]"); a.set_ylabel("C_e [MJ/(m³ K)]"); a.set_title("(d) Electron heat capacity C_e(T_e)"); a.legend(fontsize=7)
fig.suptitle("Model inputs of the reference case against measured or first-principles values for titanium", fontsize=10.5)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "input_validation.png"), dpi=150); plt.close(fig)

# ------------------------------------------------------------------ 4. first-order effect of the measured rho and sigma
# The mass per area of the film is set by the deposited film: 71.4 nm of solid Ti at 4506 kg/m3. The
# liquid keeps that mass and expands, h_liq = h * 4506 / rho_liq(T). The model keeps rho = 4515 and h.
h, RHO_S = 71.4e-9, 4506.0
def lit(Tk): return 4193 - 0.23762 * (Tk - 1941), 1.613 - 0.2049e-3 * (Tk - 1941)
rows = []
for Tk in (2000, 2500, 3000):
    r, s = lit(Tk); r0, s0 = 4515.0, 1.64; hl = h * RHO_S / r
    rows.append((Tk, [100 * (r0 / RHO_S - 1),                                   # acceleration p/(rho h)
                      100 * ((r * hl / s) / (r0 * h / s0) - 1),                 # Weber number rho v^2 h / sigma
                      100 * (np.sqrt(s / (r * hl)) / np.sqrt(s0 / (r0 * h)) - 1),   # Taylor-Culick speed
                      100 * (np.sqrt(r * hl ** 3 / s) / np.sqrt(r0 * h ** 3 / s0) - 1)]))  # capillary time
labels = ["launch acceleration p/(ρh)\n(same pressure)", "Weber number ρv²h/σ\n(same speed)", "sheet retraction speed\n√(2σ/(ρh))", "capillary time\n√(ρh³/σ)"]
fig, ax = plt.subplots(figsize=(9, 4.2)); x = np.arange(len(labels)); wdt = 0.26
for k, (Tk, v) in enumerate(rows):
    b = ax.bar(x + (k - 1) * wdt, v, wdt, label="metal at %d K" % Tk, color=["#90CAF9", "#42A5F5", "#1565C0"][k])
    for xi, vi in zip(x + (k - 1) * wdt, v): ax.text(xi, vi + (0.4 if vi >= 0 else -1.6), "%+.0f%%" % vi, ha="center", fontsize=7)
ax.axhline(0, color="k", lw=0.7); ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("change with the measured ρ(T), σ(T) [%]")
ax.set_title("First-order effect of the measured liquid density and surface tension of titanium, film mass fixed\n(the liquid film expands to keep the deposited mass; scaling estimates, not a rerun; density fit extrapolated above 2,090 K)", fontsize=9)
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(os.path.join(OUT, "input_effect_estimate.png"), dpi=150); plt.close(fig)
for Tk, v in rows: print(Tk, ["%+.1f" % x for x in v], "rho=%.0f sigma=%.3f" % lit(Tk))
print("figures in", os.path.abspath(OUT))
