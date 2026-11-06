#!/usr/bin/env python3
"""Figures of the model inputs of the reference case, for docs/MODEL.md.

    python3 tools/plot_model_inputs.py

Every curve is computed from the input values of cases/fsLIFT_332mJ_axisymmetric
(repeated below) with the formulas the solver implements. Output: figures/model/*.png
"""
import math, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures", "model"); os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3})

# ---------------------------------------------------------------- inputs of the reference case
E_FILM = 4.1089e-8          # pulse energy reaching the film [J]
FWHM = 214e-15              # pulse duration, full width at half maximum [s]
W0 = 3e-6                   # 1/e^2 beam radius [m]
R_REFL = 0.42832            # reflectivity, glass-titanium interface
A_ABS = 1.03469e8           # absorption coefficient [1/m]
H = 71.4e-9                 # film thickness [m]
NFILM = 16                  # cells across the film
GAMMA, CL, KL, KE0 = 630.0, 2.5e6, 4.5, 18.9
G_TAB = np.array([(300, 1e18), (1000, 3e18), (2000, 6e18), (3000, 1e19), (5000, 2e19),
                  (7000, 3e19), (10000, 5e19), (15000, 7e19), (20000, 1e20)])
HV, RV, TV, PREF, BETA, TCAP = 9.1e6, 174.0, 3560.0, 101325.0, 0.18, 10000.0
RHO, SIGMA, MU = 4515.0, 1.64, 2.35e-3
TMELT, TCRIT = 1941.0, 9500.0

sig_t = FWHM / (2 * math.sqrt(2 * math.log(2)))
I0 = E_FILM / ((math.pi * W0**2 / 2) * FWHM * math.sqrt(math.pi / (4 * math.log(2))))   # peak intensity on the film
Q0 = (1 - R_REFL) * I0 * A_ABS                                                         # peak absorbed power density

# ---------------------------------------------------------------- 1. laser
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
t = np.linspace(-600e-15, 600e-15, 400)
ax[0].plot(t * 1e15, np.exp(-t**2 / (2 * sig_t**2)), color="#C62828")
ax[0].axhline(0.5, color="k", lw=0.6, ls=":"); ax[0].annotate("FWHM 214 fs", (110, 0.52), fontsize=8)
ax[0].set_xlabel("time from the pulse peak [fs]"); ax[0].set_ylabel("I(t) / I₀")
ax[0].set_title("(a) Pulse in time: Gaussian, FWHM 214 fs")
d = np.linspace(0, H, 400); faces = np.linspace(0, H, NFILM + 1)
ax[1].semilogy(d * 1e9, Q0 * np.exp(-A_ABS * d), color="#6D28D9", label="exact  (1−R) I₀ a e^(−a d)")
cell = Q0 * (np.exp(-A_ABS * faces[:-1]) - np.exp(-A_ABS * faces[1:])) / (A_ABS * np.diff(faces))
ax[1].step(faces * 1e9, np.r_[cell, cell[-1]], where="post", color="k", lw=1, label="cell averages deposited (16 cells)")
ax[1].set_xlabel("depth below the glass-side surface [nm]"); ax[1].set_ylabel("absorbed power at the pulse peak [W/m³]")
ax[1].set_title("(b) Absorption across the 71.4 nm film"); ax[1].legend(fontsize=7)
r = np.linspace(0, 8e-6, 300)
ax[2].plot(r * 1e6, np.exp(-2 * r**2 / W0**2), color="#0F8F8C")
ax[2].axvline(3, color="k", lw=0.6, ls=":"); ax[2].annotate("w₀ = 3 µm (1/e²)", (3.1, 0.4), fontsize=8)
ax[2].set_xlabel("distance from the beam axis r [µm]"); ax[2].set_ylabel("I(r) / I(0)")
ax[2].set_title("(c) Beam profile in radius")
fig.suptitle("Laser source: peak intensity on the film I₀ = %.3g W/m², peak absorbed power %.3g W/m³" % (I0, Q0), fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "laser_source.png"), dpi=150); plt.close(fig)

# ---------------------------------------------------------------- 2. two-temperature model
Te = np.linspace(300, 20000, 600)
G = np.interp(Te, G_TAB[:, 0], G_TAB[:, 1]); Ce = GAMMA * Te
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
ax[0].plot(Te, Ce / 1e6, label="electrons  C_e = 630 T_e", color="#C62828")
ax[0].axhline(CL / 1e6, color="#1565C0", label="lattice  C_l = 2.5e6 (constant)")
ax[0].set_xlabel("electron temperature T_e [K]"); ax[0].set_ylabel("heat capacity [MJ/(m³ K)]")
ax[0].set_title("(a) Heat capacities"); ax[0].legend(fontsize=7)
ax[1].semilogy(Te, G, color="k"); ax[1].plot(G_TAB[:, 0], G_TAB[:, 1], "o", color="k", ms=4, label="table points (linear interpolation)")
ax[1].set_xlabel("electron temperature T_e [K]"); ax[1].set_ylabel("G [W/(m³ K)]")
ax[1].set_title("(b) Electron–lattice coupling G(T_e)"); ax[1].legend(fontsize=7)
tau_eq = 1 / (G * (1 / Ce + 1 / CL))
ax[2].semilogy(Te, Ce / G * 1e12, label="electron cooling  C_e / G", color="#C62828")
ax[2].semilogy(Te, CL / G * 1e12, label="lattice heating  C_l / G", color="#1565C0")
ax[2].semilogy(Te, tau_eq * 1e12, label="equilibration  1 / [G (1/C_e + 1/C_l)]", color="k", ls="--")
ax[2].set_xlabel("electron temperature T_e [K]"); ax[2].set_ylabel("time [ps]")
ax[2].set_title("(c) Characteristic coupling times"); ax[2].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "two_temperature_properties.png"), dpi=150); plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 3.8))
for Tl, c in ((300, "#7C8798"), (3000, "#0F8F8C"), (10000, "#6D28D9")):
    ax.plot(Te, KE0 * Te / Tl, color=c, label="Drude  k_e = 18.9 T_e/T_l,  T_l = %d K" % Tl)
ax.plot(Te, GAMMA * Te * 1e-4, color="#C62828", ls="--", label="older form  k_e = C_e D_e (D_e = 1e-4 m²/s)")
ax.plot([300, 20000], [KE0, KE0], color="k", lw=0.8, ls=":", label="Drude at equilibrium (T_e = T_l): 18.9")
ax.set_yscale("log"); ax.set_xlabel("electron temperature T_e [K]"); ax.set_ylabel("electron conductivity k_e [W/(m K)]")
ax.set_title("Electron heat conductivity"); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "electron_conductivity.png"), dpi=150); plt.close(fig)

# ---------------------------------------------------------------- 3. phase change
T = np.linspace(2000, 16000, 800); Tc = np.minimum(T, TCAP)
psat = PREF * np.exp(HV / RV * (1 / TV - 1 / Tc))
j = BETA * psat / np.sqrt(2 * math.pi * RV * Tc)
fig, ax = plt.subplots(1, 2, figsize=(12, 4))
ax[0].semilogy(T, psat / 1e6, color="k", label="saturation pressure p_sat")
ax[0].semilogy(T, 0.54 * psat / 1e6, color="#C62828", label="recoil on a free surface 0.54 p_sat")
ax[0].axhline(524, color="#6D28D9", ls="--", lw=1, label="peak trapped-vapour pressure in the run (524 MPa)")
ax[1].semilogy(T, j, color="k", label="evaporation flux j")
ax2 = ax[1].twinx(); ax2.semilogy(T, j / RHO * 1e9 / 1e9, alpha=0)   # same curve, recession speed axis
ax2.set_yscale("log"); lo, hi = ax[1].get_ylim(); ax2.set_ylim(lo / RHO, hi / RHO); ax2.set_ylabel("surface recession speed j/ρ [m/s]  (1 m/s = 1 nm/ns)")
ax2.grid(False)
for a in ax:
    for x, lab, ha in ((TMELT, "melting 1941 K", "right"), (TV, "boiling 3560 K", "right"), (TCRIT, "critical (case value 9500 K)", "right"), (TCAP, "T capped at 10,000 K", "left")):
        a.axvline(x, color="grey", lw=0.6, ls=":")
        a.text(x + (60 if ha == "left" else -60), 0.78 if x < 5000 else 0.45, lab, rotation=90, fontsize=6.5, va="center", ha=ha, color="dimgrey",
               transform=a.get_xaxis_transform())
    a.set_xlabel("surface temperature [K]")
ax[0].set_ylabel("pressure [MPa]"); ax[0].set_title("(a) Clausius–Clapeyron saturation pressure and recoil"); ax[0].legend(fontsize=7, loc="lower right", framealpha=1)
ax[1].set_ylabel("mass flux [kg/(m² s)]"); ax[1].set_title("(b) Hertz–Knudsen evaporation flux, β = 0.18"); ax[1].legend(fontsize=7, loc="lower right")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "phase_change.png"), dpi=150); plt.close(fig)

# ---------------------------------------------------------------- 4. trapped-vapour model (schematic)
fig, ax = plt.subplots(figsize=(10, 4.2)); ax.set_xlim(-0.3, 10.3); ax.set_ylim(-0.6, 4.4); ax.axis("off")
ax.add_patch(Rectangle((0, 3.0), 10, 1.0, color="#B3D7E8")); ax.text(5, 3.5, "donor glass", ha="center", va="center")
prof = lambda x: 2.35 - 0.95 * np.exp(-(x / 3.3) ** 2)                 # dome shape of the film (illustrative)
xs = np.linspace(0, 10, 300)
ax.fill_between(xs, prof(xs) - 0.25, prof(xs), color="#F2C14E"); ax.text(8.6, 2.05, "titanium film", fontsize=8)
for k in range(8):                                                    # rings = vapour columns
    x0, x1 = k * 0.9, (k + 1) * 0.9; xm = 0.5 * (x0 + x1)
    ax.add_patch(Rectangle((x0, prof(xm)), 0.9, 3.0 - prof(xm), facecolor="#E8E0F5", edgecolor="#6D28D9", lw=0.8))
    ax.text(xm, 2.85, "b=%d" % k, ha="center", fontsize=7, color="#6D28D9")
    if k < 4:
        ax.add_patch(FancyArrowPatch((xm, 2.75), (xm, prof(xm) + 0.05), arrowstyle="-|>", mutation_scale=9, color="#6D28D9"))
ax.annotate("", (0, -0.25), (10, -0.25), arrowprops=dict(arrowstyle="<-", color="k")); ax.text(10.1, -0.3, "r", fontsize=10)
ax.plot([0, 0], [-0.4, 4.2], "k-.", lw=0.8); ax.text(0.05, 4.25, "beam axis", fontsize=8)
ax.text(5.0, 1.0, "gap (gas): recoil 0.54 p_sat on the free surface", ha="center", fontsize=8)
ax.text(4.3, 0.35,
        "each ring b (width 100 nm) keeps its own vapour mass m_b, film displacement d_b and pressure p_b;\n"
        "p_b = min(p_sat,b , m_b R T_v,b / V_b) acts on the glass-side surface of the film below the ring",
        ha="center", fontsize=8)
ax.set_title("Trapped-vapour model: one vapour column per radial ring (schematic, not to scale)", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "trapped_vapour_schematic.png"), dpi=150); plt.close(fig)

# ---------------------------------------------------------------- 5. mesh along the beam axis
def section(L, n, er):
    """OpenFOAM simpleGrading section: length L, n cells, expansion ratio er = last/first cell."""
    if abs(er - 1) < 1e-12:
        return np.full(n, L / n)
    q = er ** (1.0 / (n - 1)); first = L * (1 - q) / (1 - q ** n)
    return first * q ** np.arange(n)
def block(L, spec):
    fl = sum(s[0] for s in spec)
    return np.concatenate([section(L * s[0] / fl, s[1], s[2]) for s in spec])
dy = np.concatenate([block(8.0, [(1, 24, 0.05)]),
                     block(12.0, [(8.5, 170, 1), (0.5, 24, 0.12), (3.0, 500, 1)]),
                     block(0.0714, [(1, 16, 1)]),
                     block(12.0, [(0.5, 26, 10), (11.5, 68, 8)])])
yc = np.cumsum(dy) - dy / 2
fig, ax = plt.subplots(figsize=(9, 3.6))
ax.semilogy(yc, dy * 1e3, color="k", lw=1.2)
for y0, y1, c, lab in ((0, 8, "#D9D9D9", "receiver"), (8, 20, "#FFFFFF", "gap"), (20, 20.0714, "#F2C14E", "film"), (20.0714, 32.0714, "#B3D7E8", "glass (gas)")):
    ax.axvspan(y0, y1, color=c, alpha=0.5); ax.text(0.5 * (y0 + y1) if lab != "film" else 20.3, 600, lab, ha="center" if lab != "film" else "left", fontsize=8)
ax.set_xlabel("height y along the beam axis [µm]"); ax.set_ylabel("cell height [nm]")
ax.set_title("Cell height along the beam axis (%d cells; radially 160 cells of 50 nm)" % len(dy), fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "mesh_cell_height.png"), dpi=150); plt.close(fig)

print("I0 = %.4g W/m2, peak absorbed power = %.4g W/m3, cells along y = %d" % (I0, Q0, len(dy)))
print("film cell %.3f nm; smallest gap cell %.3f nm; first glass cell %.3f nm" % (dy[24 + 694] * 1e3, dy[24 + 194:24 + 694].min() * 1e3, dy[24 + 694 + 16] * 1e3))
for TT in (3560, 5000, 7000, 10000):
    p = PREF * math.exp(HV / RV * (1 / TV - 1 / TT)); jj = BETA * p / math.sqrt(2 * math.pi * RV * TT)
    print("T = %5d K: p_sat = %.4g Pa, j = %.4g kg/m2s, recession %.4g m/s" % (TT, p, jj, jj / RHO))
for TT in (1000, 5000, 10000, 16000):
    g = np.interp(TT, G_TAB[:, 0], G_TAB[:, 1]); ce = GAMMA * TT
    print("Te = %5d K: Ce = %.3g, G = %.3g, Ce/G = %.3g ps, Cl/G = %.3g ps, tau_eq = %.3g ps" % (TT, ce, g, ce / g * 1e12, CL / g * 1e12, 1e12 / (g * (1 / ce + 1 / CL))))
print("figures in", os.path.abspath(OUT))
