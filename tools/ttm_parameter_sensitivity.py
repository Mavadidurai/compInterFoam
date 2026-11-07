#!/usr/bin/env python3
"""Sensitivity of the laser heating to the two-temperature parameters of titanium.

Solves the 1-D two-temperature model across the 71.4 nm film (tools/ttm_fd_1d.py) with the
physical laser inputs of cases/fsLIFT_332mJ_axisymmetric (no evaporation, no flow) for three
parameter sets:
  A  the reference case: C_e = 630 T_e, G(T_e) table rising to 1e20 W/(m3 K)
  B  as A, but G capped at 3.7e18 W/(m3 K): the peak of Lin, Zhigilei & Celli (2008) for Ti as
     implied by Hobmaier et al. (2026), who report their own peak (2.5e18) as about 1.5 times lower
  C  first-principles values of Hobmaier et al. (2026) for Ti: gamma = 318 J/(m3 K2), G capped at
     2.5e18 W/(m3 K). C_e is kept linear, so above about 7,200 K it exceeds their computed C_e,
     which peaks at 2.29 MJ/(m3 K) near 12,000 K.
The electron conductivity is the Drude form k_e = 18.9 T_e/T_l in all three.

    python3 tools/ttm_parameter_sensitivity.py            # run the three cases (about 5 min each) and plot
    python3 tools/ttm_parameter_sensitivity.py --plot      # plot from the stored data only
Output: figures/literature/ttm_parameter_sensitivity.png, figures/literature/data/*.csv
"""
import math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures", "literature"); DATA = os.path.join(OUT, "data")
os.makedirs(DATA, exist_ok=True)
CASES = {"A": "A: reference case (G table up to 1e20, C_e = 630 T_e)",
         "B": "B: G capped at 3.7e18 (Lin et al. 2008 level for Ti)",
         "C": "C: Hobmaier et al. 2026 Ti (G ≤ 2.5e18, C_e = 318 T_e)"}


def run_case(case):
    sys.path.insert(0, HERE)
    import ttm_fd_1d as fd
    p = dict(fd.P); p.update(alpha=1.03469e8, R=0.42832, E_pulse=4.1089e-8, tau=214e-15, ke0=18.9)
    table = fd.G_TABLE.copy()
    if case == "B":
        fd.G_TABLE = np.c_[table[:, 0], np.minimum(table[:, 1], 3.7e18)]
    elif case == "C":
        fd.G_TABLE = np.c_[table[:, 0], np.minimum(table[:, 1], 2.5e18)]; p["gamma"] = 318.0
    r = fd.run(N=64, dt=1e-15, t_end=20e-12, p=p, ke_model="drude", record_every=10)
    fd.G_TABLE = table
    h = r["hist"]; t = h["t"] - 3 * p["tau"]          # time after the pulse peak
    np.savetxt(os.path.join(DATA, "ttm_sensitivity_case_%s.csv" % case),
               np.c_[t, h["Te_front"], h["Tl_front"], h["Tl_back"]], delimiter=",",
               header="t_after_pulse_peak_s,Te_glass_side_K,Tl_glass_side_K,Tl_gap_side_K")


def summary(D):
    t, te, tl, tb = D.T
    d = np.abs(te - tl) / tl; k = np.flatnonzero((d < 0.05) & (np.arange(len(t)) > te.argmax()))
    i2 = np.argmin(abs(t - 2e-12))
    return dict(Te_max=te.max(), Tl_max=tl.max(), t_Tl_max=t[tl.argmax()] * 1e12,
                t_meet=t[k[0]] * 1e12 if k.size else float("nan"), Tl_2ps=tl[i2], Tl_20ps=tl[-1], Tback_20ps=tb[-1],
                t_above_10kK=(t[np.flatnonzero(tl > 1e4)[-1]] * 1e12) if (tl > 1e4).any() else 0.0)


if __name__ == "__main__":
    if "--plot" not in sys.argv:
        for c in CASES:
            run_case(c)
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.3))
    for (c, lab), col in zip(CASES.items(), ("k", "#1565C0", "#C62828")):
        D = np.loadtxt(os.path.join(DATA, "ttm_sensitivity_case_%s.csv" % c), delimiter=",")
        t = D[:, 0] * 1e12; s = summary(D)
        ax[0].plot(t, D[:, 1] / 1e3, color=col, ls="--", lw=1)
        ax[0].plot(t, D[:, 2] / 1e3, color=col, lw=1.6, label=lab)
        ax[1].plot(t, D[:, 2] / 1e3, color=col, lw=1.6, label=lab)
        print("%s: Te max %.0f K, Tl max %.0f K at %.2f ps, meet %.2f ps, Tl at 2 ps %.0f K, at 20 ps %.0f K (gap side %.0f K), above 10 kK until %.1f ps"
              % (c, s["Te_max"], s["Tl_max"], s["t_Tl_max"], s["t_meet"], s["Tl_2ps"], s["Tl_20ps"], s["Tback_20ps"], s["t_above_10kK"]))
    ax[0].set_xlim(-0.5, 5); ax[0].set_xlabel("time after the pulse peak [ps]"); ax[0].set_ylabel("temperature at the glass-side surface [kK]")
    ax[0].set_title("First 5 ps (dashed: electrons, solid: lattice)", fontsize=10); ax[0].legend(fontsize=7.5)
    ax[1].axhline(10, color="grey", ls=":", lw=1); ax[1].text(0.3, 10.2, "temperature cap of p_sat and evaporation (10 kK)", fontsize=7, color="grey")
    ax[1].set_xlabel("time after the pulse peak [ps]"); ax[1].set_ylabel("lattice temperature at the glass-side surface [kK]")
    ax[1].set_title("Lattice temperature to 20 ps", fontsize=10)
    fig.suptitle("1-D two-temperature model across the 71.4 nm film, physical laser inputs (no evaporation, no flow)", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "ttm_parameter_sensitivity.png"), dpi=150)
    print("figure in", os.path.abspath(OUT))
