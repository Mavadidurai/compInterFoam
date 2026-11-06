#!/usr/bin/env python3
"""
1D finite-difference two-temperature model (TTM) for a thin Ti donor film
irradiated through the glass carrier (fs-LIFT geometry).

    Ce(Te) dTe/dt = d/dz( ke(Te,Tl) dTe/dz ) - G(Te) (Te - Tl) + S(z,t)
    Cl     dTl/dt = d/dz( kl        dTl/dz ) + G(Te) (Te - Tl)

    Ce = gamma*Te,  ke = ke0*Te/Tl (Drude, valid for Te << T_Fermi)
    S  = (1-R) F0 * alpha*exp(-alpha z) * f(t),  f = Gaussian, integral 1

Discretisation:
  * cell-centred grid, N cells across the film, zero-flux at both faces
  * backward Euler in time, Te and Tl solved SIMULTANEOUSLY
    (2x2 block-tridiagonal system, block-Thomas algorithm)
  * Picard iterations on the nonlinear coefficients, with
    Ce evaluated as gamma*(Te^n + Te^{n+1})/2 so the converged step conserves
    the electron energy gamma*Te^2/2 exactly
  * the source is averaged EXACTLY over each cell (Beer-Lambert) and over each
    time step (erf of the Gaussian), so the discrete absorbed energy equals
    the analytical one.

No phase change / hydrodynamics: lattice temperatures above T_melt are
"energy markers" only. Requires numpy (matplotlib optional for plots).
"""
import argparse
import math

import numpy as np

# ---------------------------------------------------------------- parameters
# Values below are the inputs of cases/column1D_laser_TTM.
P = dict(
    L=71.4e-9,          # film thickness [m]
    alpha=1.07e7,       # absorption coefficient [1/m] (as in laserProperties)
    R=0.35,             # reflectivity [-]
    E_pulse=6e-8,       # pulse energy [J]
    w0=3e-6,            # 1/e^2 beam radius [m]
    tau=200e-15,        # FWHM pulse duration [s]
    gamma=630.0,        # Sommerfeld coefficient [J/m^3/K^2] (as in the case)
    Cl=2.5e6,           # lattice heat capacity [J/m^3/K]
    kl=4.5,             # lattice conductivity [W/m/K]
    ke0=17.4,           # electron conductivity at Te=Tl=300 K [W/m/K]
    T0=300.0,           # initial temperature [K]
)
# G(Te) table [K, W/m^3/K], as in the case controlDict
G_TABLE = np.array([
    (300, 1.0e18), (1000, 3.0e18), (2000, 6.0e18), (3000, 1.0e19),
    (5000, 2.0e19), (7000, 3.0e19), (10000, 5.0e19), (15000, 7.0e19),
    (20000, 1.0e20)])


def G_of(Te):
    return np.interp(Te, G_TABLE[:, 0], G_TABLE[:, 1])


def peak_fluence(p):
    """Peak fluence of a Gaussian beam: F0 = 2E/(pi w0^2)."""
    return 2.0 * p["E_pulse"] / (math.pi * p["w0"] ** 2)


def block_thomas(Lo, D, Up, b):
    """Solve block-tridiagonal system with 2x2 blocks.
    Lo[i], Up[i] are diagonal 2-vectors (couple same species only),
    D[i] is a full 2x2 block, b[i] a 2-vector. Returns x (N,2)."""
    N = D.shape[0]
    Dp = D.copy()
    bp = b.copy()
    for i in range(1, N):
        # M = Lo[i] * inv(Dp[i-1])  (Lo diagonal)
        a, bb, c, d = Dp[i - 1, 0, 0], Dp[i - 1, 0, 1], Dp[i - 1, 1, 0], Dp[i - 1, 1, 1]
        det = a * d - bb * c
        inv = np.array([[d, -bb], [-c, a]]) / det
        M = Lo[i][:, None] * inv
        Dp[i] -= M * Up[i - 1][None, :]       # M @ diag(Up[i-1])
        bp[i] -= M @ bp[i - 1]
    x = np.empty_like(b)
    x[-1] = np.linalg.solve(Dp[-1], bp[-1])
    for i in range(N - 2, -1, -1):
        x[i] = np.linalg.solve(Dp[i], bp[i] - Up[i] * x[i + 1])
    return x


def run(N=40, dt=1e-15, t_end=20e-12, p=P, picard_tol=1e-10, max_picard=30,
        ke_model="drude", record_every=1):
    L, alpha, gamma, Cl, kl = p["L"], p["alpha"], p["gamma"], p["Cl"], p["kl"]
    dz = L / N
    zf = np.linspace(0.0, L, N + 1)                   # faces
    # exact cell average of alpha*exp(-alpha z)
    shape = (np.exp(-alpha * zf[:-1]) - np.exp(-alpha * zf[1:])) / dz
    Fa = (1.0 - p["R"]) * peak_fluence(p)             # absorbed-side fluence
    t0 = 3.0 * p["tau"]                               # pulse peak time
    a = 2.0 * math.sqrt(math.log(2.0)) / p["tau"]

    def pulse_fraction(t1, t2):                        # int_t1^t2 f dt
        return 0.5 * (math.erf(a * (t2 - t0)) - math.erf(a * (t1 - t0)))

    def ke_of(Te, Tl):
        if ke_model == "drude":
            return p["ke0"] * Te / Tl
        if ke_model == "repo":                         # ke = Ce*De, De=1e-4
            return gamma * Te * 1e-4
        raise ValueError(ke_model)

    def harm(k):                                        # face values, zero-flux BC
        kf = np.zeros(N + 1)
        kf[1:-1] = 2.0 * k[:-1] * k[1:] / (k[:-1] + k[1:])
        return kf

    Te = np.full(N, p["T0"])
    Tl = np.full(N, p["T0"])
    nsteps = int(round(t_end / dt))
    E0 = dz * np.sum(0.5 * gamma * Te ** 2 + Cl * Tl)
    E_in = 0.0
    E_total = Fa * (1.0 - math.exp(-alpha * L))   # analytic absorbed [J/m^2]
    hist = {"t": [], "Te_front": [], "Te_back": [], "Tl_front": [], "Tl_back": [],
            "Te_max": [], "energy_error": [], "picard": []}
    klf = harm(np.full(N, kl))
    for n in range(nsteps):
        t = n * dt
        S = Fa * shape * pulse_fraction(t, t + dt) / dt   # step-averaged [W/m^3]
        Te_k, Tl_k = Te.copy(), Tl.copy()
        for it in range(max_picard):
            Ce = 0.5 * gamma * (Te + Te_k)
            kef = harm(ke_of(Te_k, Tl_k))
            G = G_of(Te_k)
            D = np.zeros((N, 2, 2))
            D[:, 0, 0] = Ce / dt + (kef[:-1] + kef[1:]) / dz ** 2 + G
            D[:, 0, 1] = -G
            D[:, 1, 0] = -G
            D[:, 1, 1] = Cl / dt + (klf[:-1] + klf[1:]) / dz ** 2 + G
            Lo = np.stack([-kef[:-1], -klf[:-1]], axis=1) / dz ** 2
            Up = np.stack([-kef[1:], -klf[1:]], axis=1) / dz ** 2
            b = np.stack([Ce / dt * Te + S, Cl / dt * Tl], axis=1)
            x = block_thomas(Lo, D, Up, b)
            change = np.max(np.abs(x[:, 0] - Te_k) / x[:, 0])
            Te_k, Tl_k = x[:, 0], x[:, 1]
            if change < picard_tol:
                break
        Te, Tl = Te_k, Tl_k
        E_in += dz * np.sum(S) * dt
        if n % record_every == 0 or n == nsteps - 1:
            E = dz * np.sum(0.5 * gamma * Te ** 2 + Cl * Tl)
            hist["t"].append(t + dt)
            hist["Te_front"].append(Te[0]); hist["Te_back"].append(Te[-1])
            hist["Tl_front"].append(Tl[0]); hist["Tl_back"].append(Tl[-1])
            hist["Te_max"].append(Te.max())
            hist["energy_error"].append((E - E0 - E_in) / E_total)  # rel. to total absorbed
            hist["picard"].append(it + 1)
    hist = {k: np.array(v) for k, v in hist.items()}
    return dict(Te=Te, Tl=Tl, z=0.5 * (zf[:-1] + zf[1:]), hist=hist,
                E_in=E_in, E_in_exact=Fa * (1 - math.exp(-alpha * L)))


def convergence(t_end=3e-12):
    print("\n--- temporal convergence (N=20, Tl_back at t=%.1f ps) ---" % (t_end * 1e12))
    ref = run(N=20, dt=0.125e-15, t_end=t_end)["Tl"][-1]
    prev = None
    for dt in (4e-15, 2e-15, 1e-15, 0.5e-15):
        err = abs(run(N=20, dt=dt, t_end=t_end)["Tl"][-1] - ref)
        order = "" if prev is None else "  order=%.2f" % math.log2(prev / err)
        print("dt=%5.2f fs  |error|=%.4e K%s" % (dt * 1e15, err, order))
        prev = err
    print("\n--- spatial convergence (dt=1 fs, Tl_back at t=%.1f ps) ---" % (t_end * 1e12))
    ref = run(N=320, dt=1e-15, t_end=t_end)["Tl"][-1]
    prev = None
    for N in (10, 20, 40, 80):
        err = abs(run(N=N, dt=1e-15, t_end=t_end)["Tl"][-1] - ref)
        order = "" if prev is None else "  order=%.2f" % math.log2(prev / err)
        print("N=%4d  dz=%5.2f nm  |error|=%.4e K%s" % (N, P["L"] / N * 1e9, err, order))
        prev = err


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=40)
    ap.add_argument("--dt", type=float, default=1e-15)
    ap.add_argument("--t_end", type=float, default=20e-12)
    ap.add_argument("--gamma", type=float, default=P["gamma"])
    ap.add_argument("--convergence", action="store_true")
    ap.add_argument("--plot", default="")
    args = ap.parse_args()
    p = dict(P, gamma=args.gamma)

    F0 = peak_fluence(p)
    print("Peak fluence F0 = 2E/(pi w0^2) = %.1f J/m^2 = %.1f mJ/cm^2" % (F0, F0 / 10))
    print("Absorbed in film (1-R)F0(1-exp(-alpha L)) = %.1f J/m^2  (%.1f%% of (1-R)F0)"
          % ((1 - p["R"]) * F0 * (1 - math.exp(-p["alpha"] * p["L"])),
             100 * (1 - math.exp(-p["alpha"] * p["L"]))))
    r = run(N=args.N, dt=args.dt, t_end=args.t_end, p=p, record_every=10)
    h = r["hist"]
    i_pk = int(np.argmax(h["Te_max"]))
    # equilibration: first time front Te-Tl < 5% of Tl
    eq = np.where(np.abs(h["Te_front"] - h["Tl_front"]) < 0.05 * h["Tl_front"])[0]
    eq = eq[eq > i_pk]
    print("gamma = %.0f J/m^3/K^2, N=%d, dt=%.2f fs" % (p["gamma"], args.N, args.dt * 1e15))
    print("Peak Te          = %.0f K at t = %.3f ps (pulse peak at %.3f ps)"
          % (h["Te_max"][i_pk], h["t"][i_pk] * 1e12, 3 * p["tau"] * 1e12))
    if eq.size:
        print("Te~Tl (5%%) after t = %.2f ps" % (h["t"][eq[0]] * 1e12))
    print("Final Tl front/back = %.0f / %.0f K at t = %.1f ps"
          % (r["Tl"][0], r["Tl"][-1], args.t_end * 1e12))
    print("Energy: discrete absorbed %.6e J/m^2, analytic %.6e J/m^2"
          % (r["E_in"], r["E_in_exact"]))
    print("Max relative energy-balance error = %.2e" % np.max(np.abs(h["energy_error"])))
    print("Max Picard iterations per step     = %d" % h["picard"].max())

    if args.convergence:
        convergence()

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        tp = h["t"] * 1e12
        ax[0].plot(tp, h["Te_front"], label="$T_e$, glass side (z=0)")
        ax[0].plot(tp, h["Tl_front"], label="$T_l$, glass side (z=0)")
        ax[0].plot(tp, h["Te_back"], "--", label="$T_e$, free side (z=L)")
        ax[0].plot(tp, h["Tl_back"], "--", label="$T_l$, free side (z=L)")
        ax[0].axhline(1941, color="grey", lw=0.8, ls=":")
        ax[0].text(tp[-1], 1941, " Ti $T_m$", va="bottom", ha="right", fontsize=8)
        ax[0].set_xscale("log"); ax[0].set_xlim(0.1, tp[-1])
        ax[0].set_xlabel("t [ps]"); ax[0].set_ylabel("T [K]")
        ax[0].set_title("Ti film 71.4 nm, F0=%.0f mJ/cm$^2$, 200 fs (no phase change)"
                        % (F0 / 10), fontsize=9)
        ax[0].legend(fontsize=8)
        ax[1].plot(r["z"] * 1e9, r["Te"], label="$T_e$")
        ax[1].plot(r["z"] * 1e9, r["Tl"], label="$T_l$")
        ax[1].set_xlabel("z [nm] (0 = glass/Ti interface)"); ax[1].set_ylabel("T [K]")
        ax[1].set_title("Profiles at t = %.0f ps" % (args.t_end * 1e12), fontsize=9)
        ax[1].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(args.plot, dpi=150)
        print("saved", args.plot)


if __name__ == "__main__":
    main()
