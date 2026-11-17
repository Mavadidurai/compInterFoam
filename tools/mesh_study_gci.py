#!/usr/bin/env python3
"""Grid-convergence analysis of the mesh study (tools/run_mesh_study.sh).

    python3 tools/mesh_study_gci.py <work dir> [time in s, default 8.02e-10]

Needs <work dir>/L1, L2, L3, each with its logs and mesh_study_fields.csv
(xvfb-run -a pvpython tools/mesh_study_fields.py <work dir>/L<n>).

For each quantity it applies the procedure of Celik et al. (2008, J. Fluids Eng. 130, 078001),
which extends Roache's grid convergence index to non-constant refinement ratios:
  h_i = (1/N_i)^(1/2) for a 2-D mesh of N_i cells (1 = finest), r21 = h2/h1, r32 = h3/h2,
  apparent order p from  p = |ln|e32/e21| + q(p)| / ln r21,  q = ln((r21^p - s)/(r32^p - s)),
  s = sign(e32/e21), e21 = phi2 - phi1, e32 = phi3 - phi2;
  extrapolated value phi_ext = (r21^p phi1 - phi2)/(r21^p - 1);
  GCI_fine = 1.25 |phi1 - phi2|/|phi1| / (r21^p - 1).
If e32/e21 < 0 the convergence is oscillatory and if |e32/e21| < 1 (the change on the finer pair
is the larger one) it diverges; no order is then
reported, and the spread of the three values is given instead.
Output: figures/fsLIFT_332mJ_Ti_meshStudy/data/gci_<t>.csv and a Markdown table on stdout.
"""
import glob, math, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures", "fsLIFT_332mJ_Ti_meshStudy"); os.makedirs(os.path.join(OUT, "data"), exist_ok=True)
work = sys.argv[1]; T = float(sys.argv[2]) if len(sys.argv) > 2 else 8.02e-10
LEVELS = (3, 2, 1)  # finest first


def ncells(d):
    s = open(os.path.join(d, "log.blockMesh")).read()
    return int(re.search(r"nCells:\s*(\d+)", s).group(1))


def logs(d):
    return [f for st in ("stage1", "stage2", "stage3") for f in sorted(glob.glob(os.path.join(d, "log.%s.part*" % st))) + [os.path.join(d, "log." + st)] if os.path.exists(f)]


def pocket(d, tmax):
    """peak trapped-vapour pressure p_max and its time, up to tmax"""
    best = (0.0, 0.0)
    pat = re.compile(r"vapourPocket\(local\): t = ([0-9.eE+-]+) .*p_max = ([0-9.eE+-]+) Pa")
    for f in logs(d):
        for line in open(f, errors="ignore"):
            m = pat.search(line)
            if m:
                t, p = float(m.group(1)), float(m.group(2))
                if t <= tmax * (1 + 1e-6) and p > best[0]:
                    best = (p, t)
    return best


def tracking(d, tmax):
    import csv
    te = tl = 0.0
    with open(os.path.join(d, "liftProcessTracking.csv")) as fh:
        for row in csv.DictReader(fh):
            try:
                t = float(row["time_ps"]) * 1e-12
            except ValueError:
                continue
            if t <= tmax * (1 + 1e-6):
                te = max(te, float(row["Te_max_K"])); tl = max(tl, float(row["Tl_max_K"]))
    return te, tl


def fields(d, t):
    D = np.genfromtxt(os.path.join(d, "mesh_study_fields.csv"), delimiter=",", names=True)
    i = np.argmin(abs(D["t_s"] - t))
    if abs(D["t_s"][i] - t) > 0.02 * t:
        raise SystemExit("%s: no written time near %g s" % (d, t))
    return D[i]


def gci(phi1, phi2, phi3, r21, r32):
    e21, e32 = phi2 - phi1, phi3 - phi2
    if e21 == 0 or e32 == 0:
        return dict(kind="no change", p=np.nan, ext=phi1, gci=0.0)
    ratio = e32 / e21
    if ratio < 0:
        return dict(kind="oscillatory", p=np.nan, ext=np.nan, gci=np.nan)
    s = 1.0
    p = abs(math.log(abs(ratio))) / math.log(r21)
    for _ in range(200):
        q = math.log((r21 ** p - s) / (r32 ** p - s))
        pn = abs(math.log(abs(ratio)) + q) / math.log(r21)
        if abs(pn - p) < 1e-10:
            break
        p = pn
    if abs(ratio) < 1:
        return dict(kind="diverging", p=p, ext=np.nan, gci=np.nan)
    ext = (r21 ** p * phi1 - phi2) / (r21 ** p - 1)
    g = 1.25 * abs((phi1 - phi2) / phi1) / (r21 ** p - 1)
    return dict(kind="monotonic", p=p, ext=ext, gci=g)


dirs = [os.path.join(work, "L%d" % L) for L in LEVELS]
N = [ncells(d) for d in dirs]
h = [N_ ** -0.5 for N_ in N]
r21, r32 = h[1] / h[0], h[2] / h[1]
Q = {}
for d in dirs:
    F = fields(d, T); pk = pocket(d, T); te, tl = tracking(d, T)
    for k, v in (("speed of the metal, r < 1 um [m/s]", -F["v_r1_ms"]),
                 ("speed of the metal, r < 2.5 um [m/s]", -F["v_r25_ms"]),
                 ("axial momentum of the metal [kg m/s]", -F["Py_kgms"]),
                 ("kinetic energy of the metal [J]", F["KE_J"]),
                 ("displacement of the film on the axis [nm]", (20.0 - F["y_bot_axis_um"]) * 1e3),
                 ("metal thickness on the axis [nm]", F["h_axis_nm"]),
                 ("peak trapped-vapour pressure [MPa]", pk[0] / 1e6),
                 ("time of the peak pressure [ps]", pk[1] * 1e12),
                 ("peak electron temperature [K]", te),
                 ("peak lattice temperature [K]", tl),
                 ("metal volume [um3]", F["V_metal_um3"])):
        Q.setdefault(k, []).append(v)
print("cells (fine, medium, coarse): %s;  r21 = %.4f, r32 = %.4f;  t = %g s\n" % (N, r21, r32, T))
print("| Quantity | fine (L3) | medium (L2) | coarse (L1) | convergence | apparent order p | extrapolated | GCI fine |")
print("|---|---|---|---|---|---|---|---|")
rows = []
for k, (f1, f2, f3) in Q.items():
    g = gci(f1, f2, f3, r21, r32)
    fmt = (lambda x: "%.4g" % x)
    print("| %s | %s | %s | %s | %s | %s | %s | %s |" % (k, fmt(f1), fmt(f2), fmt(f3), g["kind"],
          "–" if np.isnan(g["p"]) else "%.2f" % g["p"], "–" if np.isnan(g["ext"]) else fmt(g["ext"]),
          "–" if np.isnan(g["gci"]) else "%.2g %%" % (100 * g["gci"])))
    rows.append((k, f1, f2, f3, g["kind"], g["p"], g["ext"], g["gci"]))
with open(os.path.join(OUT, "data", "gci_%gs.csv" % T), "w") as fh:
    fh.write("# cells fine,medium,coarse = %s; r21 = %.5f; r32 = %.5f\n" % (N, r21, r32))
    fh.write("quantity,fine,medium,coarse,convergence,p,extrapolated,GCI_fine\n")
    for r_ in rows:
        fh.write('"%s",%.6g,%.6g,%.6g,%s,%.4g,%.6g,%.4g\n' % r_)
