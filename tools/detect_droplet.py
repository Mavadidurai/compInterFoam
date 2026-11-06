#!/usr/bin/env python3
"""Scan postProcessing/axisLine/<t>/axis_*.xy and report metal segments on the axis.

A droplet counts as detached when, on the beam axis, a metal segment (alpha>0.5)
inside the gap (8.05 um < y < 20.1 um) is separated from every other metal
segment by gas (alpha<0.5) over at least `gap_nm` nanometres.
Usage: detect_droplet.py CASE [--all]
Prints one line per sampled time (with --all) or only the latest, and exits 0 if
a detached droplet is present at the latest time, 1 otherwise.

Report-only mode: if the file CASE/droplet_report_only exists, detections are still
printed (marked "<-- DETACHED") and appended to CASE/droplet_detections.txt, but the
exit status is always 1, so a calling script never shortens the run. Use this mode
whenever the film may be smeared: a smeared film can show two "segments" that are
not really separated. Always confirm a detection in ParaView.
"""
import glob, os, sys

case = sys.argv[1]
show_all = "--all" in sys.argv
GAP_NM = 20.0


def segments(path):
    # raw set files are named axis_<f1>_<f2>_..._.xy; scalar columns follow y in that order
    names = os.path.basename(path)[len("axis_"):-len(".xy")].split("_")
    ia = 1 + names.index("alpha.metal")
    ys, a = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.split()
            ys.append(float(c[0])); a.append(float(c[ia]))
    segs, start = [], None
    for y, v in zip(ys, a):
        if v > 0.5 and start is None:
            start = y
        if v <= 0.5 and start is not None:
            segs.append([start, prev]); start = None
        prev = y
    if start is not None:
        segs.append([start, prev])
    # keep only segments in the receiver gap + film (8.05-20.1 um): drops the receiver
    # substrate and any fragment on the donor (glass) side, which must not count as a droplet
    segs = [s for s in segs if s[1] > 8.05e-6 and s[0] < 20.1e-6]
    merged = []
    for s in segs:
        if merged and (s[0] - merged[-1][1]) * 1e9 < GAP_NM:
            merged[-1][1] = s[1]
        else:
            merged.append(s)
    return merged


dirs = sorted(glob.glob(os.path.join(case, "postProcessing", "axisLine", "*")),
              key=lambda d: float(os.path.basename(d)))
if not dirs:
    print("no samples yet"); sys.exit(1)
detached = False
for d in (dirs if show_all else dirs[-1:]):
    f = glob.glob(os.path.join(d, "axis_*alpha.metal*.xy"))
    if not f:
        continue
    m = segments(f[0])
    detached = len(m) >= 2
    desc = "; ".join("%.4f-%.4f um (%.0f nm)" % (s[0]*1e6, s[1]*1e6, (s[1]-s[0])*1e9) for s in m)
    print("t=%.4g ps  metal segments on axis: %d  [%s]%s" % (
        float(os.path.basename(d))*1e12, len(m), desc, "  <-- DETACHED" if detached else ""))
if os.path.exists(os.path.join(case, "droplet_report_only")):
    if detached:
        with open(os.path.join(case, "droplet_detections.txt"), "a") as log:
            log.write("t=%.6g ps detached segments: %s\n" % (float(os.path.basename(dirs[-1]))*1e12, desc))
    sys.exit(1)
sys.exit(0 if detached else 1)
