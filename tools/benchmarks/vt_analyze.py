#!/usr/bin/env python3
"""Decay rate of the J0 mode in a viscous test case (make_vt.py) vs the exact nu k^2.  Usage: vt_analyze.py <case>"""
import sys, os, re, numpy as np
c = sys.argv[1]; nk2 = 2.35e-3 / 4515 * (8.653727912911013 / 0.6e-6)**2
def readU(f):
    s = open(f).read(); m = re.search(r'internalField\s+nonuniform\s+List<vector>\s*\n?(\d+)\s*\n?\(', s); n = int(m.group(1))
    return np.array([list(map(float, t.split())) for t in re.findall(r'\(([^)]*)\)', s[m.end():])[:n]])[:, 1]
U0 = readU(c + '/0/U')
for t in sorted([d for d in os.listdir(c) if re.match(r'^[0-9.e-]+$', d) and d != '0'], key=float):
    A = (readU(c + '/' + t + '/U') * U0).sum() / (U0 * U0).sum(); tt = float(t)
    print('t=%.3g ns  amplitude %.5f  exact %.5f  rate/exact %.4f' % (tt * 1e9, A, np.exp(-nk2 * tt), -np.log(A) / tt / nk2))
