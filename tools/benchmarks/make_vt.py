#!/usr/bin/env python3
"""Viscous-operator verification: single-phase Ti, axisymmetric axial flow
u_y(r,t) = U0 J0(k r) exp(-nu k^2 t), J0(k R) = 0 at the no-slip wall r = R = 0.6 um (3rd zero).
Exact solution of the axisymmetric unsteady Stokes AND Navier-Stokes equations (u.grad u = 0).
Usage (case built by make_drop.sh, writeFormat ascii, after `postProcess -func writeCellCentres -time 0`):
    python3 make_vt.py <case>"""
import sys, re, numpy as np
c = sys.argv[1]
def read(f):
    s = open(f).read(); m = re.search(r'internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\n?\(', s)
    return np.array(s[m.end():].split(')')[0].split(), float)
r = np.hypot(read(c + '/0/Cx'), read(c + '/0/Cz')); n = len(r)
th = np.linspace(0, np.pi, 4001)
trap = getattr(np, 'trapezoid', None) or np.trapz   # NumPy < 2.0 has only trapz
J0 = lambda v: trap(np.cos(np.outer(v, np.sin(th))), th, axis=1) / np.pi
k = 8.653727912911013 / 0.6e-6; u = J0(k * r)
s = open(c + '/0/U').read()
s = re.sub(r'internalField\s+uniform\s*\([^)]*\)\s*;', 'internalField   nonuniform List<vector>\n%d\n(\n' % n + '\n'.join('(0 %.10e 0)' % v for v in u) + '\n);\n', s, count=1)
for pat, rep in ((r'right\s*\{[^}]*\}', 'right { type noSlip; }'), (r'substrate\s*\{[^}]*\}', 'substrate { type zeroGradient; }'), (r'donor\s*\{[^}]*\}', 'donor { type zeroGradient; }')):
    s = re.sub(pat, rep, s, count=1)
open(c + '/0/U', 'w').write(s)
p = open(c + '/0/p_rgh').read()
for pat, rep in ((r'right\s*\{[^}]*\}', 'right { type fixedFluxPressure; value uniform 1e5; }'), (r'substrate\s*\{[^}]*\}', 'substrate { type fixedValue; value uniform 1e5; }'), (r'donor\s*\{[^}]*\}', 'donor { type fixedValue; value uniform 1e5; }')):
    p = re.sub(pat, rep, p, count=1)
open(c + '/0/p_rgh', 'w').write(p)
for f, v in (('alpha.metal', 1), ('alpha.air', 0)):
    s = open(c + '/0/' + f).read(); s = re.sub(r'internalField\s+uniform\s+\d+;', 'internalField   uniform %d;' % v, s, count=1)
    s = re.sub(r'right\s*\{[^}]*\}', 'right { type zeroGradient; }', s, count=1); open(c + '/0/' + f, 'w').write(s)
print(n, 'cells; nu k^2 = %.4e 1/s' % (2.35e-3 / 4515 * k**2))
