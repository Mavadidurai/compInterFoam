"""Figures of the mesh-study case with fixes 8 and 9 (docs/RESULTS.md section 1).

    python3 tools/plot_mesh_study_flight.py

Reads figures/fsLIFT_332mJ_Ti_meshStudy/data/ and writes into figures/fsLIFT_332mJ_Ti_meshStudy/:
  mesh_study_315ps.png   film momentum and speed within 1 um of the axis on L1/L2/L3, 0-315 ps
  fix_effect_315ps.png   L1 momentum at 315 ps before the fixes, with fix 8, with fixes 8 and 9
  flight_L1_history.png  L1 continued after 315 ps to 30 ns: speed, momentum (to 7.6 ns) and the lowest metal
"""
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F = os.path.join(ROOT, 'figures', 'fsLIFT_332mJ_Ti_meshStudy')
D = os.path.join(F, 'data')
plt.rcParams.update({'font.size': 12})


def load(p):
    return np.genfromtxt(p, delimiter=',', names=True)


# 1. mesh study at 315 ps
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for L, c, n in ((1, '#9ecae1', '62,010'), (2, '#3182bd', '124,200'), (3, '#08306b', '248,252')):
    d = load(os.path.join(D, 'fixes8and9', 'fields_L%d.csv' % L))
    t = d['t_s'] * 1e12
    ax[0].plot(t, -d['Py_kgms'] * 1e13, 'o-', color=c, ms=3.5, label='L%d (%s cells)' % (L, n))
    ax[1].plot(t, -d['v_r1_ms'], 'o-', color=c, ms=3.5)
ax[0].set_ylabel('film momentum towards the receiver\n(10$^{-13}$ kg m/s, full circle)')
ax[1].set_ylabel('speed of the metal within 1 µm\nof the axis (m/s)')
ax[0].legend(frameon=False)
for a in ax:
    a.set_xlabel('time (ps)')
    a.set_xlim(0, 320)
    a.set_ylim(bottom=0)
    a.spines[['top', 'right']].set_visible(False)
fig.suptitle('Mesh study with fixes 8 and 9: momentum within 0.6 %, speed within 0.03 % (L2-L3) at 315 ps', fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(F, 'mesh_study_315ps.png'), dpi=150)
plt.close(fig)

# 2. effect of fixes 8 and 9 on L1 (values from docs/INVESTIGATION_LOG.md sections 2, 6 and 9)
fig, ax = plt.subplots(figsize=(7, 4.2))
lab = ['before the fixes\n(Phase A)', 'fix 8: temperatures carried\nwith the metal (Phase C)', 'fixes 8 and 9: pocket on the\ndonor film only (Phase D)']
val = [5.010, 7.388, 9.116]
b = ax.bar(lab, val, color=['0.7', '#6baed6', '#08519c'], width=0.6)
ax.bar_label(b, fmt='%.3f')
ax.set_ylabel('L1 film momentum at 315 ps\n(10$^{-13}$ kg m/s)')
ax.set_ylim(0, 10.5)
ax.spines[['top', 'right']].set_visible(False)
ax.tick_params(axis='x', labelsize=10)
fig.tight_layout()
fig.savefig(os.path.join(F, 'fix_effect_315ps.png'), dpi=150)
plt.close(fig)

# 3. flight on L1
d = load(os.path.join(D, 'flightL1', 'fields.csv'))
fr = load(os.path.join(D, 'flightL1', 'front.csv'))
fc = load(os.path.join(D, 'flightL1', 'full_run_curves.csv'))
t = d['t_s'] * 1e9
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
ax[0].plot(t, -d['v_r1_ms'], color='#238b45', lw=2, label='within 1 µm of the axis')
ax[0].plot(fr['t_ns'], -fr['v_below18p5_ms'], 's', color='#74c476', label='all metal below y = 18.5 µm')
ax[0].plot(fc['t_ns'][fc['t_ns'] > 7.0], fc['v_axis_ms'][fc['t_ns'] > 7.0], color='#00441b', lw=1.5, ls='--', label='within 0.5 µm of the axis (jet)')
ax[0].set_ylabel('speed towards the receiver (m/s)')
ax[0].legend(frameon=False, fontsize=10)
ax[1].plot(t, -d['Py_kgms'] * 1e13, color='#08519c', lw=2)
ax[1].set_ylabel('film momentum (10$^{-13}$ kg m/s)')
ax[2].plot(t, d['y_bot_axis_um'], color='k', lw=2, label='lower film surface on the axis')
fg = fc['t_ns'] > 7.0
ax[2].plot(fc['t_ns'][fg], 8.0 + fc['gap_um'][fg], color='#d94801', lw=2, label='lowest metal (α > 0.5)')
ax[2].axvline(23.7, color='#d94801', ls=':', lw=1.2)
ax[2].text(23.2, 13.5, 'contact\n23.6–23.7 ns', color='#d94801', fontsize=10, ha='right')
ax[2].axhline(8.0, color='0.4', ls='--')
ax[2].text(0.3, 8.15, 'receiver surface', color='0.4', fontsize=10)
ax[2].axhspan(8.0, 16.5, color='#fdd0a2', alpha=0.35, lw=0)
ax[2].text(0.3, 12.0, '71 nm cells on L1:\ncap thinner than one cell', color='#a63603', fontsize=10)
ax[2].set_ylabel('height above the receiver base (µm)')
ax[2].set_ylim(7, 20.5)
ax[2].legend(frameon=False, fontsize=10, loc='upper right')
for a in ax:
    a.set_xlabel('time (ns)')
    a.set_xlim(0, 31)
    a.axvspan(3.96, 4.26, color='0.85', lw=0)
    a.spines[['top', 'right']].set_visible(False)
ax[0].set_ylim(0, 480)
ax[1].set_ylim(0, 11.5)
ax[1].text(4.4, 1.0, 'grey band: first\nbreak of the dome wall', fontsize=10, color='0.35')
fig.suptitle('Coarse mesh L1 (62,010 cells) continued after 315 ps: launch, tear, flight and impact', fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(F, 'flight_L1_history.png'), dpi=150)
plt.close(fig)
print('written to', F)
