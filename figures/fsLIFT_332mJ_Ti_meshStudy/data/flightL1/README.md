# Coarse mesh L1 continued after 315 ps (fixes 8 and 9)

Run: mesh level 1 of `cases/fsLIFT_332mJ_Ti_meshStudy` (62,010 cells) from t = 0 with both fixes of
docs/CODE_REVIEW.md §1 and `energyDiagnostics true`. Stages: 0–2 ps at maxDeltaT 1e-14 s, to 30 ps
at 5.95e-13 s, to 315 ps at 5.95e-13 s (written every 2.5 ps), to 1 ns at 1e-12 s (every 20 ps),
then at 2.5e-12 s (every 0.1 ns) to 30 ns (last write 29.95 ns). It matches the Phase D L1 run within 0.04 % at all common times.

`fields.csv`: `tools/mesh_study_fields.py` output (full-circle values), 1 ps to 7.595 ns. After about
7.2 ns the cap on the axis is thinner than one cell (71 nm cells for y < 16.5 µm), so
`y_bot_axis_um` is NaN from 7.195 ns and `h_axis_nm` loses its meaning.

`front.csv`: `tools/flight_front.py` output, about every 1 ns from 7 ns to 18.7 ns. The column
`alpha_max_axis_below18p5` stays at 0.78 from about 9 ns: no cell of the flying cap near the axis
is full of metal, so the cap is not resolved on this mesh and these rows show where the metal is,
not its shape.

Main values: speed within 1 µm of the axis 432.9 m/s at 1.9 ns (highest), 413 m/s at 6.5 ns;
first break of the dome wall between the writes at 3.96 and 4.26 ns; lowest metal at 10.52 µm
(2.52 µm above the receiver) at 18.7 ns.

`momentum_by_alpha_and_radius.csv`: the film momentum (10⁻¹³ kg m/s) split by α class and by radius, 2.96–7.9 ns.
`force_decomposition_*.csv`: `tools/force_decomposition.py` output before (0.4–3.9 ns) and after (3.96–7.7 ns) the tear; see docs/INVESTIGATION_LOG.md §14. The momentum rise after the tear is mostly a numerical net surface-tension force plus the unvented pocket, so the momentum after about 4 ns is not physical.

`full_run_curves.csv`: `tools/full_run_curves.py` output at the 362 frame times of the whole-run
animation (1 ps to 29.95 ns): gap between the lowest metal and the receiver, α-weighted speed and
lattice temperature of the metal within 0.5 µm of the axis. `gap_um` = 0.106 means contact (the cell
layer touching the receiver is not counted). The jet tip reaches the receiver between 23.63 and
23.73 ns at about 420 m/s; afterwards the liquid receiver is cratered, which is not physical
(docs/INVESTIGATION_LOG.md §15).
