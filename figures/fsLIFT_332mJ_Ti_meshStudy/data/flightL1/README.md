# Coarse mesh L1 continued after 315 ps (fixes 8 and 9)

Run: mesh level 1 of `cases/fsLIFT_332mJ_Ti_meshStudy` (62,010 cells) from t = 0 with both fixes of
docs/CODE_REVIEW.md §1 and `energyDiagnostics true`. Stages: 0–2 ps at maxDeltaT 1e-14 s, to 30 ps
at 5.95e-13 s, to 315 ps at 5.95e-13 s (written every 2.5 ps), to 1 ns at 1e-12 s (every 20 ps),
then at 2.5e-12 s (every 0.1 ns). It matches the Phase D L1 run within 0.04 % at all common times.

`fields.csv`: `tools/mesh_study_fields.py` output (full-circle values), 1 ps to 7.595 ns. After about
7.2 ns the cap on the axis is thinner than one cell (71 nm cells for y < 16.5 µm), so
`y_bot_axis_um` is NaN from 7.195 ns and `h_axis_nm` loses its meaning.

`front.csv`: `tools/flight_front.py` output, about every 1 ns from 7 ns to 18.7 ns. The column
`alpha_max_axis_below18p5` stays at 0.78 from about 9 ns: no cell of the flying cap near the axis
is full of metal, so the cap is not resolved on this mesh and these rows show where the metal is,
not its shape.

Main values: speed within 1 µm of the axis 432.9 m/s at 1.9 ns (highest), 413 m/s at 6.5 ns;
first break of the dome wall between the writes at 3.96 and 4.26 ns; lowest metal at 10.52 µm
(2.52 µm above the receiver) at 18.7 ns. The run continues to 30 ns; this directory will be
extended when it finishes.
