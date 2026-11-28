# Run C (refined flight path): data

`full_run_curves_C.csv`: `tools/full_run_curves.py` on every write of run C (`flightC/`), same columns and definitions as `../flightL1/full_run_curves.csv` (run A).
Interim: written on 27 Nov at 21:14 UTC, from 2.9625 ns (the mapped L1 state) to 7.4774 ns. It will be replaced when C reaches 35 ns.

Comparison with run A at the same times: [docs/INVESTIGATION_LOG.md §13](../../../../docs/INVESTIGATION_LOG.md#13-run-c-refined-flight-path-22-nov).

Update 28 Nov (C at 24.47 ns): `full_run_curves_C.csv` now runs to 24.47 ns. `front_tip_C.csv` is `tools/front_tip_stats.py` at 12 times from 12 to 24.47 ns (columns in the tool's docstring). On mesh C the first counted cell above the receiver has its centre at 8.0525 µm, so `gap_um` = 0.0525 means contact (0.106 on L1). The tip reaches the receiver at 20.6 ns (23.6–23.7 ns on L1).
