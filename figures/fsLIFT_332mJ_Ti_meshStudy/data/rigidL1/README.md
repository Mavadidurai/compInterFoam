# Rigid receiver: run A's jet hitting a wall (coarse mesh L1)

Run: run A (`../flightL1/`) restarted at 21.5 ns, when the lowest metal was 1.17 µm above the receiver. The 1,802 receiver cells (y < 8 µm) were removed, and the exposed faces became a wall:
- no slip;
- near 300 K (`T` fixedValue 300 K; `Tl`, `Te` mixed with valueFraction 0.9);
- zero-gradient `alpha.metal`, i.e. a 90° contact angle.

The run went to 29.98 ns, written every 0.1 ns, on 60,208 cells. Set-up and results: docs/INVESTIGATION_LOG.md §16. The metal does not solidify in this model.

`impact_rigid.csv`: `tools/impact_rigid.py` output at every write, 21.5–29.98 ns. Metal means α > 0.5, and volumes are full-circle µm³.
- Contact is between 23.9 and 24.0 ns: `gap_um` = 0.035 means the first cell layer.
- `p_wall_max_MPa` peaks at 508 MPa at 24.1 ns.
- `r_wall_um` returns to 0 after 27.3 ns: the film has become thinner than half a cell, so no wall cell is more than half metal. Use `wall_film_profile.csv` for its extent.
- `V_wall_um3` includes about 0.007–0.009 µm³ of cold receiver metal left above the cut at 21.5 ns.

`wall_film_profile.csv`: `tools/wall_film_profile.py` output at the 85 animation frames, 21.6–29.98 ns.
- Each row is `t_ns` followed by the equivalent film thickness h in nm, at the radial cell centres (µm) given in the header.
- h is the sum of α Δy over the two cell layers next to the wall, with metal below 400 K (the leftover receiver metal) excluded.

`full_run_curves_rigid.csv`: the curves of the whole-run animation, from `../flightL1/full_run_curves.csv` up to 21.5 ns, then `tools/full_run_curves.py` on the rigid run.
