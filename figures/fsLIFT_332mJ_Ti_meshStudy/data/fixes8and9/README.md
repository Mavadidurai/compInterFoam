# Mesh study with fixes 8 and 9

Runs: L1/L2/L3 (62,010 / 124,200 / 248,252 cells, refinement ratio sqrt(2)),
all from t = 0 with both fixes of docs/CODE_REVIEW.md §1:
temperatures carried by the metal flux (fix 8) and the confined-vapour
pocket acting on the donor film only (fix 9). 0-2 ps at maxDeltaT 1e-14 s,
then maxDeltaT 5.95e-13 s to 3.148e-10 s.

`fields_L*.csv`: tools/mesh_study_fields.py output (full-circle values).
`gci_3.148e-10s.csv`: tools/mesh_study_gci.py at 315 ps. The peak pocket
pressure and its time come from the `vapourPocket` log lines, which these
runs print only a few dozen times in 315 ps; those two rows are limited by
that sampling on every level and are not a mesh comparison.

At 315 ps the axial momentum is 9.116 / 9.104 / 9.061 e-13 kg m/s
(L1/L2/L3): spread 0.6 %, but the change grows on the finer pair, so no
order or GCI can be given. The momentum budget closes to 0.969 (L2) and
0.977 (L3) of the film's momentum change over 2-315 ps (L1: 0.959 over 0-305 ps).

This supersedes `../metalFlux/` (fix 8 only), which counted the receiver push
that fix 9 removed and is about 20-25 % low.
