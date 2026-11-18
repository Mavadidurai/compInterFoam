# Mesh study with temperatures carried by the metal flux

Runs: L1/L2/L3 (62,010 / 124,200 / 248,252 cells, refinement ratio sqrt(2)),
maxDeltaT 5.95e-13 s, `temperatureAdvectionFlux metal`, restarted from the
2 ps states of the earlier study, run to 3.148e-10 s.

`fields_L*.csv`: tools/mesh_study_fields.py output (full-circle values).
`gci_3.148e-10s.csv`: three-level Richardson/GCI at 315 ps.

The files one level up (`fields_L*.csv`, `gci_1e-09s.csv`) were made with
the earlier mixture-flux transport and are superseded: that transport cooled
the hottest metal cells at a mesh-dependent rate and biased momentum low by
about 30-50 %.
