# fsLIFT_332mJ_axisymmetric: main results

| Quantity | Value |
|---|---|
| Simulated time | 19.97 ns of 20 ns (finished) |
| Highest speed of the film on the beam axis | 191.7 m/s at 0.54 ns |
| Speed of the metal on the beam axis at the end (after about 10 ns the gas behind the film pushes the detached piece, so this late value is too high; see docs/RESULTS.md) | 126.4 m/s |
| Weber number at the highest speed (film inertia divided by surface tension, ρv²h/σ) | 7.2  (= 4515 · 191.7² · 71.4e-9 / 1.64) |
| Kinetic energy per area divided by the energy of the two new surfaces (½ρhv² / 2σ) | 5.92 / 3.28 J/m² = 1.81 |
| Distance moved by the film on the beam axis at 1 / 5 / 20 ns | 0.18 µm / 0.90 µm / 2.16 µm (at 19.97 ns) |
| Separate metal pieces crossed by the beam axis at the end (the axis detector only sees a gas gap along the axis; a piece that separates by tearing off the axis is counted in docs/RESULTS.md) | 1 |
| Detections by the axis droplet detector | 0 |
| Change of total metal volume (includes metal turned to vapour) | -0.0135 % |
| Pressure of the vapour trapped on the glass side: peak / at the end | 524.0 MPa / 0.000 MPa |

The volumes and speeds of the separated pieces come from tools/count_metal_pieces.py; see docs/RESULTS.md.
