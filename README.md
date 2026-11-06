# compInterFoam-fsLIFT

An OpenFOAM v2406 solver for **femtosecond laser-induced forward transfer (fs-LIFT)** of a thin metal film. It is written for, and tested on, a 71.4 nm titanium film on glass.

In LIFT, a laser pulse is fired through a transparent donor plate (here glass) into a thin metal film coated on it. The pulse heats the film, and metal is pushed off the plate towards a receiver surface a few micrometres away. Whether the metal leaves as one clean droplet, as a spray, or not at all depends on the laser energy and on how the film heats, melts, evaporates and moves. This solver simulates that sequence from the arrival of the pulse (femtoseconds) to the flight of the ejected metal (tens of nanoseconds).

The solver is built on OpenFOAM's `compressibleInterFoam`, which tracks two immiscible phases (here metal and gas) with the volume-of-fluid method. It adds:

| Added model | What it represents |
|---|---|
| Laser heat source | A Gaussian pulse in time and space, absorbed in the first ~10 nm of the film |
| Two-temperature model | Separate temperatures for the electrons (heated first by the laser) and for the atoms of the metal (heated by the electrons within picoseconds) |
| Evaporation and recoil pressure | Metal evaporates from hot surfaces; the departing vapour pushes back on the surface |
| Trapped-vapour model | Vapour formed between the glass and the film cannot escape; its pressure pushes the film away from the glass |

---

## Status: what has been checked and what has not

**Verified:** each added equation has been checked against an exact solution, an independent code or a hand calculation. The tests, their results and how to repeat them are in [docs/VERIFICATION.md](docs/VERIFICATION.md). Examples:
- laser energy deposited in the film: within 0.2 % of the requested energy;
- the two-temperature model: within about 1 % of an independent finite-difference code;
- surface tension: pressure jump of a static drop within 2.75 % of the Laplace law at 10 cells per drop radius, and within 4.1 % at 40 cells (the error settles near −4 % rather than vanishing, a known limitation of the method);
- a film moving as a whole at 170 m/s keeps its velocity and pressure (Galilean invariance), and the metal volume is conserved exactly.

**Not validated:** the results have **not** been compared with measurements in this repository. "Verified" means the code solves its equations correctly. It does not mean the equations describe the real experiment, and the model leaves out known physics (see [Limitations](docs/RESULTS.md#limitations)).

## The reference result in short

Case `cases/fsLIFT_332mJ_axisymmetric`: titanium film 71.4 nm thick, laser 343 nm wavelength, 214 fs pulse, beam radius 3 µm, peak fluence 332 mJ/cm² before the glass. The gap to the receiver is 12 µm. Full description and figures: [docs/RESULTS.md](docs/RESULTS.md).

1. **0–2 ps:** the laser heats the electrons in the first ~10 nm of the film. They peak at about 16,400 K at 0.54 ps, and the atoms reach the same temperature by about 0.6 ps. This is above the critical temperature of titanium, where a liquid–vapour description no longer strictly holds (see [Limitations](docs/RESULTS.md#limitations)).
2. **Up to 0.6 ns:** the film surface facing the glass evaporates. The trapped vapour reaches a peak pressure of 524 MPa at 26 ps and has fallen below 1 MPa by 0.59 ns.
3. **0.5 ns:** the centre of the film has been pushed into a dome moving away from the glass at up to **192 m/s**. Its kinetic energy is 1.8 times the energy needed to create the two new surfaces of a detached piece.
4. **6–13 ns:** the dome deepens into a cup whose side walls thin from about 107 nm to 5 nm.
5. **About 13 ns:** the walls tear, and the cup bottom flies on as **one main droplet** at 70–90 m/s, with small ring-shaped fragments around it. At 19.8 ns the main droplet holds 0.25 µm³ (the volume of a 0.79 µm sphere). Nothing reaches the receiver within 20 ns, so there is no splash.

The tearing happens where the wall is thinner than one computational cell. The time of separation and the volumes therefore depend on the mesh and are **not yet mesh-converged**. After about 10 ns, gas flowing in behind the film (a model artefact) also pushes the separated metal. Read [docs/RESULTS.md](docs/RESULTS.md) before quoting any number.

![Metal surface at 1, 6, 13.4 and 19.8 ns, drawn in 3-D by revolving the axisymmetric result](figures/fsLIFT_332mJ_axisymmetric/gallery/montage_3d.png)

*The simulation is two-dimensional (axisymmetric, see below). The 3-D pictures are the 2-D result revolved about the beam axis.*

---

## Quick start

Requirements:
- OpenFOAM v2406 (the ESI/OpenCFD release);
- Python 3 with NumPy and Matplotlib for the analysis scripts;
- ParaView 5.10 or later (`pvpython`) for the images.

```bash
source /usr/lib/openfoam/openfoam2406/etc/bashrc     # adjust to your installation

# 1. build the solver (about 1 minute; installs into $FOAM_USER_APPBIN)
cd src && wmake && cd ..

# 2. run the reference case on 4 processes to 20 ns (about 20 hours on 4 cores)
tools/run_case.sh cases/fsLIFT_332mJ_axisymmetric 4 2e-8
RESUME=1 tools/run_case.sh cases/fsLIFT_332mJ_axisymmetric 4 2e-8    # continue after an interruption

# 3. results: film position, speed and temperature on the beam axis, results table, images
python3 tools/axis_timeseries.py cases/fsLIFT_332mJ_axisymmetric
python3 tools/results_package.py cases/fsLIFT_332mJ_axisymmetric
xvfb-run -a pvpython tools/render_gallery.py cases/fsLIFT_332mJ_axisymmetric
```

A short check that the build works is the 1-D laser column (about 10 minutes on one core). It is described in [docs/USER_GUIDE.md](docs/USER_GUIDE.md#3-the-one-dimensional-test-cases).

## Geometry of the reference case

The case is **axisymmetric**: every quantity is assumed to be the same at every angle around the beam axis. It is computed on a 5° wedge, a thin slice of the full cylinder:
- **y** runs along the beam axis;
- **r** is the distance from the axis;
- the laser travels in the −y direction.

```
 y [µm]
 32.07 ┌──────────────────────────────────────────┐
       │  donor glass (modelled as a gas)         │  ← laser enters here, travelling −y
 20.07 ├══════════════════════════════════════════┤  film surface on the glass side (vapour is trapped here)
       │  titanium film, 71.4 nm                  │
 20.00 ├══════════════════════════════════════════┤  film surface on the gap side (free surface)
       │  gap, 12 µm                              │
  8.00 ├──────────────────────────────────────────┤  receiver surface
       │  receiver                                │
  0.00 └──────────────────────────────────────────┘
       r = 0 (beam axis)                 r = 8 µm (open boundary)
```

## Repository layout

| Path | Content |
|---|---|
| `src/` | Solver source code and build files (`Make/`) |
| `cases/fsLIFT_332mJ_axisymmetric/` | Reference case: input dictionaries and initial fields (`0.orig`) |
| `cases/column1D_laser_TTM/` | 1-D column: laser heating and the two-temperature model, compared with an independent code |
| `cases/column1D_galilean_test/` | 1-D column moving as a whole at −170 m/s: velocity and pressure must not change |
| `tools/` | Run script, analysis and plotting scripts, independent finite-difference two-temperature code |
| `tools/benchmarks/` | Drop benchmarks (surface tension, drop oscillation, viscosity) |
| `figures/` | Figures and data of the reference result and of the verification tests |
| `docs/` | Documentation (below) |

## Documentation

| Document | Read it for |
|---|---|
| [docs/MODEL.md](docs/MODEL.md) | The physics and equations, the material and laser inputs, and what the solver does in each time step |
| [docs/USER_GUIDE.md](docs/USER_GUIDE.md) | Building, running, restarting, the case settings explained, the analysis tools, and how to set up a new case |
| [docs/VERIFICATION.md](docs/VERIFICATION.md) | Every test of the code, its result and how to repeat it, and the defects found and corrected during testing |
| [docs/RESULTS.md](docs/RESULTS.md) | The reference result in detail, with figures, and the limitations of the model |

## Licence

GNU General Public License version 3 or later (see `LICENSE`), the licence of OpenFOAM, from which this solver is derived. OpenFOAM is a registered trademark of OpenCFD Ltd. This software is not approved or endorsed by OpenCFD Ltd.

## Author

Mavadidurai Gopinathan.
