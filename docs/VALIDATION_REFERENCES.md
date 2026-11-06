# Validation and verification references

This page records candidate literature from the supplied `fsLIFT_literature_register.xlsx` for planning the next validation and verification work on compInterFoam-fsLIFT. It is a curated index, not a claim that the referenced datasets have been obtained or that every citation has been independently verified.

## Current project status

See [VERIFICATION.md](VERIFICATION.md) for the verification cases already documented in this repository. It reports code checks and selected drop benchmarks, while identifying two relevant gaps: no mesh-refinement study of the production fs-LIFT case and no comparison with experimental measurements. This page focuses on references that can help address those gaps.

## Numerical verification and solver benchmarks

| Reference | Use in this project | Source |
|---|---|---|
| Deshpande, Anumolu & Trujillo (2012), “Evaluating the performance of the two-phase flow solver interFoam” | InterFoam advection, high-Weber-number and surface-tension benchmark context. Numerical verification, not fs-LIFT experimental validation. | [DOI: 10.1088/1749-4699/5/1/014016](https://doi.org/10.1088/1749-4699/5/1/014016) |
| Scheufler & Roenby (2023), “TwoPhaseFlow: A Framework for Developing Two Phase Flow Solvers in OpenFOAM” | OpenFOAM two-phase solver methods, curvature and phase-change benchmarking context. | [DOI: 10.51560/ofj.v3.80](https://doi.org/10.51560/ofj.v3.80) |
| Lippert et al. (2022), “A benchmark for surface-tension-driven incompressible two-phase flows” | Static, translating and oscillating drop cases; useful independent benchmark targets for interface and surface-tension numerics. | [arXiv record](https://doi.org/10.48550/arXiv.2212.02904) |
| Vachaparambil & Einarsrud (2019), “Comparison of Surface Tension Models for the Volume of Fluid Method” | Compare CSF-family surface-tension treatments and spurious-current behavior. | [DOI: 10.3390/pr7080542](https://doi.org/10.3390/pr7080542) |
| Roache (1994), “Perspective: A Method for Uniform Reporting of Grid Refinement Studies” | Reporting framework for grid convergence and solution verification. | [DOI: 10.1115/1.2910291](https://doi.org/10.1115/1.2910291) |

## Experimental LIFT evidence

| Reference | Potential comparison observables | Source |
|---|---|---|
| Kuznetsov, Unger, Koch et al. (2012), “Laser-induced jet formation and droplet ejection from thin metal films” | Jet and droplet formation, timing and morphology; confirm material and pulse conditions before quantitative comparison. | [DOI: 10.1007/s00339-011-6747-8](https://doi.org/10.1007/s00339-011-6747-8) |
| Klini, Loukakos, Gray et al. (2008), “Laser Induced Forward Transfer of metals by temporally shaped femtosecond laser pulses” | Transfer morphology and pulse-shaping effects; compare only when pulse and donor conditions match. | [DOI: 10.1364/OE.16.011300](https://doi.org/10.1364/OE.16.011300) |
| Banks, Grivas, Mills, Eason & Zergioti (2006), “Nanodroplets deposited in microarrays by femtosecond Ti:sapphire laser-induced forward transfer” | Transferred feature size and morphology for femtosecond Cr-film LIFT; check availability of underlying measurements. | [DOI: 10.1063/1.2386921](https://doi.org/10.1063/1.2386921) |
| Pohl, Visser, Römer et al. (2015), “Ejection Regimes in Picosecond Laser-Induced Forward Transfer of Metals” | Regime-map and ejection observations; picosecond conditions differ from fs-LIFT and should be treated as qualitative context unless matched. | [DOI: 10.1103/PhysRevApplied.3.024001](https://doi.org/10.1103/PhysRevApplied.3.024001) |
| Alloncle, Bouffaron & Hermann et al. (2006), “Experimental study of front and back ablation of metal thin film using ultrashort laser pulses” | Ablation mode, fluence and film-thickness context; precursor physics rather than direct droplet validation. | [DOI: 10.1117/12.674683](https://doi.org/10.1117/12.674683) |

## Candidate observables and comparison controls

The supplied register recommends tracking film displacement versus time, axial velocity, dome height, radial extent, sheet or neck thickness, separation time, droplet diameter and transferred volume. For each comparison, record the measurement definition and uncertainty alongside the simulation extraction method.

Before comparing quantitative values, record at minimum: material and film thickness, wavelength, pulse duration and shape, fluence definition and value, donor-receiver gap, ambient pressure, geometry, and time origin. Do not combine values across different conditions as if they represented one validation case.

## Data acquisition status

The workbook supplies DOI and publisher or discovery links; it does not itself contain experimental time-series data. A DOI landing page is a citation/source record, not evidence that raw measurements are downloadable. No external data files are included here because availability and reuse terms have not been established from the register alone.

For each candidate paper, check the publisher page and any linked supplementary information or public repository. If a reusable dataset is available, add the actual data files or a stable repository link here with its license, version, units, and mapping to the observable table above. If the paper provides only figures, digitized values should be labeled as figure-derived and accompanied by the digitization method and uncertainty.

## Source workbook notes

The supplied register's QA sheet says DOI fields are populated only when supported by known or publisher bibliographic records, while several rows remain explicitly unverified. Confirm exact title, author list, year and publication status before citing a row in a report. The register also cautions against treating titanium ablation or VOF-method papers as direct experimental validation of titanium fs-LIFT.
