# SEM Grating Pitch Benchmark Dataset

## Overview
This dataset accompanies the manuscript *"Self-consistent SEM magnification linearity assessment using a holographic grating: Fourier analysis, detector comparison, and open dataset"*.

The dataset comprises 1500 SEM images in total (750 per detector) of a Thorlabs GH13-06U reflective holographic diffraction grating (600 lines/mm nominal), acquired on a JEOL JSM-IT700HR field-emission SEM in four montage series using both secondary-electron (SED) and backscattered-electron (BED-S) detectors. The series are labelled 400×, 500×, 750× and 1000× after the survey magnification that defines the montage area: 400×: 19×19 grid (361 SED + 361 BED-S), 500×: 15×15 (225 + 225), 750×: 10×10 (100 + 100), 1000×: 8×8 (64 + 64).

Every image of every series was recorded with the same scan settings: field of view 17.07 µm × 12.80 µm (1280 × 960 px), which the image metadata record as an image magnification of 7500×. The series therefore differ in the scanned area and in the acquisition date, not in the image magnification, and the pixel scale of 0.01333594 µm/px applies to all images (see `revision/analysis/sem_metadata.csv`).

The raw SEM images (~2.5 GB) are not hosted in this repository owing to their size; they are available from the corresponding author upon reasonable request. This repository contains the complete analysis results and the R/Python analysis pipelines.

## Contents

```
# Per-image and aggregate results
per_image_pitch_mag{mag}_{detector}.csv   # Per-image 1D pitch + sub-bin profile-FFT pitch (synced from fft_subpixel_per_image.csv)
fft_subpixel_per_image.csv                # Sub-bin-refined 2D FFT pitch, all 1500 images
pitch1d_physical.csv                      # Per-image 1D pitch + stage coordinates, both detectors
linearity_physical_units.csv              # Linearity slopes in um/um (SED-only, superseded)
bootstrap_summary.csv                     # Bootstrap BCa confidence intervals
bland_altman_bootstrap_summary.csv        # Bootstrapped Bland-Altman statistics
groove_summary_mag{mag}.csv               # Summary statistics per magnification
groove_summary_mag{mag}_{detector}.csv    # Per-detector summary statistics
parallelism_summary_mag{mag}.csv          # Groove parallelism results
simulation_validation.csv                 # Pipeline validation on synthetic images

# R pipeline
Diffraction_grating.R                     # Main 1D analysis pipeline
bootstrap_fft_analysis.R                  # Bootstrap + 2D FFT analysis
gp_linearity.R                            # Gaussian process regression
simulate_validate.R                       # Synthetic image validation

# Python analyses (reproduce the published figures/numbers)
python/sem_fft_subpixel.py                # Sub-bin FFT refinement (16× zero-pad + parabolic)
python/sem_fft_figures.py                 # FFT invariance + 1D-FFT agreement figures
sem_detector_forest.py                    # Bootstrapped Bland-Altman forest plot (SED vs BED-S)
sem_ba_fig.py                             # Per-image Bland-Altman plot (1000x, n=56)
sem_magerr_fig.py                         # Magnification discrepancy + spacing figure (not used in the revised manuscript)
audit_numbers.py, audit_part2.py          # Full numeric audit of the manuscript vs CSVs
sem_gp.py                                 # Numerically stable GP fit (LRT, lengthscale, maps; replaces gp_linearity.R, whose chol() failed)
python/sem_pitch1d_physical.py            # 1D pitch + physical-unit linearity (Table 4)
python/sem_linearity_fig.py               # linearity_combined.pdf (stage coordinates)
python/method_comparison_pitch.py         # Cross-method comparison figure
wli/wli_batch.py                          # WLI batch analysis (Optiv DualZ 763)
wli/wli_example_2d.py                     # WLI example figure (3 panels)

# Additional analyses of the revised manuscript (revision 1)
revision/audit_metadata.py                     # Acquisition metadata of every image (field of view, pixel scale, stage position, time)
revision/analyse_profiles.R                    # 1D profile-minima analysis with the positions of all minima
revision/pipeline/03_smoothing_sensitivity.R   # Sensitivity of the 1D pitch to the smoothing parameter
revision/pipeline/04_referee_analyses.py       # Acceptance criterion, detector comparison, regressions on stage position and time
revision/pipeline/05_figures.py                # Revised figures (pitch vs stage position with insets, method comparison, acquisition schematic)
revision/pipeline/06_region_identity.py        # Test whether two series imaged the same area of the grating
revision/pipeline/10_nonperiodic_component.py  # Non-periodic part of the column-averaged profiles
revision/analysis/sem_metadata.csv             # Per-image acquisition metadata (summary: sem_metadata_summary.csv)
revision/analysis/profiles_per_image.csv       # Per-image minima positions, spacings and acceptance
revision/analysis/smoothing_variants_per_image.csv  # Per-image results for all smoothing settings
revision/analysis/referee_*.csv, referee_*.txt # Numerical results of the additional analyses
```

The scripts in `revision/` locate their input relative to their own position. `revision/pipeline/04_referee_analyses.py` runs on the files of this repository alone; the scripts that read the raw images (`audit_metadata.py`, `analyse_profiles.R`, `06_region_identity.py`, `10_nonperiodic_component.py` and parts of `05_figures.py`) expect the image folders `mag{mag}_SED-BED/{detector}/` at the top level, and `03_smoothing_sensitivity.R` uses the profile file written by `analyse_profiles.R`.

> **Deprecation note:** `linearity_helpers.R` and the `linearity_tiles_mag{mag}.csv` / `linearity_summary_overview.csv` outputs belong to the original tile-index pipeline. They are superseded by `python/sem_pitch1d_physical.py` + `pitch1d_physical.csv`, which evaluate linearity against the physical stage coordinates in µm/µm.

## Methods Summary

- **Instrument:** JEOL JSM-IT700HR FE-SEM, 8 kV, high vacuum, 1280 × 960 px, 0.0133 µm/px (field of view 17.07 µm × 12.80 µm for every image)
- **Grating:** Thorlabs GH13-06U, 600 lines/mm nominal (1.667 µm pitch; a design value without a stated uncertainty, not used as a reference)
- **1D method:** column-averaged profile, B-spline smoothing, minima detection
- **FFT method:** column-averaged intensity profile FFT with 16× zero-padding and sub-bin parabolic peak refinement (no window function; mandatory: integer bins quantise the pitch in 0.155 µm steps)
- **Variation across the scanned area:** per-image pitch regressed on stage coordinates (`$CM_STAGE_POSITION`) in physical units (µm/µm); within a series the stage y coordinate is collinear with the acquisition time
- **WLI cross-check:** Hexagon Optiv DualZ 763 multisensor CMM (1904 × 1200 px, 0.4745 µm/px), block-median 2D FFT with sub-bin refinement
- **Optical reference:** red diode laser (λ ≈ 635 nm), reflection geometry

## Key Findings

- Sub-bin-refined FFT pitch: **1.644 µm** for all series and both detectors (series medians 1.642–1.648 µm, SNR > 150), i.e. the instrument-reported scale was reproducible between the series within 0.4 %; consistency between different image magnifications is not tested by this dataset
- Naive integer-bin FFT: 1.707 µm (analytical artefact, not physical)
- 1D profile-minima pitch: 1.643–1.648 µm (bias vs FFT +0.0045 µm, 95% LoA [−0.062, +0.071] µm, n = 1432)
- WLI pitch: **1.6651 µm** (30 repeats: SD 0.07 nm; 36 positions: SD 0.13 nm)
- Optical diffraction: 1.65 ± 0.21 µm (k=3)
- Variation of the measured pitch across the scanned areas (stage coordinates): **≤ 0.9 %**, with a decrease along the acquisition (stage y) direction in all series that cannot be separated from a time-dependent change of the image scale
- SED–BED-S detector bias: ≤ 0.013 µm (bootstrapped Bland–Altman)

## License
CC-BY 4.0

## Citation
Kucharski, D. (2026). SEM Grating Pitch Benchmark Dataset [Data set]. Zenodo.
https://doi.org/10.5281/zenodo.20641960 (all versions; latest: v1.1,
https://doi.org/10.5281/zenodo.22045044)
