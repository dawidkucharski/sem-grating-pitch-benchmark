# SEM Grating Pitch Benchmark Dataset

## Overview
This dataset accompanies the manuscript *"Self-consistent SEM magnification linearity assessment using a holographic grating: 2D Fourier analysis, detector comparison, and open benchmark dataset"*.

The dataset comprises 1500 SEM images in total (750 per detector) of a Thorlabs GH13-06U reflective holographic diffraction grating (600 lines/mm nominal), acquired on a JEOL JSM-IT700HR field-emission SEM at four magnifications (400×, 500×, 750×, 1000×) using both secondary-electron (SED) and backscattered-electron (BED-S) detectors: 400×: 19×19 grid (361 SED + 361 BED-S), 500×: 15×15 (225 + 225), 750×: 10×10 (100 + 100), 1000×: 8×8 (64 + 64).

## Contents

```
# Per-image and aggregate results
per_image_pitch_mag{mag}_{detector}.csv   # Per-image 1D + 2D FFT pitch values
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
python/sem_pitch1d_physical.py            # 1D pitch + physical-unit linearity (Table 4)
python/sem_linearity_fig.py               # linearity_combined.pdf (stage coordinates)
python/method_comparison_pitch.py         # Cross-method comparison figure
wli/wli_batch.py                          # WLI batch analysis (Optiv DualZ 763)
wli/wli_example_2d.py                     # WLI example figure (3 panels)
```

> **Deprecation note:** `linearity_helpers.R` and the `linearity_tiles_mag{mag}.csv` / `linearity_summary_overview.csv` outputs belong to the original tile-index pipeline. They are superseded by `python/sem_pitch1d_physical.py` + `pitch1d_physical.csv`, which evaluate linearity against the physical stage coordinates in µm/µm.

## Methods Summary

- **Instrument:** JEOL JSM-IT700HR FE-SEM, 8 kV, high vacuum, 1280 × 960 px, 0.0133 µm/px
- **Grating:** Thorlabs GH13-06U, 600 lines/mm nominal (1.667 µm pitch)
- **1D method:** column-averaged profile, B-spline smoothing, minima detection
- **2D FFT method:** Hann-windowed 2D FFT with 16× zero-padding and sub-bin parabolic peak refinement (mandatory: integer bins quantise the pitch in 0.155 µm steps)
- **Linearity:** per-image pitch regressed on stage coordinates (`$CM_STAGE_POSITION`) in physical units (µm/µm)
- **WLI cross-check:** Hexagon Optiv DualZ 763 multisensor CMM (1904 × 1200 px, 0.4745 µm/px), block-median 2D FFT with sub-bin refinement
- **Optical reference:** red diode laser (λ ≈ 635 nm), reflection geometry

## Key Findings

- Sub-bin-refined 2D FFT pitch: **1.644 µm**, invariant across magnifications and detectors (per-magnification medians 1.642–1.648 µm, SNR > 150)
- Naive integer-bin FFT: 1.707 µm (analytical artefact, not physical)
- 1D profile-minima pitch: 1.643–1.648 µm (agrees with FFT to +0.002 µm)
- WLI pitch: **1.6651 µm** (30 repeats: SD 0.07 nm; 36 positions: SD 0.13 nm)
- Optical diffraction: 1.65 ± 0.21 µm (k=3)
- Field linearity (stage coordinates): **≤ 0.5%** across every scanned field (≤ 0.9% for a single detector)
- SED–BED-S detector bias: ≤ 0.013 µm (bootstrapped Bland–Altman)

## License
CC-BY 4.0

## Citation
Kucharski, D. (2026). SEM Grating Pitch Benchmark Dataset [Data set]. Zenodo.
https://doi.org/10.5281/zenodo.20641960 (all versions; latest: v1.1,
https://doi.org/10.5281/zenodo.22045044)
