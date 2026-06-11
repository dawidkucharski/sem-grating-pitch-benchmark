# SEM Grating Pitch Benchmark Dataset

## Overview
This dataset accompanies the manuscript *"Grating-based SEM magnification assessment using a commercial holographic diffraction grating: detector comparison, field linearity, and uncertainty budgeting"* submitted to Measurement Science and Technology.

The dataset comprises 1500 SEM images of a Thorlabs GH13-06U reflective holographic diffraction grating (600 lines/mm nominal), acquired on a JEOL JSM-IT700HR field-emission SEM at four magnifications (400×, 500×, 750×, 1000×) using both secondary-electron (SED) and backscattered-electron (BED-S) detectors.

## Contents

```
R_output_dir/
  per_image_pitch_mag{mag}_{detector}.csv   # Per-image 1D and 2D FFT pitch values
  bootstrap_summary.csv                     # Bootstrap BCa confidence intervals
  bland_altman_bootstrap_summary.csv        # Bootstrapped Bland-Altman statistics
  linearity_tiles_mag{mag}.csv              # Tile-level linearity data
  linearity_summary_overview.csv            # Linearity regression summary
  gp_linearity_summary.csv                  # Gaussian process regression summary
  groove_summary_mag{mag}.csv               # Summary statistics per magnification
  parallelism_summary_mag{mag}.csv          # Groove parallelism results
  simulation_validation.csv                 # Pipeline validation on synthetic images

R scripts:
  Diffraction_grating.R                     # Main analysis pipeline
  linearity_helpers.R                       # Linearity assessment
  bootstrap_fft_analysis.R                  # Bootstrap + 2D FFT analysis
  gp_linearity.R                            # Gaussian process regression
  simulate_validate.R                       # Synthetic image validation
```

## Methods Summary

- **Instrument:** JEOL JSM-IT700HR FE-SEM, 8 kV, high vacuum
- **Grating:** Thorlabs GH13-06U, 600 lines/mm nominal (1.667 µm pitch)
- **Image size:** 1280 × 960 pixels, 0.0133 µm/px
- **1D method:** B-spline smoothed profile minima detection
- **2D FFT method:** Hann-windowed 2D FFT, fundamental peak at kx=0
- **Optical reference:** Red diode laser (λ ≈ 635 nm), grating equation

## Key Findings

- 2D FFT pitch: 1.707 µm (invariant across magnifications and detectors)
- Optical diffraction pitch: 1.70 ± 0.11 µm (k=3)
- Agreement between optical and FFT: 0.2%
- Both methods indicate the actual grating pitch is ~2.4% above nominal

## License
CC-BY 4.0

## Citation
Kucharski, D. (2026). SEM Grating Pitch Benchmark Dataset [Data set]. Zenodo. https://doi.org/XXXXX
