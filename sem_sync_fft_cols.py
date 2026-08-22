#!/usr/bin/env python3
"""Sync the FFT columns of the per-image pitch CSVs with the canonical
sub-bin FFT results (fft_subpixel_per_image.csv). The old columns carried
integer-bin values (1.707 um)."""
import os
import pandas as pd
import numpy as np

BASE = "/Users/dawid/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/SEM_measurements/Calibration"
S = 0.01333594
fft = pd.read_csv(f"{BASE}/R_output_dir/fft_subpixel_per_image.csv")
fft["base"] = fft.file.str.replace(r"\.(png|tif|txt)$", "", regex=True)
lookup = {(r.mag, r.det, r.base): (r.pitch, r.snr) for _, r in fft.iterrows()}

for d in [f"{BASE}/R_output_dir", f"{BASE}/sem_grating_benchmark"]:
    for m in [400, 500, 750, 1000]:
        for det in ["SED", "BED-S"]:
            p = f"{d}/per_image_pitch_mag{m}_{det}.csv"
            df = pd.read_csv(p)
            df["base"] = df.image_path.str.extract(r"(Grp_\d+_[A-Za-z0-9-]+_X\d+_Y\d+)")[0]
            n_old = df.pitch_fft_um.notna().sum()
            df["pitch_fft_um"] = [lookup.get((m, det, b), (np.nan, np.nan))[0]
                                  for b in df.base]
            df["fft_snr"] = [lookup.get((m, det, b), (np.nan, np.nan))[1]
                             for b in df.base]
            df["fft_peak_freq"] = S / df.pitch_fft_um
            df["fft_quality"] = "sub-bin"
            df = df.drop(columns=["base"])
            df.to_csv(p, index=False)
            print(f"updated {p}  (old fft rows: {n_old}, new: {df.pitch_fft_um.notna().sum()})")
print("done")
