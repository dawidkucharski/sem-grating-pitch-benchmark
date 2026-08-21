#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sub-bin re-analysis of the SEM FFT pitch on the existing image dataset
(fixes review point M1: bin-quantised FFT pitch).

For each image (1280x960 px): column-averaged intensity profile along x,
zero-padded FFT (16x), parabolic sub-bin refinement of the peak in the
band corresponding to 1.2-2.8 um pitch. Outputs per-image pitch values
and per magnification/detector summaries.
"""
import os
import glob
import numpy as np
import pandas as pd
from PIL import Image

BASE = "/Users/dawid/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/" \
       "SEM_measurements/Calibration"
S = 0.01333594            # um/px (17.07 um / 1280 px)
PAD = 16
K0 = 10                   # nominal fundamental bin (1.707 um)
MAGS = [1000, 750, 500, 400]
DETS = ["SED", "BED-S"]


def fft_pitch(img_path):
    im = Image.open(img_path).convert("L")
    a = np.asarray(im, dtype=np.float64)
    prof = a.mean(axis=0)                       # grooves vertical -> x profile
    n = len(prof)
    npad = n * PAD
    F = np.fft.rfft(prof - prof.mean(), n=npad)
    P = np.abs(F)
    k_lo = int(S * npad / 2.8)      # pitch band 1.2-2.8 um -> k = S*npad/pitch
    k_hi = int(S * npad / 1.2)
    k_lo = max(1, min(k_lo, npad // 2 - 1))
    k_hi = min(k_hi, npad // 2 - 1)
    if k_hi <= k_lo:
        return np.nan, np.nan, np.nan
    kmax = k_lo + int(np.argmax(P[k_lo:k_hi + 1]))
    if kmax <= 0 or kmax >= npad // 2:
        return np.nan, np.nan, np.nan
    # parabolic refinement on magnitude
    p0, p1, p2 = P[kmax - 1], P[kmax], P[kmax + 1]
    den = p0 - 2 * p1 + p2
    delta = 0.0 if den == 0 else np.clip(0.5 * (p0 - p2) / den, -0.5, 0.5)
    k_est = kmax + delta
    pitch = S * npad / k_est
    snr = p1 / np.median(P[k_lo:k_hi + 1]) if np.median(P[k_lo:k_hi + 1]) > 0 \
        else np.nan
    quant_pitch = S * npad / kmax            # integer-bin (old method)
    return pitch, quant_pitch, snr


def main():
    rows = []
    for mag in MAGS:
        folder = os.path.join(BASE, f"mag{mag}_SED-BED")
        for det in DETS:
            files = sorted(glob.glob(os.path.join(folder, det, "*.png")))
            print(f"mag{mag} {det}: {len(files)} images")
            for f in files:
                p, q, snr = fft_pitch(f)
                rows.append(dict(mag=mag, det=det, file=os.path.basename(f),
                                 pitch=p, quant=q, snr=snr))
    df = pd.DataFrame(rows)
    out = os.path.join(BASE, "R_output_dir", "fft_subpixel_per_image.csv")
    df.to_csv(out, index=False)
    print(f"\nSaved: {out}\n")

    for (mag, det), g in df.groupby(["mag", "det"]):
        p = g.pitch
        print(f"mag{mag} {det}: n={len(g)} | sub-bin pitch median "
              f"{p.median():.5f} um, sd {p.std(ddof=1):.6f}, range "
              f"{p.min():.5f}-{p.max():.5f} | integer-bin pitch "
              f"{g.quant.median():.5f}")
    p = df.pitch
    print(f"\nALL (n={len(p)}): median {p.median():.5f} um | sd {p.std(ddof=1):.6f} um "
          f"({p.std(ddof=1)/p.median()*100:.3f}% rel) | range {p.min():.5f}-{p.max():.5f}")
    # quantised spread for comparison
    q = df.quant
    print(f"quantised: unique values {np.unique(q)}")
    # relative spread of refined pitch across mags
    for mag in MAGS:
        g = df[df.mag == mag].pitch
        print(f"mag{mag}: median {g.median():.5f} sd {g.std(ddof=1):.6f}")


if __name__ == "__main__":
    main()
