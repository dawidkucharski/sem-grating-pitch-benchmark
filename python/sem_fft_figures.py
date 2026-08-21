#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Regenerate fft_invariance.pdf and fft_1d_agreement.pdf from the sub-bin
FFT re-analysis and the per-image 1D values (linearity_tiles CSVs).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/dawid/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/" \
       "SEM_measurements/Calibration"
R = os.path.join(BASE, "R_output_dir")

fft = pd.read_csv(os.path.join(R, "fft_subpixel_per_image.csv"))
fft["base"] = fft["file"].str.replace(".png", "", regex=False)

tiles = []
for mag in [1000, 750, 500, 400]:
    t = pd.read_csv(os.path.join(R, f"linearity_tiles_mag{mag}.csv"))
    t["base"] = t["image_path"].str.split("/").str[-1].str.replace(
        ".png", "", regex=False)
    tiles.append(t)
lin = pd.concat(tiles)
lin["det"] = lin["image_path"].str.contains("BED-S").map(
    {True: "BED-S", False: "SED"})

df = fft.merge(lin[["base", "local_pitch_um", "det"]].drop_duplicates("base"),
               on="base", how="inner")
print(f"merged per-image rows: {len(df)}")

# ---- Figure 1: FFT invariance ----
fig, ax = plt.subplots(figsize=(6.5, 4.2))
mags = [1000, 750, 500, 400]
med, sd = [], []
for mag in mags:
    p = df[df.mag == mag].pitch
    med.append(p.median())
    sd.append(p.std(ddof=1))
    ax.scatter(np.full(len(p), mag), p, s=6, alpha=0.25, color="tab:blue")
ax.errorbar(mags, med, yerr=sd, fmt="o", color="tab:red", ms=8, capsize=5,
            lw=1.6, label="median ± 1 SD")
ax.axhline(1.667, color="grey", ls="--", lw=1.2)
ax.text(1005, 1.6673, "nominal 1.667 µm", color="grey", fontsize=8.5)
ax.set_xlim(380, 1020)
ax.set_ylim(1.60, 1.70)
ax.set_xlabel("magnification")
ax.set_ylabel("pitch (sub-bin 2D FFT) [µm]")
ax.set_title("Sub-bin FFT pitch across magnifications (IFR)\n"
             f"per-mag medians {med[0]:.4f}, {med[1]:.4f}, "
             f"{med[2]:.4f}, {med[3]:.4f} µm", fontsize=10)
ax.legend(loc="lower right")
ax.grid(ls=":", alpha=0.5)
fig.tight_layout()
fig.savefig(os.path.join(BASE, "fft_invariance.pdf"), dpi=200)
print("Saved fft_invariance.pdf")

# ---- Figure 2: 1D vs FFT agreement ----
d = df.local_pitch_um - df.pitch
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
ax = axes[0]
for mag, c in zip(mags, ["tab:blue", "tab:orange", "tab:green", "tab:red"]):
    m = df.mag == mag
    ax.scatter(df.loc[m, "local_pitch_um"], df.loc[m, "pitch"], s=8,
               alpha=0.35, color=c, label=f"{mag}×")
lims = [1.60, 1.68]
ax.plot(lims, lims, "k--", lw=1)
ax.set_xlim(lims); ax.set_ylim(lims)
ax.set_xlabel("1D profile-minima pitch [µm]")
ax.set_ylabel("sub-bin 2D FFT pitch [µm]")
ax.set_title("Per-image method agreement")
ax.legend()
ax.grid(ls=":", alpha=0.5)
ax = axes[1]
mn = (df.local_pitch_um + df.pitch) / 2
ax.scatter(mn, d, s=8, alpha=0.35, color="tab:blue")
md = d.mean(); sd_d = d.std(ddof=1)
ax.axhline(md, color="red", lw=1.4,
           label=f"bias = {md:+.4f} µm")
ax.axhline(md + 1.96 * sd_d, color="grey", ls="--", lw=1)
ax.axhline(md - 1.96 * sd_d, color="grey", ls="--", lw=1,
           label=f"95% LoA = [{md-1.96*sd_d:+.4f}, {md+1.96*sd_d:+.4f}]")
ax.set_xlabel("mean pitch [µm]")
ax.set_ylabel("1D − FFT [µm]")
ax.set_title("Bland–Altman (1D vs sub-bin FFT)")
ax.legend(fontsize=8)
ax.grid(ls=":", alpha=0.5)
fig.tight_layout()
fig.savefig(os.path.join(BASE, "fft_1d_agreement.pdf"), dpi=200)
print("Saved fft_1d_agreement.pdf")
print(f"n={len(d)} | bias={md:+.5f} µm | sd={sd_d:.5f} | "
      f"LoA=[{md-1.96*sd_d:+.5f}, {md+1.96*sd_d:+.5f}]")
