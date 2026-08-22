#!/usr/bin/env python3
"""Regenerate summary_magnification_errors.pdf:
(a) apparent magnification discrepancy eps_M vs magnification,
(b) median groove spacing with k=3 expanded uncertainty bars and
reference lines at nominal 1.667 and sub-bin FFT 1.644.
All labels inside the axes."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/dawid/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/SEM_measurements/Calibration"
GB = os.path.join(BASE, "sem_grating_benchmark")
MAGS = [400, 500, 750, 1000]
NOM, FFT, U = 1.667, 1.644, 0.09

rows = []
for m in MAGS:
    d = []
    for det in ["SED", "BED-S"]:
        df = pd.read_csv(f"{GB}/per_image_pitch_mag{m}_{det}.csv")
        d.append(df.pitch_1d_um.dropna())
    d = pd.concat(d)
    rows.append((m, d.median(), (d.median() - NOM) / NOM * 100))
df = pd.DataFrame(rows, columns=["mag", "d", "eps"])

fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.6))

# (a) apparent magnification discrepancy
ax = axes[0]
ax.plot(df.mag, df.eps, "o-", color="#1f77b4", ms=6, lw=1.4)
for _, r in df.iterrows():
    ax.annotate(f"{r.eps:.2f}%", (r.mag, r.eps), textcoords="offset points",
                xytext=(0, 8), ha="center", fontsize=8.5, color="#1f77b4")
ax.set_xlabel("Magnification", fontsize=9.5)
ax.set_ylabel("$\\varepsilon_M$ [%]", fontsize=9.5)
ax.set_title("(a) Apparent magnification discrepancy", fontsize=9.5)
ax.set_xticks(MAGS)
ax.set_xlim(370, 1030)
ax.set_ylim(-1.60, -0.95)
ax.tick_params(labelsize=8.5)

# (b) median spacing with U bars and references (labels inside the axes)
ax = axes[1]
ax.errorbar(df.mag, df.d, yerr=U, fmt="o-", color="#1f77b4", ms=6,
            capsize=4, lw=1.4, label="median spacing $\\pm U$ ($k=3$)")
ax.axhline(NOM, color="0.35", ls="--", lw=1.0)
ax.axhline(FFT, color="crimson", ls="--", lw=1.0)
ax.text(395, NOM + 0.004, "nominal 1.667 \u00b5m", color="0.35",
        fontsize=8, ha="left", va="bottom")
ax.text(395, FFT - 0.006, "sub-bin FFT 1.644 \u00b5m", color="crimson",
        fontsize=8, ha="left", va="top")
ax.set_xlabel("Magnification", fontsize=9.5)
ax.set_ylabel("Groove spacing [\u00b5m]", fontsize=9.5)
ax.set_title("(b) Median groove spacing", fontsize=9.5)
ax.set_xticks(MAGS)
ax.set_xlim(370, 1030)
ax.set_ylim(1.50, 1.80)
ax.tick_params(labelsize=8.5)
ax.legend(fontsize=7.5, loc="lower right", frameon=False)

fig.tight_layout()
fig.savefig(f"{BASE}/summary_magnification_errors.pdf")
print("saved summary_magnification_errors.pdf")
print(df.to_string(index=False))
