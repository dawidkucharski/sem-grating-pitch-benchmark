#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Method comparison of grating pitch: SEM 1D, SEM FFT, optical diffraction,
nominal, and the WLI campaign (Hexagon Optiv DualZ 763, 30 same-position
repeats + 36 positions). WLI error bar: empirical between-scan scatter
(0.00013 um, k=1); the annotation notes that absolute accuracy remains
bounded by the declared machine MPE (ISO 10360-7).
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "method_comparison_pitch.pdf")

methods = ["Nominal\n(600 lines/mm)",
           "SEM 1D\n(4 magnifications)",
           "SEM 2D FFT\n(sub-bin)",
           "Optical\ndiffraction",
           "WLI (CMM)\n36 positions"]

pitch = [1.667, 1.6455, 1.6442, 1.650, 1.66512]
yerr = [0.0, 0.0025, 0.0, 0.21, 0.00013]   # SEM 1D: range 1.643-1.648
colours = ["black", "tab:blue", "tab:blue", "tab:green", "tab:orange"]
markers = ["s", "o", "o", "D", "o"]

fig, ax = plt.subplots(figsize=(8, 4.6))
ax.axhline(1.667, color="grey", ls="--", lw=1, zorder=1)
ax.text(4.3, 1.6674, "nominal 1.667 µm", color="grey", fontsize=9, va="bottom")

for i, (m, p, e, c, mk) in enumerate(zip(methods, pitch, yerr, colours, markers)):
    ax.errorbar(i, p, yerr=e, fmt=mk, color=c, ms=9, capsize=6, lw=1.6,
                zorder=2)

ax.set_xticks(range(len(methods)))
ax.set_xticklabels(methods, fontsize=9)
ax.set_ylabel("pitch [µm]")
ax.set_ylim(1.40, 1.90)
ax.set_title("Grating pitch: method comparison\n"
             "WLI (Optiv DualZ 763): empirical scatter ±0.0001 µm"
             " between scans", fontsize=11)
ax.grid(axis="y", ls=":", alpha=0.5)
fig.subplots_adjust(right=0.68)
ax.annotate("WLI (rightmost point):\n"
            "error bar = ±1 SD between scans\n"
            "(0.00013 µm, n = 36 positions).\n"
            "Same-position repeatability:\n"
            "±0.00007 µm (n = 30 scans).\n"
            "No separate accuracy specification\n"
            "is declared for the WLI sensor.",
            xy=(4, 1.66512), xytext=(1.04, 0.72),
            textcoords="axes fraction",
            fontsize=8.5, color="0.30", ha="left", va="center",
            clip_on=False)
fig.tight_layout()
fig.savefig(OUT, dpi=200)
fig.savefig(OUT.replace(".pdf", ".png"), dpi=150)
print(f"Saved: {OUT}")
