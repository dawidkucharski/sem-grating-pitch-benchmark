#!/usr/bin/env python3
"""Forest plot: bootstrapped Bland-Altman SED vs BED-S (10 000 replicates, 95% BCa).

Reads sem_grating_benchmark/bland_altman_bootstrap_summary.csv and writes
detector_forest.pdf. Title is split over two lines so it never clips.
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/dawid/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/SEM_measurements/Calibration"

rows = []
with open(f"{BASE}/sem_grating_benchmark/bland_altman_bootstrap_summary.csv") as f:
    for r in csv.DictReader(f):
        rows.append(r)

mags = [int(r["magnification"]) for r in rows]
bias = np.array([float(r["bias_um"]) for r in rows])
b_lo = np.array([float(r["bias_ci_low"]) for r in rows])
b_hi = np.array([float(r["bias_ci_upp"]) for r in rows])
loa_l = np.array([float(r["loa_lower_um"]) for r in rows])
loa_u = np.array([float(r["loa_upper_um"]) for r in rows])
l_lo = np.array([float(r["loa_lower_ci_low"]) for r in rows])
l_hi = np.array([float(r["loa_lower_ci_upp"]) for r in rows])
u_lo = np.array([float(r["loa_upper_ci_low"]) for r in rows])
u_hi = np.array([float(r["loa_upper_ci_upp"]) for r in rows])
n_pairs = [int(r["n_pairs"]) for r in rows]

fig, ax = plt.subplots(figsize=(7.0, 3.4))
y = np.arange(len(mags))[::-1]  # 1000x on top

# LoA intervals (grey) with their BCa CI whiskers
for yi, (ll, lu, llo, lhi, ulo, uhi) in enumerate(zip(loa_l, loa_u, l_lo, l_hi, u_lo, u_hi)):
    yy = y[yi]
    ax.plot([ll, lu], [yy, yy], color="0.6", lw=1.4, zorder=2)
    ax.plot([llo, lhi], [yy, yy], color="0.6", lw=0.8, marker="|",
            ms=4, mew=0.8, mec="0.6", zorder=2)
    ax.plot([ulo, uhi], [yy, yy], color="0.6", lw=0.8, marker="|",
            ms=4, mew=0.8, mec="0.6", zorder=2)

# Bias (dark red) with BCa CI
ax.errorbar(bias, y, xerr=[bias - b_lo, b_hi - bias], fmt="D", ms=6,
            color="#8B0000", ecolor="#8B0000", elinewidth=1.4,
            capsize=3, zorder=3)

ax.axvline(0, color="k", lw=0.8, ls="--", zorder=1)
ax.set_yticks(y)
ax.set_yticklabels([f"{m}$\\times$" for m in mags])
ax.set_xlabel("SED $-$ BED-S pitch difference [\u00b5m]", fontsize=10)
ax.set_title("SED $-$ BED-S: bootstrapped Bland\u2013Altman\n"
             "(10 000 replicates, 95% BCa)", fontsize=11)
ax.tick_params(labelsize=9)

# right-hand annotations
xmax = max(np.max(u_hi), np.max(loa_u)) * 1.18
ax.set_xlim(min(np.min(l_lo), np.min(b_lo)) * 1.25, xmax)
for yi, (b, n) in enumerate(zip(bias, n_pairs)):
    ax.text(xmax * 0.985, y[yi] + 0.22, f"Bias: {b:+.4f}", ha="right",
            va="bottom", fontsize=8, color="#8B0000")
    ax.text(xmax * 0.985, y[yi] - 0.22, f"$n$ = {n}", ha="right",
            va="top", fontsize=8, color="0.3")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(f"{BASE}/detector_forest.pdf")
print("saved detector_forest.pdf")
