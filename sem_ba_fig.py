#!/usr/bin/env python3
"""Regenerate elsarticle/bland_altman_1000x.pdf from the canonical
per-image pitch data (tile pairing, n=56 at 1000x)."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

BASE = "/Users/dawid/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/SEM_measurements/Calibration"
GB = os.path.join(BASE, "sem_grating_benchmark")

s = pd.read_csv(f"{GB}/per_image_pitch_mag1000_SED.csv").dropna(subset=["pitch_1d_um"])
b = pd.read_csv(f"{GB}/per_image_pitch_mag1000_BED-S.csv").dropna(subset=["pitch_1d_um"])
m = s.merge(b, on=["tile_x", "tile_y"], suffixes=("_S", "_B"))
diff = (m.pitch_1d_um_S - m.pitch_1d_um_B).values
meanv = (m.pitch_1d_um_S + m.pitch_1d_um_B).values / 2
bias = diff.mean()
sd = diff.std(ddof=1)
loa_lo, loa_hi = bias - 1.96 * sd, bias + 1.96 * sd
p = stats.ttest_rel(m.pitch_1d_um_S, m.pitch_1d_um_B).pvalue

fig, ax = plt.subplots(figsize=(5.8, 4.6))
ax.scatter(meanv, diff, s=26, color="#3366cc", alpha=0.55, edgecolors="none")
ax.axhline(bias, color="crimson", lw=1.2)
ax.axhline(loa_lo, color="0.35", ls="--", lw=1.0)
ax.axhline(loa_hi, color="0.35", ls="--", lw=1.0)

# axes limits with margins so all lines and labels are visible
xlo = meanv.min() - 0.004
xhi = meanv.max() + 0.004
ylo = min(diff.min(), loa_lo) - 0.035
yhi = max(diff.max(), loa_hi) + 0.035
ax.set_xlim(xlo, xhi)
ax.set_ylim(ylo, yhi)

# line labels inside, on the left edge of the axes
ax.text(xlo + 0.0015, bias + 0.008, f"bias = {bias:+.4f} \u00b5m",
        ha="left", va="bottom", color="crimson", fontsize=8.5)
ax.text(xlo + 0.0015, loa_hi + 0.004, f"+1.96 SD = {loa_hi:+.4f} \u00b5m",
        ha="left", va="bottom", color="0.25", fontsize=8)
ax.text(xlo + 0.0015, loa_lo - 0.004, f"\u22121.96 SD = {loa_lo:+.4f} \u00b5m",
        ha="left", va="top", color="0.25", fontsize=8)

# summary block inside the top-right corner of the axes
ax.text(0.985, 0.965,
        f"$n$ = {len(diff)} pairs\npaired $t$-test $p$ = {p:.2f}",
        transform=ax.transAxes, ha="right", va="top", fontsize=8,
        color="0.2",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.8", lw=0.5))

ax.set_xlabel("Mean groove spacing (SED + BED-S)/2 [\u00b5m]", fontsize=9.5)
ax.set_ylabel("Difference SED \u2212 BED-S [\u00b5m]", fontsize=9.5)
ax.set_title("Bland\u2013Altman: SED vs BED-S at 1000\u00d7", fontsize=11)
ax.tick_params(labelsize=8.5)
fig.tight_layout()
fig.savefig(f"{BASE}/elsarticle/bland_altman_1000x.pdf")
print(f"saved bland_altman_1000x.pdf  n={len(diff)} bias={bias:.4f} sd={sd:.4f} "
      f"LoA=[{loa_lo:.4f}, {loa_hi:.4f}] p={p:.3f}")
