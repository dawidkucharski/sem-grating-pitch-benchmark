#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Regenerate elsarticle/linearity_combined.pdf with physical stage coordinates:
local groove spacing vs stage x [um] for both detectors, 4 magnifications.
Slopes b_x in um/um annotated.
"""
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = ("/Users/dawid/Library/Mobile Documents/com~apple~CloudDocs/"
        "Dokumenty/SEM_measurements/Calibration")
df = pd.read_csv(f"{BASE}/R_output_dir/pitch1d_physical.csv").dropna(
    subset=["pitch"])

fig, axes = plt.subplots(2, 2, figsize=(8.5, 7.2), sharex=False)
for ax, mag in zip(axes.ravel(), [400, 500, 750, 1000]):
    g = df[df.mag == mag]
    for det, c, m in [("SED", "#1f77b4", "o"), ("BED-S", "#d62728", "s")]:
        gg = g[g.det == det]
        ax.plot(gg.stage_x_um, gg.pitch, m, ms=4, mfc="none",
                mec=c, mew=0.7, alpha=0.75, label=det)
    bx = stats.linregress(g.stage_x_um, g.pitch)
    xx = np.linspace(g.stage_x_um.min(), g.stage_x_um.max(), 50)
    ax.plot(xx, bx.intercept + bx.slope * xx, "k-", lw=1.2)
    dbar = g.pitch.mean()
    ax.axhline(dbar, color="0.45", ls="--", lw=0.8)
    xr = g.stage_x_um.max() - g.stage_x_um.min()
    ax.set_title(f"${mag}\\times$   "
                 f"$b_x = {bx.slope:.1e}$ µm/µm ($p = {bx.pvalue:.2f}$)", fontsize=10)
    ax.text(0.03, 0.95,
            f"max. rel. change $= {100*abs(bx.slope)*xr/dbar:.2f}$\\%",
            transform=ax.transAxes, va="top", fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.25", fc="0.95", ec="0.7", lw=0.4))
    ax.set_xlabel("Stage position $x$ [µm]")
    ax.set_ylabel("Local groove spacing [µm]")
    ax.tick_params(labelsize=8)
    if mag == 400:
        ax.legend(loc="lower right", fontsize=8, frameon=False)

fig.tight_layout()
fig.savefig(f"{BASE}/elsarticle/linearity_combined.pdf")
print("saved", f"{BASE}/elsarticle/linearity_combined.pdf")

# print table values for manuscript
print(f"{'mag':>5} {'n':>4} {'dx':>6} {'dy':>6} {'bx':>9} {'p_x':>6} "
      f"{'by':>9} {'p_y':>6} {'maxrel%':>8}")
for mag, g in df.groupby("mag"):
    xr = g.stage_x_um.max() - g.stage_x_um.min()
    yr = g.stage_y_um.max() - g.stage_y_um.min()
    bx = stats.linregress(g.stage_x_um, g.pitch)
    by = stats.linregress(g.stage_y_um, g.pitch)
    dbar = g.pitch.mean()
    mr = max(100*abs(bx.slope)*xr/dbar, 100*abs(by.slope)*yr/dbar)
    print(f"{mag:>5} {len(g):>4} {xr:>6.0f} {yr:>6.0f} {bx.slope:>9.1e} "
          f"{bx.pvalue:>6.3f} {by.slope:>9.1e} {by.pvalue:>6.3f} {mr:>8.2f}")
