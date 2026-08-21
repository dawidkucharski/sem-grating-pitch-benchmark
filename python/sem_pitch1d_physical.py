#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-image 1D profile-minima pitch for both detectors (same procedure as the
R pipeline: column profile, B-spline, minima, regular-period filter), joined
with stage coordinates -> linearity slopes in physical units (um/um).
"""
import os
import glob
import numpy as np
import pandas as pd
from PIL import Image
from scipy.interpolate import UnivariateSpline
from scipy import stats

BASE = "/Users/dawid/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/" \
       "SEM_measurements/Calibration"
S = 0.01333594


def pitch_1d(arr):
    prof = arr.mean(axis=0)                       # mean over y, along x
    n = len(prof)
    x = np.arange(n)
    spl = UnivariateSpline(x, prof,
                           s=0.5 * np.sum((prof - prof.mean()) ** 2))
    sp = spl(x)
    d = np.diff(sp)
    minima = np.where((d[:-1] < 0) & (d[1:] > 0))[0] + 1
    if minima.size < 2:
        return np.nan
    dist = np.diff(minima) * S
    med = np.median(dist)
    keep = (dist >= 0.7 * med) & (dist <= 1.3 * med)
    if keep.sum() == 0:
        return np.nan
    return float(dist[keep].mean())


def stage_pos(path_txt):
    try:
        with open(path_txt, errors="ignore") as fh:
            for line in fh:
                if line.startswith("$CM_STAGE_POSITION"):
                    v = line.split()[1:]
                    return float(v[0]), float(v[1])
    except Exception:
        return None
    return None


rows = []
for mag in [1000, 750, 500, 400]:
    folder = os.path.join(BASE, f"mag{mag}_SED-BED")
    for det in ["SED", "BED-S"]:
        pngs = sorted(glob.glob(os.path.join(folder, det, "*.png")))
        coords = {}
        for f in glob.glob(os.path.join(folder, det, "*.txt")):
            coords[os.path.basename(f).replace(".txt", "")] = stage_pos(f)
        for p in pngs:
            b = os.path.basename(p).replace(".png", "")
            sp = coords.get(b)
            if sp is None:
                continue
            a = np.asarray(Image.open(p).convert("L"), dtype=np.float64)
            rows.append(dict(mag=mag, det=det, base=b, pitch=pitch_1d(a),
                             stage_x_um=sp[0] * 1000, stage_y_um=sp[1] * 1000))
        print(f"mag{mag} {det}: {len(pngs)} done")

df = pd.DataFrame(rows)
df.to_csv(os.path.join(BASE, "R_output_dir", "pitch1d_physical.csv"),
          index=False)

print(f"\n{'mag':>5} {'det':>6} {'n':>4} {'median':>8} "
      f"{'bx[um/um]':>10} {'p_x':>8} {'by[um/um]':>10} {'p_y':>8} "
      f"{'max_rel%':>9}")
for (mag, det), g in df.groupby(["mag", "det"]):
    xr = g.stage_x_um.max() - g.stage_x_um.min()
    yr = g.stage_y_um.max() - g.stage_y_um.min()
    bx = stats.linregress(g.stage_x_um, g.pitch)
    by = stats.linregress(g.stage_y_um, g.pitch)
    dbar = g.pitch.mean()
    max_rel = max(abs(bx.slope * xr) / dbar, abs(by.slope * yr) / dbar) * 100
    print(f"{mag:>5} {det:>6} {len(g):>4} {g.pitch.median():>8.4f} "
          f"{bx.slope:>10.2e} {bx.pvalue:>8.3f} {by.slope:>10.2e} "
          f"{by.pvalue:>8.3f} {max_rel:>9.2f}")
