#!/usr/bin/env python3
"""Non-periodic part of the column-averaged profile (Referee 2, comments 2-3).

For every image the profile is folded at its FFT period; the mean period
shape is subtracted and the local RMS of the residual (window: one period)
is evaluated.  The script reports
  frac      variance of the residual / variance of the profile
  peak      largest local RMS of the residual / median local RMS
  peak_px   position of the largest local RMS
and tests whether, in images with irregular inter-minimum spacings, the
largest residual coincides with the anomalous spacing.

Output: revision/analysis/referee_nonperiodic.csv and .txt
"""
import os

import numpy as np
import pandas as pd
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(HERE)
ROOT = os.path.dirname(REV)
OUT = os.path.join(REV, "analysis")
# Submitted result tables: R_output_dir/ in the working folder, top level in
# the public repository.
RES = os.path.join(ROOT, "R_output_dir")
if not os.path.isdir(RES):
    RES = ROOT
S = 0.01333594
BINS = 48
lines = []


def say(*a):
    text = " ".join(str(v) for v in a)
    print(text)
    lines.append(text)


p = pd.read_csv(os.path.join(OUT, "profiles_per_image.csv"))
p["file"] = p.image_path.str.split("/").str[-1]
fft = pd.read_csv(os.path.join(RES, "fft_subpixel_per_image.csv")).rename(
    columns={"mag": "group_label", "det": "detector", "pitch": "fft"})
p = p.merge(fft[["group_label", "detector", "file", "fft"]], on=["group_label", "detector", "file"])

rows = []
for r in p.itertuples():
    prof = np.asarray(Image.open(os.path.join(ROOT, r.image_path)).convert("L"), dtype=float).mean(axis=0)
    n = len(prof)
    period = r.fft / S
    phase = (np.arange(n) / period) % 1.0
    idx = np.minimum((phase * BINS).astype(int), BINS - 1)
    template = np.array([prof[idx == b].mean() for b in range(BINS)])
    centres = (np.arange(BINS) + 0.5) / BINS
    model = np.interp(phase, np.r_[centres[-1] - 1, centres, centres[0] + 1], np.r_[template[-1], template, template[0]])
    res = prof - model
    w = int(round(period))
    local = np.sqrt(np.convolve(res ** 2, np.ones(w) / w, mode="same"))
    core = slice(w // 2, n - w // 2)
    k = int(np.argmax(local[core])) + w // 2
    m = np.array([float(v) for v in str(r.minima_pixels_1based).split(";")]) - 1
    s = np.diff(m)
    bad = np.where((s / np.median(s) < 0.85) | (s / np.median(s) > 1.15))[0]
    hit, cover = np.nan, np.nan
    if len(bad):
        mask = np.zeros(n, bool)
        for j in bad:
            mask[max(0, int(m[j] - period / 2)):min(n, int(m[j + 1] + period / 2) + 1)] = True
        hit, cover = bool(mask[k]), float(mask[core].mean())
    rows.append(dict(image_path=r.image_path, group_label=r.group_label, detector=r.detector,
                     irregular=bool(len(bad)), frac=float(res.var() / prof.var()),
                     peak=float(local[core].max() / np.median(local[core])), peak_px=k, hit=hit, cover=cover))
d = pd.DataFrame(rows)
d.to_csv(os.path.join(OUT, "referee_nonperiodic.csv"), index=False)

say("== Non-periodic part of the column-averaged profile ==")
for flag, g in d.groupby("irregular"):
    say(f"irregular={flag}: n={len(g)}  variance fraction median {100 * g.frac.median():.2f} % "
        f"(quartiles {100 * g.frac.quantile(.25):.2f}-{100 * g.frac.quantile(.75):.2f} %)  "
        f"peak/median local RMS median {g.peak.median():.2f} (quartiles {g.peak.quantile(.25):.2f}-{g.peak.quantile(.75):.2f})")
say(d.groupby(["detector", "irregular"]).agg(n=("frac", "size"), frac_median_pct=("frac", lambda v: 100 * v.median()),
                                             peak_median=("peak", "median")).round(2).to_string())
irr = d[d.irregular]
say(f"irregular images in which the largest local residual lies at the anomalous spacing: {int(irr.hit.sum())} of {len(irr)} "
    f"(expected by chance: {irr.cover.sum():.1f}); by detector:",
    irr.groupby("detector").hit.agg(["sum", "size"]).to_dict("index"))
thr = d[~d.irregular].frac.quantile(.95)
say(f"irregular images above the 95th percentile of the regular images in the variance fraction ({100 * thr:.2f} %): "
    f"{int((irr.frac > thr).sum())} of {len(irr)}")
out = pd.read_csv(os.path.join(OUT, "referee_fig6_outliers.csv"))
o = out.merge(d, on="image_path")
ref = d[~d.irregular].frac.median()
say("Fig. 6 outliers (variance fraction in %, ratio to the median of the regular images, peak ratio):")
say(o.assign(frac_pct=100 * o.frac, ratio=o.frac / ref)[["mag", "det", "base", "pitch", "frac_pct", "ratio", "peak"]]
    .round(2).to_string(index=False))
with open(os.path.join(OUT, "referee_nonperiodic.txt"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
