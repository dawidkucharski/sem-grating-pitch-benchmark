#!/usr/bin/env python3
"""Analyses requested by the referees of MST-139559 (no change of the
submitted estimators, pixel scale or acceptance rule).

All per-image values are those of the submitted analysis:
  1D   revision/analysis/profiles_per_image.csv (reproduces R_output_dir)
  FFT  R_output_dir/fft_subpixel_per_image.csv
  Fig. 6 data  R_output_dir/pitch1d_physical.csv
  Fig. 5 data  R_output_dir/linearity_tiles_mag*.csv

Output: revision/analysis/referee_analyses.txt and referee_*.csv
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

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
MAGS = [400, 500, 750, 1000]
lines = []


def say(*a):
    text = " ".join(str(v) for v in a)
    print(text)
    lines.append(text)


p = pd.read_csv(os.path.join(OUT, "profiles_per_image.csv"))
meta = pd.read_csv(os.path.join(OUT, "sem_metadata.csv"))[["image_path", "timestamp"]]
p["time"] = pd.to_datetime(p.timestamp)
p["t_min"] = p.groupby("group_label").time.transform(lambda t: (t - t.min()).dt.total_seconds() / 60)
p["file"] = p.image_path.str.split("/").str[-1]
fft = pd.read_csv(os.path.join(RES, "fft_subpixel_per_image.csv"))
p = p.merge(fft.rename(columns={"mag": "group_label", "det": "detector", "pitch": "pitch_fft_um"})[
    ["group_label", "detector", "file", "pitch_fft_um"]], on=["group_label", "detector", "file"])
assert len(p) == 1500

# ---------------------------------------------------------------- A
say("== A. Inter-minimum spacings relative to the image mean (all 1500 images) ==")
rows = []
for r in p.itertuples():
    m = np.array([float(v) for v in str(r.minima_pixels_1based).split(";")])
    s = np.diff(m)
    for k, v in enumerate(s):
        rows.append((r.image_path, r.group_label, r.detector, k, len(s), 0.5 * (m[k] + m[k + 1]),
                     v, v / s.mean(), v / np.median(s)))
L = pd.DataFrame(rows, columns="image_path group detector k n centre_px px r_mean r_med".split())
bins = [0, 0.3, 0.6, 0.8, 0.9, 1.1, 1.2, 1.4, 1.6, 2.0, 3.0]
tab = pd.crosstab(pd.cut(L.r_mean, bins), L.detector)
tab["all"] = tab.sum(axis=1)
say(tab.to_string())
say("spacings:", len(L), "| within 0.9-1.1 of image mean:", int(((L.r_mean > .9) & (L.r_mean <= 1.1)).sum()),
    f"({100 * ((L.r_mean > .9) & (L.r_mean <= 1.1)).mean():.2f} %)")
say("relative to the image MEDIAN: <0.85:", int((L.r_med < .85).sum()), "| 0.85-1.15:",
    int(((L.r_med >= .85) & (L.r_med <= 1.15)).sum()), "| 1.15-1.6:", int(((L.r_med > 1.15) & (L.r_med <= 1.6)).sum()),
    "| >1.6:", int((L.r_med > 1.6).sum()), "| smallest ratio above 1.15:", round(L.r_med[L.r_med > 1.15].min(), 3),
    "| largest ratio in core:", round(L.r_med[L.r_med <= 1.15].max(), 3))
core = L[(L.r_med >= .85) & (L.r_med <= 1.15)]
say("SD of core spacings relative to image median [%]:", (core.groupby("detector").r_med.std() * 100).round(2).to_dict())
p["min_med"] = L.groupby("image_path").r_med.min().reindex(p.image_path).values
p["max_med"] = L.groupby("image_path").r_med.max().reindex(p.image_path).values
p["two_sided"] = (p.min_med >= .85) & (p.max_med <= 1.15)
p["kind"] = np.where(p.two_sided, "regular", np.where(p.max_med > 1.15, "missed minimum", "additional minimum"))
say(pd.crosstab([p.group_label, p.detector], [p.accepted, p.kind]).to_string())
t = p.pivot_table(index=["group_label", "tile_x", "tile_y"], columns="detector", values="two_sided", aggfunc="first")
say("tiles with irregular spacing: both detectors", int((~t["SED"] & ~t["BED-S"]).sum()), "| SED only",
    int((~t["SED"] & t["BED-S"]).sum()), "| BED-S only", int((t["SED"] & ~t["BED-S"]).sum()))
irr = p[~p.two_sided]
say("FFT pitch of the", len(irr), "irregular images: median", round(irr.pitch_fft_um.median(), 4), "range",
    round(irr.pitch_fft_um.min(), 4), "-", round(irr.pitch_fft_um.max(), 4), "| regular images: median",
    round(p[p.two_sided].pitch_fft_um.median(), 4), "range", round(p[p.two_sided].pitch_fft_um.min(), 4), "-",
    round(p[p.two_sided].pitch_fft_um.max(), 4))
ao = p[p.accepted & ~p.two_sided]
say("accepted by the submitted rule but irregular:", len(ao), ao.groupby("detector").size().to_dict(),
    "| their 1D pitch range", round(ao.pitch_1d_um.min(), 3), "-", round(ao.pitch_1d_um.max(), 3))

# ---------------------------------------------------------------- B
say("\n== B. Influence of the acceptance rule on the 1D medians (um) ==")
rules = [("none", 0.0, np.inf), ("min>=0.5", 0.5, np.inf), ("min>=0.6 (submitted)", 0.6, np.inf),
         ("min>=0.7", 0.7, np.inf), ("min>=0.8", 0.8, np.inf), ("0.6<=r<=1.4", 0.6, 1.4),
         ("two-sided +-15 % of median", None, None)]
ref = p[p.accepted].groupby(["group_label", "detector"]).pitch_1d_um.median()
out = []
for name, lo, hi in rules:
    ok = p.two_sided if lo is None else (p.min_spacing_ratio >= lo) & (p.max_spacing_ratio <= hi)
    g = p[ok].groupby(["group_label", "detector"]).pitch_1d_um
    med = g.median()
    dev = (p[ok].pitch_1d_um - p[ok].groupby(["group_label", "detector"]).pitch_1d_um.transform("median")).abs()
    out.append(dict(rule=name, accepted=int(ok.sum()), rejected=int((~ok).sum()),
                    max_abs_median_change_nm=round(float((med - ref).abs().max() * 1e3), 2),
                    pooled_median=round(float(p[ok].pitch_1d_um.median()), 5),
                    max_abs_dev_from_median_um=round(float(dev.max()), 3),
                    n_dev_gt_0p09=int((dev > 0.09).sum())))
B = pd.DataFrame(out)
B.to_csv(os.path.join(OUT, "referee_rule_sensitivity.csv"), index=False)
say(B.to_string(index=False))

# ---------------------------------------------------------------- C
say("\n== C. Smoothing sensitivity (submitted acceptance rule) ==")
v = pd.read_csv(os.path.join(OUT, "smoothing_variants_per_image.csv"))
v = v.merge(p[["image_path", "group_label", "detector"]], on="image_path")
v["ok"] = v.min_over_mean >= 0.6
v["pitch"] = v.period_px * S
base = v[(v.setting == "spar_0.10") & v.ok].groupby(["group_label", "detector"]).pitch.median()
out = []
for name, g in v.groupby("setting"):
    med = g[g.ok].groupby(["group_label", "detector"]).pitch.median()
    d = (med - base)
    out.append(dict(setting=name, edf=round(float(g.df.median()), 0) if g.df.notna().any() else np.nan,
                    accepted=int(g.ok.sum()), accepted_SED=int(g[g.detector == "SED"].ok.sum()),
                    accepted_BED=int(g[g.detector == "BED-S"].ok.sum()),
                    max_dmed_SED_nm=round(float(d.xs("SED", level="detector").abs().max() * 1e3), 2),
                    max_dmed_BED_nm=round(float(d.xs("BED-S", level="detector").abs().max() * 1e3), 2)))
C = pd.DataFrame(out).set_index("setting").loc[
    ["none", "spar_0.01", "spar_0.05", "spar_0.10", "spar_0.20", "spar_0.40", "spar_0.70", "GCV"]].reset_index()
C.to_csv(os.path.join(OUT, "referee_smoothing_sensitivity.csv"), index=False)
say(C.to_string(index=False))

# ---------------------------------------------------------------- D
say("\n== D. Figure 6 / Table 5 of the submitted version (pitch1d_physical.csv) ==")
f6 = pd.read_csv(os.path.join(RES, "pitch1d_physical.csv"))
for mag in MAGS:
    g = f6[f6.mag == mag].dropna(subset=["pitch"])
    bx, by = stats.linregress(g.stage_x_um, g.pitch), stats.linregress(g.stage_y_um, g.pitch)
    say(f"{mag}x: N={len(g)} of {int((f6.mag == mag).sum())}  bx={bx.slope:+.1e} (p={bx.pvalue:.2g})  by={by.slope:+.1e} "
        f"(p={by.pvalue:.2g})  median={g.pitch.median():.4f}  SD={g.pitch.std():.4f}  min={g.pitch.min():.3f} max={g.pitch.max():.3f}")
f6["dev"] = f6.pitch - f6.groupby("mag").pitch.transform("median")
f6["image_path"] = "mag" + f6.mag.astype(str) + "_SED-BED/" + f6.det + "/" + f6.base + ".png"
o = f6[f6.dev.abs() > 0.03].merge(
    p[["image_path", "pitch_1d_um", "pitch_fft_um", "accepted", "kind", "min_med", "max_med"]], on="image_path")
o = o.sort_values("dev", key=np.abs, ascending=False)
o.to_csv(os.path.join(OUT, "referee_fig6_outliers.csv"), index=False)
say("points of Fig. 6 deviating by more than 0.03 um from the series median:", len(o),
    o.groupby(["mag", "det"]).size().to_dict())
say(o[["mag", "det", "base", "pitch", "pitch_1d_um", "pitch_fft_um", "accepted", "kind"]].round(4).to_string(index=False))

say("\n== E. Figure 5 selection (linearity_tiles_mag*.csv, SED frames, 3 scaled MAD) ==")
for mag in MAGS:
    d = pd.read_csv(os.path.join(RES, f"linearity_tiles_mag{mag}.csv"))
    med = d.local_pitch_um.median()
    mad = (d.local_pitch_um - med).abs().median() * 1.4826
    say(f"{mag}x: SED frames {len(d)}, kept {int(((d.local_pitch_um - med).abs() < 3 * mad).sum())}; 3 MAD = {3 * mad * 1e3:.1f} nm")

# ---------------------------------------------------------------- F
say("\n== F. Stage-position and time regressions of per-tile values ==")
tile = p.groupby(["group_label", "tile_x", "tile_y"]).agg(
    x=("stage_x_um", "mean"), y=("stage_y_um", "mean"), t=("t_min", "mean"),
    fft=("pitch_fft_um", "mean")).reset_index()
reg1d = p[p.two_sided].groupby(["group_label", "tile_x", "tile_y"]).pitch_1d_um.mean().rename("d1").reset_index()
tile = tile.merge(reg1d, on=["group_label", "tile_x", "tile_y"], how="left")
out = []
for mag in MAGS:
    g = tile[tile.group_label == mag]
    dur = g.t.max() - g.t.min()
    for col, lab in [("fft", "FFT, all frames"), ("d1", "1D, regular frames")]:
        gg = g.dropna(subset=[col])
        bx, by, bt = (stats.linregress(gg[c], gg[col]) for c in ["x", "y", "t"])
        out.append(dict(series=mag, estimator=lab, tiles=len(gg), dx_um=round(np.ptp(gg.x)), dy_um=round(np.ptp(gg.y)),
                        duration_min=round(dur, 1), corr_t_y=round(np.corrcoef(g.t, g.y)[0, 1], 3),
                        corr_t_x=round(np.corrcoef(g.t, g.x)[0, 1], 3),
                        bx=bx.slope, px=bx.pvalue, by=by.slope, py=by.pvalue,
                        rel_x_pct=100 * abs(bx.slope) * np.ptp(gg.x) / gg[col].mean(),
                        rel_y_pct=100 * abs(by.slope) * np.ptp(gg.y) / gg[col].mean(),
                        rate_nm_per_h=bt.slope * 6e4, sd_nm=gg[col].std() * 1e3,
                        sd_detrended_nm=(gg[col] - bt.intercept - bt.slope * gg.t).std() * 1e3))
F = pd.DataFrame(out)
F.to_csv(os.path.join(OUT, "referee_spatial_time_regressions.csv"), index=False)
pd.set_option("display.width", 250)
say(F.assign(bx=F.bx.map("{:+.1e}".format), by=F.by.map("{:+.1e}".format), px=F.px.map("{:.2g}".format),
             py=F.py.map("{:.2g}".format)).round(2).to_string(index=False))
say("pooled per-image regressions, submitted rule (1D accepted) for comparison:")
for mag in MAGS:
    g = p[(p.group_label == mag) & p.accepted]
    bx, by = stats.linregress(g.stage_x_um, g.pitch_1d_um), stats.linregress(g.stage_y_um, g.pitch_1d_um)
    say(f"  {mag}x: N={len(g)} bx={bx.slope:+.1e} (p={bx.pvalue:.2g}) by={by.slope:+.1e} (p={by.pvalue:.2g})")
first = tile.sort_values("t").groupby("group_label").head(3)[["group_label", "tile_x", "tile_y", "t"]]
say("first three tiles in time order per series:", first.values.tolist())

say("\n== G. Common central area (stage window of the 1000x series) ==")
a = tile[tile.group_label == 1000]
x0, x1, y0, y1 = a.x.min() - 2, a.x.max() + 2, a.y.min() - 2, a.y.max() + 2
for mag in MAGS:
    g = tile[(tile.group_label == mag) & tile.x.between(x0, x1) & tile.y.between(y0, y1)]
    allg = tile[tile.group_label == mag]
    say(f"{mag}x: tiles in window {len(g)}; FFT median in window {g.fft.median():.4f} (whole series {allg.fft.median():.4f}); "
        f"centre of series ({allg.x.mean():.0f}, {allg.y.mean():.0f}) um")

# ---------------------------------------------------------------- H
say("\n== H. Position of the spacing within the frame (regular spacings, all series) ==")
core = core.assign(bin=pd.cut(core.centre_px, [0, 320, 640, 960, 1280]))
say(((core.groupby(["bin", "detector"], observed=True).r_med.mean() - 1) * 100).round(3).unstack().to_string())
for det, g in core.groupby("detector"):
    r = stats.linregress(g.centre_px, g.px)
    lo, hi = r.slope - 1.96 * r.stderr, r.slope + 1.96 * r.stderr
    say(f"{det}: n={len(g)} change of spacing across 1280 px = {100 * r.slope * 1280 / g.px.mean():+.2f} % "
        f"[{100 * lo * 1280 / g.px.mean():+.2f}, {100 * hi * 1280 / g.px.mean():+.2f}]")

# ---------------------------------------------------------------- I
say("\n== I. Paired SED - BED-S differences of the 1D pitch (um) ==")
for label, flag in [("submitted rule", "accepted"), ("two-sided rule", "two_sided")]:
    w = p[p[flag]].pivot_table(index=["group_label", "tile_x", "tile_y"], columns="detector", values="pitch_1d_um")
    w = w.dropna()
    d = w["SED"] - w["BED-S"]
    for mag in MAGS:
        dd = d.xs(mag, level="group_label")
        say(f"{label} {mag}x: pairs={len(dd)} bias={dd.mean():+.4f} SD={dd.std():.4f} "
            f"LoA=[{dd.mean() - 1.96 * dd.std():+.4f}, {dd.mean() + 1.96 * dd.std():+.4f}] median={dd.median():+.4f} "
            f"p={stats.ttest_1samp(dd, 0).pvalue:.2f}")

# ---------------------------------------------------------------- J
say("\n== J. Anomalous spacings: type and detection by the submitted rule ==")
img = L.groupby("image_path").agg(min_med=("r_med", "min"), max_med=("r_med", "max"),
                                  min_mean=("r_mean", "min"), n=("n", "first")).reset_index()
img = img.merge(p[["image_path", "detector", "accepted"]], on="image_path")
short, long_ = img[img.min_med < .85], img[img.max_med > 1.15]
say("largest ratio to the image median:", round(L.r_med.max(), 3), "| smallest:", round(L.r_med.min(), 3))
say("images with a short spacing (<0.85 of the median):", len(short), short.groupby("detector").size().to_dict(),
    "| smallest spacing below 0.6 of the image mean (rejected):", int((short.min_mean < .6).sum()),
    "| not rejected:", int((short.min_mean >= .6).sum()))
say("images with a long spacing (>1.15 of the median):", len(long_), long_.groupby("detector").size().to_dict(),
    "| also containing a short spacing:", int((long_.min_med < .85).sum()),
    "| rejected by the submitted rule:", int((~long_.accepted).sum()))


def is_split(g):
    r = g.sort_values("k").r_med.values
    k = int(np.argmin(r))
    return any(0.85 <= r[k] + r[j] <= 1.15 for j in (k - 1, k + 1) if 0 <= j < len(r))


split = L[L.image_path.isin(short.image_path)].groupby("image_path").apply(is_split)
short = short.merge(split.rename("split").reset_index(), on="image_path")
say("short-spacing images in which the short spacing and a neighbour add up to one period (split period):",
    int(short.split.sum()), "of", len(short), "| largest value of the smallest spacing relative to the image mean in "
    "these images:", round(short[short.split].min_mean.max(), 3), "| rejected:", int((short[short.split].min_mean < .6).sum()))
say("short-spacing images that are not a split period:", int((~short.split).sum()),
    "| of these at the first or last spacing:",
    int(sum(L[L.image_path == q].sort_values("k").r_med.values.argmin() in (0, int(L[L.image_path == q].n.iloc[0]) - 1)
            for q in short[~short.split].image_path)), "| rejected:", int((short[~short.split].min_mean < .6).sum()))

# ---------------------------------------------------------------- K
say("\n== K. Row-wise test of the y trend (one mean value per tile row) ==")
for mag in MAGS:
    g = tile[tile.group_label == mag]
    for col, lab in [("fft", "FFT"), ("d1", "1D")]:
        rows = g.dropna(subset=[col]).groupby("tile_y").agg(y=("y", "mean"), v=(col, "mean"))
        r = stats.linregress(rows.y, rows.v)
        say(f"{mag}x {lab}: rows={len(rows)} slope={r.slope:+.2e} um/um  p={r.pvalue:.1e}  r={r.rvalue:+.3f}")

# ---------------------------------------------------------------- L
say("\n== L. Pairing of the detector images and overall spread of the FFT pitch ==")
w = p.pivot_table(index=["group_label", "tile_x", "tile_y"], columns="detector",
                  values=["stage_x_um", "stage_y_um", "t_min"])
say("largest stage-coordinate difference between the SED and BED-S image of a tile [um]: x",
    float((w["stage_x_um"]["SED"] - w["stage_x_um"]["BED-S"]).abs().max()), "y",
    float((w["stage_y_um"]["SED"] - w["stage_y_um"]["BED-S"]).abs().max()))
dt = (w["t_min"]["BED-S"] - w["t_min"]["SED"]) * 60
say(f"BED-S minus SED acquisition time [s]: median {dt.median():.0f}, min {dt.min():.0f}, max {dt.max():.0f}")
step = p[p.detector == "SED"].sort_values(["group_label", "time"]).groupby("group_label").time.diff().dt.total_seconds()
say(f"time between consecutive SED images [s]: median {step.median():.0f}, 1st-99th percentile "
    f"{step.quantile(.01):.0f}-{step.quantile(.99):.0f}")
say(f"FFT pitch, all 1500 images: min {p.pitch_fft_um.min():.4f} max {p.pitch_fft_um.max():.4f} "
    f"SD {p.pitch_fft_um.std():.4f} ({100 * p.pitch_fft_um.std() / p.pitch_fft_um.median():.2f} %) "
    f"range {100 * (p.pitch_fft_um.max() - p.pitch_fft_um.min()) / p.pitch_fft_um.median():.2f} %")
say("FFT medians per series:", p.groupby("group_label").pitch_fft_um.median().round(4).to_dict())

with open(os.path.join(OUT, "referee_analyses.txt"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
