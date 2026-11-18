#!/usr/bin/env python3
"""Figures for the revision of MST-139559.

Manuscript figures (new file names; submitted figures are not overwritten):
  figures/linearity_combined_rev1.pdf   Fig. 6 with insets (same data and fits)
  figures/method_comparison_rev1.pdf    Fig. 11 without the nominal reference line
  figures/acquisition_schematic.pdf     schematic requested for Sections 2.6-2.7
Figures for the response letter:
  response/figR_spacing_ratio.pdf, figR_outliers.pdf, figR_time_position.pdf
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Arc, Rectangle
from PIL import Image
from scipy import stats
from scipy.interpolate import UnivariateSpline

HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(HERE)
ROOT = os.path.dirname(REV)
AN = os.path.join(REV, "analysis")
FIG = os.path.join(REV, "figures")
RESP = os.path.join(REV, "response")
# Submitted result tables: R_output_dir/ in the working folder, top level in
# the public repository.
RES = os.path.join(ROOT, "R_output_dir")
if not os.path.isdir(RES):
    RES = ROOT
os.makedirs(FIG, exist_ok=True)
os.makedirs(RESP, exist_ok=True)
S = 0.01333594
MAGS = [400, 500, 750, 1000]
plt.rcParams.update({"font.size": 9, "pdf.fonttype": 42})


def report_text_overlaps(fig, name):
    """Print pairs of drawn text labels whose bounding boxes intersect."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = list(fig.texts) + [t for leg in fig.legends for t in leg.get_texts()]
    for ax in fig.findobj(matplotlib.axes.Axes):
        texts += [ax.title, ax.xaxis.label, ax.yaxis.label] + list(ax.texts)
        if ax.get_legend() is not None:
            texts += ax.get_legend().get_texts()
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            texts += [lab for loc, lab in zip(axis.get_majorticklocs(), axis.get_majorticklabels()) if lo <= loc <= hi]
    boxes = [(t.get_text().replace("\n", " ")[:28], t.get_window_extent(renderer))
             for t in texts if t.get_text().strip() and t.get_visible()]
    hits = [(a, b) for i, (a, ba) in enumerate(boxes) for b, bb in boxes[i + 1:] if ba.overlaps(bb)]
    print(f"{name}: {len(boxes)} labels, {len(hits)} overlapping pairs", hits[:6])

# ------------------------------------------------------------ Fig. 6
# Layout and styling of the submitted figure (sem_linearity_fig.py); only the insets are new.
df = pd.read_csv(os.path.join(RES, "pitch1d_physical.csv")).dropna(subset=["pitch"])
fig, axes = plt.subplots(2, 2, figsize=(8.6, 7.6))
handles = []
for ax, mag in zip(axes.ravel(), MAGS):
    g = df[df.mag == mag]
    bx = stats.linregress(g.stage_x_um, g.pitch)
    xx = np.linspace(g.stage_x_um.min(), g.stage_x_um.max(), 50)
    dbar = g.pitch.mean()
    xr = g.stage_x_um.max() - g.stage_x_um.min()
    targets = [ax]
    if mag in (400, 500):
        ins = ax.inset_axes([0.17, 0.40, 0.80, 0.50])
        ins.set_ylim(1.60, 1.72)
        ins.tick_params(labelsize=7)
        ins.text(0.02, 0.05, "enlarged ordinate", transform=ins.transAxes, fontsize=7)
        targets.append(ins)
    for a in targets:
        for det, c, m in [("SED", "#1f77b4", "o"), ("BED-S", "#d62728", "s")]:
            gg = g[g.det == det]
            hd = a.plot(gg.stage_x_um, gg.pitch, m, ms=4 if a is ax else 3, mfc="none", mec=c, mew=0.7, alpha=0.75)[0]
            if mag == 400 and a is ax:
                handles.append((hd, det))
        a.plot(xx, bx.intercept + bx.slope * xx, "k-", lw=1.2)
        a.axhline(dbar, color="0.45", ls="--", lw=0.8)
    ax.set_title(f"${mag}\\times$:  $b_x = {bx.slope:.1e}$ µm/µm ($p = {bx.pvalue:.2f}$)\n"
                 f"max. rel. change $= {100 * abs(bx.slope) * xr / dbar:.2f}$%", fontsize=9.5)
    ax.set_xlabel("Stage position $x$ [µm]", fontsize=9)
    ax.set_ylabel("Local groove spacing [µm]", fontsize=9)
    ax.tick_params(labelsize=8)
fig.legend([hd for hd, _ in handles], [det for _, det in handles], loc="upper center", ncol=2, fontsize=9, frameon=False)
fig.tight_layout(rect=[0, 0, 1, 0.96])
report_text_overlaps(fig, "linearity_combined_rev1")
for ax, mag in zip(axes.ravel(), MAGS):
    g = df[df.mag == mag]
    lo, hi = ax.get_ylim()
    hidden = g[(g.pitch > lo + 0.36 * (hi - lo)) & (g.pitch < lo + 0.93 * (hi - lo))] if mag in (400, 500) else g.iloc[:0]
    print(f"  {mag}x: y-limits {lo:.2f}-{hi:.2f}; points behind the inset: {len(hidden)}; "
          f"points outside the inset range 1.60-1.72: {int(((g.pitch < 1.60) | (g.pitch > 1.72)).sum())}")
fig.savefig(os.path.join(FIG, "linearity_combined_rev1.pdf"))
plt.close(fig)

# ----------------------------------------------------------- Fig. 11
labels = ["Nominal\n(design value,\nno uncertainty\nstated)", "SEM 1D\n(instrument-\nreported scale)",
          "SEM FFT\n(instrument-\nreported scale)", "Optical\ndiffraction\n(laser wavelength,\nruler, caliper)",
          "WLI\n(stage scales\nof the CMM)"]
pitch = [1.6667, 1.6448, 1.6442, 1.650, 1.66512]
yerr = [None, 0.09, 0.09, 0.21, 0.0016]
colours = ["black", "tab:blue", "tab:blue", "tab:green", "tab:orange"]
markers = ["s", "o", "o", "D", "o"]
# Title, axis label, markers and font size as in the submitted figure (method_comparison_pitch.py).
with plt.rc_context({"font.size": 10}):
    fig, ax = plt.subplots(figsize=(8, 4.9))
    for i, (v, e, c, mk) in enumerate(zip(pitch, yerr, colours, markers)):
        if e is None:
            ax.plot(i, v, mk, mfc="none", mec=c, ms=9, mew=1.4)
        else:
            ax.errorbar(i, v, yerr=e, fmt=mk, color=c, ms=9, capsize=6, lw=1.6)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("pitch [µm]")
    ax.set_ylim(1.40, 1.90)
    ax.set_title("Grating pitch: method comparison", fontsize=11)
    ax.grid(axis="y", ls=":", alpha=0.5)
    fig.tight_layout()
    report_text_overlaps(fig, "method_comparison_rev1")
    fig.savefig(os.path.join(FIG, "method_comparison_rev1.pdf"))
    plt.close(fig)

# --------------------------------------------------------- schematic
fig, (a, b) = plt.subplots(1, 2, figsize=(8.6, 3.9), gridspec_kw={"width_ratios": [1.1, 1]})
nx, ny, w, h = 8, 8, 17.0, 12.6
for i in range(nx):
    for j in range(ny):
        a.add_patch(Rectangle((-(i + 1) * w, j * h), w, h, fc="0.93" if (i + j) % 2 else "white", ec="0.4", lw=0.6))
# Recorded path: tile (0,0) first (largest stage x), x direction alternating from row to row.
x_first, x_last = -w / 2, -nx * w + w / 2
for j in range(3):
    yc = j * h + h / 2
    x0, x1 = (x_first, x_last) if j % 2 == 0 else (x_last, x_first)
    a.annotate("", xy=(x1, yc), xytext=(x0, yc), arrowprops=dict(arrowstyle="->", color="tab:red", lw=1.1))
    a.plot([x1, x1], [yc, yc + h], color="tab:red", lw=1.1, ls="-" if j < 2 else ":")
a.text(-nx * w / 2, 3 * h + h / 2, "acquisition path (first three rows)", color="tab:red", fontsize=7, ha="center",
       va="center", bbox=dict(fc="white", ec="none", pad=0.8))
a.text(0, -4, "tile (0,0): first image", ha="right", va="top", fontsize=7)
a.text(-nx * w / 2, ny * h + 3,
       "survey area covered by the montage\n(8$\\times$8 images shown; steps of about 17 µm $\\times$ 12.6 µm)",
       ha="center", fontsize=8)
a.annotate("", xy=(-nx * w - 8, ny * h), xytext=(-nx * w - 8, 0), arrowprops=dict(arrowstyle="->", lw=1))
a.text(-nx * w - 11, ny * h / 2, "stage $y$ (acquisition order of rows)", rotation=90, va="center", ha="right", fontsize=8)
a.annotate("", xy=(0, -14), xytext=(-nx * w, -14), arrowprops=dict(arrowstyle="->", lw=1))
a.text(-nx * w / 2, -17, "stage $x$", ha="center", va="top", fontsize=8)
a.add_patch(Rectangle((-4 * w, 4 * h), w, h, fc="none", ec="tab:blue", lw=1.8))
a.text(-3.5 * w, 5 * h + 1.5, "one image", color="tab:blue", fontsize=7, ha="center")
a.set_xlim(-nx * w - 30, 8)
a.set_ylim(-26, ny * h + 20)
a.set_aspect("equal")
a.set_anchor("N")
a.axis("off")
a.set_title("(a)", loc="left", fontsize=10)

W, H, per, ang = 17.07, 12.80, 1.644, 6.0
b.add_patch(Rectangle((0, 0), W, H, fc="white", ec="tab:blue", lw=1.8))
t = np.tan(np.radians(ang))
K_RED = 4
for k in range(-2, 12):
    x0 = k * per + 0.6
    xs = np.array([x0, x0 + H * t])
    if ((xs >= 0) & (xs <= W)).all() and k != K_RED:
        b.plot(xs, [H, 0], color="0.6", lw=0.9)
# Both angles are drawn at the upper end T of one groove; image y points downwards.
xt = K_RED * per + 0.6
white = dict(fc="white", ec="none", pad=0.6)
b.plot([xt, xt + H * t], [H, 0], color="tab:red", lw=1.8)
b.plot([xt, xt], [H, H - 10.5], color="k", lw=0.8, ls="--")
b.add_patch(Arc((xt, H), 4.4, 4.4, theta1=-(90 - ang), theta2=0, color="k", lw=0.9))
b.text(xt + 2.5, H - 1.75, "$\\bar{\\alpha}$", fontsize=9, va="center", bbox=white)
b.add_patch(Arc((xt, H), 17, 17, theta1=-90, theta2=-(90 - ang), color="tab:red", lw=1.1))
b.text(xt - 0.3, H - 9.3, "$\\Delta\\alpha$", color="tab:red", fontsize=9, va="center", ha="right", bbox=white)
b.annotate("", xy=(xt + 4.6, H + 0.55), xytext=(xt, H + 0.55), arrowprops=dict(arrowstyle="->", lw=0.8))
b.text(xt + 4.8, H + 0.55, "image $x$", fontsize=7.5, va="center")
b.annotate("", xy=(W + 0.6, H - 4.0), xytext=(W + 0.6, H), arrowprops=dict(arrowstyle="->", lw=0.8))
b.text(W + 0.95, H - 2.0, "image $y$", fontsize=7.5, va="center", rotation=90)
b.annotate("", xy=(W, -1.2), xytext=(0, -1.2), arrowprops=dict(arrowstyle="<->", lw=0.8))
b.text(W / 2, -1.6, "17.07 µm (1280 px): profile and FFT axis", ha="center", va="top", fontsize=7.5)
b.annotate("", xy=(-1.0, 0), xytext=(-1.0, H), arrowprops=dict(arrowstyle="<->", lw=0.8))
b.text(-1.4, H / 2, "12.80 µm (960 px): averaging direction", rotation=90, ha="right", va="center", fontsize=7.5)
b.annotate("", xy=(0.6 + per, H + 0.55), xytext=(0.6, H + 0.55), arrowprops=dict(arrowstyle="<->", lw=0.8))
b.text(0.6 + per / 2, H + 0.95, "$d$", fontsize=9, ha="center")
b.text(W, -3.4, "groove rotation exaggerated", fontsize=7, color="0.3", ha="right")
b.set_xlim(-4, W + 2.6)
b.set_ylim(-4.4, H + 2.2)
b.set_aspect("equal")
b.set_anchor("N")
b.axis("off")
b.set_title("(b)", loc="left", fontsize=10)
fig.tight_layout()
report_text_overlaps(fig, "acquisition_schematic")
fig.savefig(os.path.join(FIG, "acquisition_schematic.pdf"))
plt.close(fig)

# ----------------------------------------- response: spacing ratios
p = pd.read_csv(os.path.join(AN, "profiles_per_image.csv")).set_index("image_path")
ratios = {"SED": [], "BED-S": []}
for r in p.itertuples():
    s = np.diff([float(v) for v in str(r.minima_pixels_1based).split(";")])
    ratios[r.detector].extend(s / np.median(s))
fig, ax = plt.subplots(figsize=(7, 3.4))
edges = np.arange(0, 3.0501, 0.05)  # the largest ratio is 2.98
ax.hist([ratios["SED"], ratios["BED-S"]], bins=edges, label=["SED", "BED-S"], color=["#1f77b4", "#d62728"], log=True)
for v, lab in [(0.6, "0.6 (acceptance criterion)"), (0.85, "0.85"), (1.15, "1.15")]:
    ax.axvline(v, color="k", ls="--" if v == 0.6 else ":", lw=0.9)
    ax.text(v, 2.2e4, lab, rotation=90, va="top", ha="right", fontsize=7)
ax.set_xlabel("Inter-minimum spacing / median spacing of the same image")
ax.set_ylabel("Number of spacings")
ax.set_ylim(0.7, 3e4)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(RESP, "figR_spacing_ratio.pdf"))
plt.close(fig)


# -------------------------------------------- response: outlier frames
def fig6_algorithm(prof):
    x = np.arange(len(prof))
    sp = UnivariateSpline(x, prof, s=0.5 * np.sum((prof - prof.mean()) ** 2))(x)
    d = np.diff(sp)
    minima = np.where((d[:-1] < 0) & (d[1:] > 0))[0] + 1
    dist = np.diff(minima) * S
    keep = (dist >= 0.7 * np.median(dist)) & (dist <= 1.3 * np.median(dist)) if len(dist) else np.array([], bool)
    return sp, minima, dist, keep


out = pd.read_csv(os.path.join(AN, "referee_fig6_outliers.csv"))
PER_PAGE = 7
detail = []
for page in range(int(np.ceil(len(out) / PER_PAGE))):
    part = out.iloc[page * PER_PAGE:(page + 1) * PER_PAGE]
    fig, axes = plt.subplots(PER_PAGE, 2, figsize=(8.6, 2.05 * PER_PAGE), gridspec_kw={"width_ratios": [1, 2.3]})
    row_titles = []
    for (a, b), r in zip(axes, part.itertuples()):
        img = np.asarray(Image.open(os.path.join(ROOT, r.image_path)).convert("L"), dtype=float)
        prof = img.mean(axis=0)
        sp, minima, dist, keep = fig6_algorithm(prof)
        main = np.array([float(v) for v in str(p.loc[r.image_path, "minima_pixels_1based"]).split(";")]) - 1
        a.imshow(img, cmap="gray", vmin=0, vmax=255)
        a.set_xticks([]); a.set_yticks([])
        a.set_ylabel(f"{r.mag}$\\times$ {r.det}\n{r.base.split('_', 3)[-1]}", fontsize=7)
        b.plot(prof, color="0.6", lw=0.6)
        b.plot(sp, color="tab:blue", lw=1.0)
        b.plot(minima, sp[minima], "v", color="tab:blue", ms=5)
        lo = prof.min() - 0.08 * np.ptp(prof)
        b.plot(main, np.full(len(main), lo), "|", color="tab:red", ms=9, mew=1.2)
        b.set_xlim(0, len(prof))
        b.tick_params(labelsize=7)
        mratio = np.diff(main) / np.median(np.diff(main))
        irregular = [f"{n} {k} spacing{'s' if n > 1 else ''}"
                     for n, k in ((int((mratio < 0.85).sum()), "short"), (int((mratio > 1.15).sum()), "long")) if n]
        status = "regular" if not irregular else (
            ", ".join(irregular) + ("; accepted" if p.loc[r.image_path, "accepted"] else "; rejected"))
        row_titles.append((a, b, f"Fig. 6 value {r.pitch:.3f} µm ({int(keep.sum())} of {len(dist)} spacings kept)   |   "
                                 f"main analysis {r.pitch_1d_um:.3f} µm ({status})   |   FFT {r.pitch_fft_um:.3f} µm"))
        b.set_title(" ", fontsize=7.5)  # reserves the height of the row title
        ratio = dist / np.median(dist)
        detail.append(dict(image_path=r.image_path, fig6_minima=len(minima), fig6_spacings=len(dist),
                           fig6_kept=int(keep.sum()), n_short=int((ratio < 0.7).sum()), n_long=int((ratio > 1.3).sum()),
                           sd_kept_pct=round(float(100 * np.std(ratio[keep])), 1) if keep.sum() > 1 else np.nan,
                           main_minima=len(main),
                           spacings_um=" ".join(f"{v:.2f}" for v in dist)))
    for a, b in axes[len(part):]:
        a.axis("off"); b.axis("off")
    axes[len(part) - 1, 1].set_xlabel("Image $x$ [px]", fontsize=8)
    fig.tight_layout()
    fig.canvas.draw()  # fixes the positions of the equal-aspect image panels
    for a, b, text in row_titles:
        fig.text(a.get_position().x0, b.get_position().y1 + 0.004, text, fontsize=7.5, ha="left", va="bottom")
    fig.savefig(os.path.join(RESP, f"figR_outliers_{page + 1}.pdf"), dpi=200)
    plt.close(fig)
detail = pd.DataFrame(detail)
detail.to_csv(os.path.join(AN, "referee_fig6_outliers_detail.csv"), index=False)
print(detail.drop(columns="image_path").to_string())

# ------------- response: the three most deviating images, both detectors
TOP3 = [(400, "Grp_002", "X008_Y010"), (400, "Grp_002", "X008_Y011"), (500, "Grp_001", "X014_Y013")]
fig, axes = plt.subplots(2, 3, figsize=(8.6, 4.75))
bar = 2.0 / S
for col, (mag, grp, tile_id) in enumerate(TOP3):
    for row, det in enumerate(["SED", "BED-S"]):
        ax = axes[row, col]
        img = np.asarray(Image.open(os.path.join(ROOT, f"mag{mag}_SED-BED", det, f"{grp}_{det}_{tile_id}.png")).convert("L"))
        ax.imshow(img, cmap="gray", vmin=0, vmax=255)
        ax.set_xticks([]); ax.set_yticks([])
        ax.add_patch(Rectangle((35, img.shape[0] - 150), bar + 50, 115, fc="black", ec="none", alpha=0.65))
        ax.plot([60, 60 + bar], [img.shape[0] - 62, img.shape[0] - 62], color="white", lw=2.5, solid_capstyle="butt")
        ax.text(60 + bar / 2, img.shape[0] - 78, "2 µm", color="white", fontsize=7, ha="center", va="bottom")
        if row == 0:
            ax.set_title(f"{mag}$\\times$ series, tile {tile_id.replace('_', '/')}", fontsize=8.5)
        if col == 0:
            ax.set_ylabel(det, fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(RESP, "figR_outliers_top3.pdf"), dpi=220)
plt.close(fig)

# ------------------- contact sheets of the images with irregular spacings
irr = []
for r in p.reset_index().itertuples():
    m = np.array([float(v) for v in str(r.minima_pixels_1based).split(";")]) - 1
    s = np.diff(m)
    bad = np.where((s / np.median(s) < 0.85) | (s / np.median(s) > 1.15))[0]
    if len(bad):
        irr.append((r.image_path, r.group_label, r.detector, r.tile_x, r.tile_y, [(m[k], m[k + 1]) for k in bad]))
print("images with irregular spacings:", len(irr))
for sheet in range(int(np.ceil(len(irr) / 25))):
    fig, axes = plt.subplots(5, 5, figsize=(15, 11.6))
    for ax, item in zip(axes.ravel(), irr[sheet * 25:(sheet + 1) * 25]):
        path, grp, det, tx, ty, spans = item
        img = Image.open(os.path.join(ROOT, path)).convert("L").resize((640, 480))
        ax.imshow(img, cmap="gray")
        for x0, x1 in spans:
            ax.plot([x0 / 2, x1 / 2], [12, 12], color="red", lw=3)
        ax.set_title(f"{grp}x {det} X{tx:03d}/Y{ty:03d}", fontsize=8)
    for ax in axes.ravel():
        ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(os.path.join(AN, f"irregular_images_sheet{sheet + 1}.png"), dpi=70)
    plt.close(fig)

# --------------------------------- response: time versus stage position
pm = p.reset_index()
pm["file"] = pm.image_path.str.split("/").str[-1]
fft = pd.read_csv(os.path.join(RES, "fft_subpixel_per_image.csv")).rename(
    columns={"mag": "group_label", "det": "detector", "pitch": "fft"})
pm = pm.merge(fft[["group_label", "detector", "file", "fft"]], on=["group_label", "detector", "file"])
pm["time"] = pd.to_datetime(pm.timestamp)
pm["t"] = pm.groupby("group_label").time.transform(lambda v: (v - v.min()).dt.total_seconds() / 60)
tile = pm.groupby(["group_label", "tile_x", "tile_y"]).agg(x=("stage_x_um", "mean"), t=("t", "mean"), fft=("fft", "mean")).reset_index()


def ptext(pv):
    if pv >= 0.01:
        return f"{pv:.2f}"
    mant, expo = f"{pv:.0e}".split("e")
    return f"{mant}\\times10^{{{int(expo)}}}"


fig, axes = plt.subplots(2, 4, figsize=(11, 4.8), sharey=True)
for k, mag in enumerate(MAGS):
    g = tile[tile.group_label == mag]
    for row, (col, lab) in enumerate([("t", "Acquisition time [min]"), ("x", "Stage position $x$ [µm]")]):
        ax = axes[row, k]
        ax.plot(g[col], g.fft, ".", ms=3, color="tab:blue" if row == 0 else "tab:green")
        r = stats.linregress(g[col], g.fft)
        xx = np.array([g[col].min(), g[col].max()])
        ax.plot(xx, r.intercept + r.slope * xx, "k-", lw=1)
        ax.set_xlabel(lab, fontsize=8)
        ax.set_title(f"{mag}$\\times$ series, $p = {ptext(r.pvalue)}$", fontsize=8)
        ax.tick_params(labelsize=7)
    axes[0, 0].set_ylabel("FFT pitch per tile [µm]", fontsize=8)
    axes[1, 0].set_ylabel("FFT pitch per tile [µm]", fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(RESP, "figR_time_position.pdf"))
plt.close(fig)
print("figures written")
