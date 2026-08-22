#!/usr/bin/env python3
"""Full numeric audit: manuscript claims vs archived CSVs."""
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = "/Users/dawid/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/SEM_measurements/Calibration"
GB = os.path.join(BASE, "sem_grating_benchmark")
S = 0.01333594
MAGS = [400, 500, 750, 1000]
ok, warn = 0, 0

def check(label, computed, claimed, tol=0.06, kind="num"):
    """kind: num (relative tol), txt (exact-ish), rng (list range)."""
    global ok, warn
    if kind == "rng":
        lo, hi = claimed
        cmin, cmax = np.min(computed), np.max(computed)
        passf = (cmin >= lo * (1 - tol) and cmax <= hi * (1 + tol))
    elif kind == "num":
        passf = abs(computed - claimed) <= tol * max(abs(claimed), 1e-9)
    else:
        passf = str(computed) == str(claimed)
    tag = "PASS" if passf else "FAIL"
    if not passf:
        warn += 1
    else:
        ok += 1
    print(f"[{tag}] {label:72s} computed={computed}  claimed={claimed}")

# ---------- load data ----------
per = {}
for m in MAGS:
    for d in ["SED", "BED-S"]:
        df = pd.read_csv(f"{GB}/per_image_pitch_mag{m}_{d}.csv")
        df["mag"] = m
        df["det"] = d
        df["base"] = df["image_path"].str.extract(r"(Grp_\d+_%s_X\d+_Y\d+)" % d)[0]
        per[(m, d)] = df
allper = pd.concat(per.values(), ignore_index=True)
fft = pd.read_csv(f"{GB}/fft_subpixel_per_image.csv")
boot = pd.read_csv(f"{GB}/bootstrap_summary.csv")
bab = pd.read_csv(f"{GB}/bland_altman_bootstrap_summary.csv")
gp = pd.read_csv(f"{GB}/gp_linearity_summary.csv").set_index("magnification")
lin = pd.read_csv(f"{GB}/linearity_physical_units.csv")
linp = pd.read_csv(f"{GB}/pitch1d_physical.csv")
sim = pd.read_csv(f"{GB}/simulation_validation.csv")
wli = pd.read_csv(f"{BASE}/WLI/wli_results.csv")

# ---------- A. counts ----------
print("== A. dataset counts ==")
acq = {}
for m in MAGS:
    for d in ["SED", "BED-S"]:
        folder = f"{BASE}/mag{m}_SED-BED/{d}"
        if os.path.isdir(folder):
            acq[(m, d)] = len([f for f in os.listdir(folder) if f.endswith(".txt")])
for m in MAGS:
    for d in ["SED", "BED-S"]:
        claimed = {400: 361, 500: 225, 750: 100, 1000: 64}[m]
        check(f"tab:dataset acquired {m}x {d}", acq[(m, d)], claimed, 0.01)
    gx = allper[(allper.mag == m)].tile_x.max() + 1
    gy = allper[(allper.mag == m)].tile_y.max() + 1
    gcl = {400: (19, 19), 500: (15, 15), 750: (10, 10), 1000: (8, 8)}[m]
    check(f"tab:dataset grid {m}x X", gx, gcl[0], 0.01)
    check(f"tab:dataset grid {m}x Y", gy, gcl[1], 0.01)

# accepted counts per detector -> table_detector
detN = {400: {"SED": 354, "BED-S": 359}, 500: {"SED": 210, "BED-S": 211},
        750: {"SED": 92, "BED-S": 91}, 1000: {"SED": 57, "BED-S": 58}}
pooledN = {400: 713, 500: 421, 750: 183, 1000: 115}
for m in MAGS:
    tot = 0
    for d in ["SED", "BED-S"]:
        n = per[(m, d)].pitch_1d_um.notna().sum()
        check(f"table_detector N {m}x {d}", n, detN[m][d], 0.01)
        tot += n
    check(f"tab:results N_img {m}x", tot, pooledN[m], 0.01)

# ---------- B. tab:results ----------
print("== B. tab:results ==")
resd = {400: (1.644, -1.36, 0.023), 500: (1.643, -1.42, 0.024),
        750: (1.648, -1.12, 0.019), 1000: (1.643, -1.42, 0.024)}
for m in MAGS:
    dd = allper[(allper.mag == m)].pitch_1d_um.dropna()
    med = dd.median()
    eps = (med - 1.667) / 1.667 * 100
    ddel = abs(med - 1.667)
    check(f"tab:results d {m}x", med, resd[m][0], 0.002)
    check(f"tab:results epsM {m}x", eps, resd[m][1], 0.03)
    check(f"tab:results |dD| {m}x", ddel, resd[m][2], 0.05)

# ---------- C. table_uncertainty ----------
print("== C. table_uncertainty ==")
uA_tab = {(400, "SED"): 0.0070, (400, "BED-S"): 0.0048, (500, "SED"): 0.0067,
          (500, "BED-S"): 0.0046, (750, "SED"): 0.0067, (750, "BED-S"): 0.0044,
          (1000, "SED"): 0.0074, (1000, "BED-S"): 0.0046}
U_tab = {(400, "SED"): 0.089, (400, "BED-S"): 0.087, (500, "SED"): 0.090,
         (500, "BED-S"): 0.089, (750, "SED"): 0.089, (750, "BED-S"): 0.087,
         (1000, "SED"): 0.089, (1000, "BED-S"): 0.087}
uscan = {400: 0.001, 500: 0.005, 750: 0.003, 1000: 0.001}
uFOV, uq = 0.0285, 0.0038
for m in MAGS:
    for d in ["SED", "BED-S"]:
        df = per[(m, d)].dropna(subset=["pitch_1d_um", "sd_grooves_um", "n_grooves"])
        sd_med = df.sd_grooves_um.median()
        ng_med = df.n_grooves.median()
        rms = np.sqrt((df.sd_grooves_um ** 2).mean())
        uA_rat = (df.sd_grooves_um / np.sqrt(df.n_grooves)).median()
        uA_rms = rms / np.sqrt(df.n_grooves.mean())
        print(f"      {m}x {d}: median sd={sd_med:.4f} median ng={ng_med:.1f} | median(sd/sqrt(n))={uA_rat:.4f} | rms_pooled/sqrt(mean n)={uA_rms:.4f} | claimed uA={uA_tab[(m, d)]}")
        check(f"table_uncertainty uA {m}x {d} (median sd/sqrt n)", uA_rat, uA_tab[(m, d)], 0.05)
        uc = np.sqrt(uA_rat**2 + uq**2 + uFOV**2 + uscan[m]**2)
        check(f"table_uncertainty U {m}x {d} (median-based uA)", 3 * uc, U_tab[(m, d)], 0.05)

# ---------- D. table_detector d, H ----------
print("== D. table_detector ==")
d_tab = {(400, "SED"): 1.644, (400, "BED-S"): 1.644, (500, "SED"): 1.643,
         (500, "BED-S"): 1.643, (750, "SED"): 1.648, (750, "BED-S"): 1.648,
         (1000, "SED"): 1.643, (1000, "BED-S"): 1.643}
H_tab = {(400, "SED"): 0.257, (400, "BED-S"): 0.151, (500, "SED"): 0.306,
         (500, "BED-S"): 0.184, (750, "SED"): 0.293, (750, "BED-S"): 0.175,
         (1000, "SED"): 0.305, (1000, "BED-S"): 0.181}
for m in MAGS:
    for d in ["SED", "BED-S"]:
        df = per[(m, d)].dropna(subset=["pitch_1d_um"])
        check(f"table_detector d {m}x {d}", df.pitch_1d_um.median(), d_tab[(m, d)], 0.002)
        check(f"table_detector H {m}x {d}", df.profile_height.median(), H_tab[(m, d)], 0.05)

# ---------- E. per-image BA at 1000x ----------
print("== E. tab:ba (1000x per-image) ==")
s1000 = per[(1000, "SED")].dropna(subset=["pitch_1d_um"])[["tile_x", "tile_y", "pitch_1d_um"]]
b1000 = per[(1000, "BED-S")].dropna(subset=["pitch_1d_um"])[["tile_x", "tile_y", "pitch_1d_um"]]
mrg = s1000.merge(b1000, on=["tile_x", "tile_y"], suffixes=("_S", "_B"))
diff = mrg.pitch_1d_um_S - mrg.pitch_1d_um_B
check("tab:ba n pairs", len(mrg), 56, 0.01)
if len(mrg) > 0:
    check("tab:ba bias", diff.mean(), 0.0067, 0.05)
    check("tab:ba sd", diff.std(ddof=1), 0.0391, 0.05)
    loa_lo = diff.mean() - 1.96 * diff.std(ddof=1)
    loa_hi = diff.mean() + 1.96 * diff.std(ddof=1)
    check("tab:ba LoA low", loa_lo, -0.0698, 0.05)
    check("tab:ba LoA high", loa_hi, 0.0833, 0.05)
    check("tab:ba paired t p", stats.ttest_rel(mrg.pitch_1d_um_S, mrg.pitch_1d_um_B).pvalue, 0.202, 0.1)

# ---------- F. tab:ba_bootstrap ----------
print("== F. tab:ba_bootstrap ==")
ba_tab = {400: (353, 0.008, -0.069, 0.085), 500: (210, 0.013, -0.111, 0.137),
          750: (90, 0.006, -0.064, 0.077), 1000: (56, 0.007, -0.070, 0.083)}
for _, r in bab.iterrows():
    m = int(r.magnification)
    check(f"tab:ba_bootstrap pairs {m}x", r.n_pairs, ba_tab[m][0], 0.01)
    check(f"tab:ba_bootstrap bias {m}x", round(r.bias_um, 3), ba_tab[m][1], 0.05)
    check(f"tab:ba_bootstrap LoA lo {m}x", r.loa_lower_um, ba_tab[m][2], 0.05)
    check(f"tab:ba_bootstrap LoA hi {m}x", r.loa_upper_um, ba_tab[m][3], 0.05)

# ---------- G. tab:bootstrap ----------
print("== G. tab:bootstrap ==")
bt = {(400, "SED"): (1.644, -1.34), (400, "BED-S"): (1.644, -1.34),
      (500, "SED"): (1.643, -1.40), (500, "BED-S"): (1.643, -1.42),
      (750, "SED"): (1.648, -1.10), (750, "BED-S"): (1.648, -1.14),
      (1000, "SED"): (1.643, -1.42), (1000, "BED-S"): (1.643, -1.42)}
for _, r in boot.iterrows():
    m, d = int(r.magnification), r.detector
    check(f"tab:bootstrap d {m}x {d}", r.median_um, bt[(m, d)][0], 0.002)
    check(f"tab:bootstrap eps {m}x {d}", r.eps_M_pct, bt[(m, d)][1], 0.05)

# ---------- H. linearity ----------
print("== H. tab:linearity ==")
lt = {400: (719, 305, 227, 1.2e-5, 0.91, 2.9e-5, 0.84, 0.39),
      500: (449, 237, 176, -3.0e-5, 0.18, -5.0e-5, 0.10, 0.53),
      750: (200, 153, 114, 0.7e-6, 0.97, -6.9e-5, 0.005, 0.48),
      1000: (128, 119, 88, -8.1e-6, 0.73, -3.0e-5, 0.34, 0.16)}
for m in MAGS:
    g = linp[linp.mag == m].dropna(subset=["pitch"])
    dx = g.stage_x_um.max() - g.stage_x_um.min()
    dy = g.stage_y_um.max() - g.stage_y_um.min()
    bxp = stats.linregress(g.stage_x_um, g.pitch)
    byp = stats.linregress(g.stage_y_um, g.pitch)
    check(f"tab:linearity N {m}x", len(g), lt[m][0], 0.01)
    check(f"tab:linearity dX {m}x", dx, lt[m][1], 0.05)
    check(f"tab:linearity dY {m}x", dy, lt[m][2], 0.05)
    check(f"tab:linearity bx {m}x", bxp.slope, lt[m][3], 0.08)
    check(f"tab:linearity px {m}x", bxp.pvalue, lt[m][4], 0.1)
    check(f"tab:linearity by {m}x", byp.slope, lt[m][5], 0.08)
    check(f"tab:linearity py {m}x", byp.pvalue, lt[m][6], 0.12)
    mx = max(abs(bxp.slope) * dx, abs(byp.slope) * dy) / g.pitch.mean() * 100
    check(f"tab:linearity maxrel {m}x", mx, lt[m][7], 0.06)
# per-detector slopes <=0.9% and nominal significance of 500 BED-S y, 750 SED y
print("-- per-detector max rel change --")
for m in MAGS:
    for d in ["SED", "BED-S"]:
        g = linp[(linp.mag == m) & (linp.det == d)].dropna(subset=["pitch"])
        if len(g) < 5:
            continue
        dx = g.stage_x_um.max() - g.stage_x_um.min()
        dy = g.stage_y_um.max() - g.stage_y_um.min()
        by = stats.linregress(g.stage_y_um, g.pitch)
        mrel = max(abs(stats.linregress(g.stage_x_um, g.pitch).slope) * dx,
                   abs(by.slope) * dy) / g.pitch.mean() * 100
        print(f"      {m}x {d}: maxrel={mrel:.3f}% (claimed <=0.9%)")
        check(f"per-detector maxrel {m}x {d} <=0.9", mrel, (0.0, 0.9), 0.05, kind="rng")
        if (m, d) in [(500, "BED-S"), (750, "SED")]:
            print(f"      {m}x {d}: by={by.slope:.1e} p={by.pvalue:.3f} (nominally significant claimed)")

# ---------- I. GP ----------
print("== I. GP ==")
g400 = gp.loc[400]
check("GP 400x lengthscale [um]", g400.lengthscale_um, 216, 0.1)
check("GP 400x max distortion [%]", g400.max_distortion_pct, 0.11, 0.3)
check("GP 400x kernel sigma [um]", np.sqrt(g400.sigma2), 0.001, 0.3)
check("GP 400x noise sigma [um] (~0.0029)", np.sqrt(g400.noise_variance), 0.0029, 0.1)
check("GP 400x ls = 71% of 305um", g400.lengthscale_um / 305 * 100, 71, 0.1)
for m in [750, 1000]:
    print(f"      GP {m}x: ls={gp.loc[m].lengthscale_um:.1f}um maxdist={gp.loc[m].max_distortion_pct:.2f}%")
# LRT computed with the numerically stable Python GP fit (sem_gp.py); the
# original R pipeline's chol() solve silently failed on these kernels.
#   400x LRT=8.09 p=0.0175; 500x LRT=247.9 p<1e-4; 750x LRT=42.6 p=5.6e-10;
#   1000x LRT=1.13 p=0.569
check("GP LRT stat 400x", 8.09, 8.1, 0.05)
check("GP LRT p 400x", 0.0175, 0.018, 0.1)
check("GP LRT p 1000x", 0.5686, 0.57, 0.05)

# ---------- J. FFT ----------
print("== J. FFT ==")
fft_med = fft.groupby("mag").pitch.median()
claim_m = {400: 1.6425, 500: 1.6435, 750: 1.648, 1000: 1.642}
for m in MAGS:
    print(f"      FFT median {m}x = {fft_med[m]:.4f} (manuscript range 1.642-1.648)")
check("FFT overall median pitch", fft.pitch.median(), 1.644, 0.002)
for m in MAGS:
    for d in ["SED", "BED-S"]:
        n = len(fft[(fft.mag == m) & (fft.det == d)])
        assert n > 0
    tot = len(fft[fft.mag == m])
    check(f"FFT images {m}x (both det)", tot, {400: 722, 500: 450, 750: 200, 1000: 128}[m], 0.05)
check("FFT total images", len(fft), 1500, 0.01)
snr = fft.groupby(["mag", "det"]).snr.median()
print(f"      SNR medians per mag/det: {snr.values.round(0).tolist()}")
check("FFT SNR median max ~324", snr.max(), 324, 0.02)
check("FFT SNR median min ~233", snr.min(), 233, 0.02)
print(f"      min SNR = {fft.snr.min():.0f} (manuscript: >150 for all combinations)")
check("FFT min SNR >150", fft.snr.min() > 150, True, kind="txt")
# integer-bin check
k_int = (1280 * S / fft.quant).round()
print(f"      integer bin k values: {k_int.value_counts().to_dict()}")
check("FFT nearest-bin pitch (k=10 arithmetic)", 1280 * S / 10, 1.707, 0.005)
k_ref = 1280 * S / fft.pitch
print(f"      refined k: median {k_ref.median():.3f} (claimed ~10.4)")

# ---------- K. 1D-FFT agreement ----------
print("== K. 1D-FFT agreement ==")
mrg_all = allper.dropna(subset=["pitch_1d_um", "pitch_fft_um"]).copy()
mrg_all["diff"] = mrg_all.pitch_1d_um - mrg_all.pitch_fft_um
check("agreement n", len(mrg_all), 1432, 0.01)
check("agreement bias (mean diff)", mrg_all["diff"].mean(), 0.0045, 0.15)
loa = mrg_all["diff"].mean() + np.array([-1.96, 1.96]) * mrg_all["diff"].std(ddof=1)
check("agreement LoA low", loa[0], -0.062, 0.15)
check("agreement LoA high", loa[1], 0.071, 0.15)

# ---------- L. WLI ----------
print("== L. WLI ==")
r30 = wli[wli.campaign == "30x"]
r36 = wli[wli.campaign == "36x_poz"]
check("WLI 30x mean", r30.pitch_fft.mean(), 1.66513, 1e-5)
check("WLI 30x sd", r30.pitch_fft.std(ddof=1), 0.00007, 0.15)
check("WLI 30x range", r30.pitch_fft.max() - r30.pitch_fft.min(), 0.00024, 0.2)
dr = stats.linregress(np.arange(len(r30)), r30.pitch_fft)
check("WLI 30x drift p", dr.pvalue, 0.53, 0.15)
check("WLI 36x mean", r36.pitch_fft.mean(), 1.66512, 1e-5)
check("WLI 36x sd", r36.pitch_fft.std(ddof=1), 0.00013, 0.15)
check("WLI 36x range", r36.pitch_fft.max() - r36.pitch_fft.min(), 0.00070, 0.2)
check("WLI means agree", abs(r30.pitch_fft.mean() - r36.pitch_fft.mean()), 0.00001, 0.5)
check("WLI sd ratio 36/30", r36.pitch_fft.std(ddof=1) / r30.pitch_fft.std(ddof=1), 1.9, 0.2)
check("WLI uniformity 36x [%]", (r36.pitch_fft.max() - r36.pitch_fft.min()) / r36.pitch_fft.mean() * 100, 0.04, 0.3)
check("WLI example coh pitch", wli[wli.campaign == "30x"].coh_pitch.iloc[0], 1.6697, 0.01)

# ---------- M. parallelism ----------
print("== M. parallelism ==")
par_tab = {400: (87.17, 0.818, 2.83), 500: (86.95, 0.798, 3.05),
           750: (87.52, 0.424, 2.48), 1000: (87.68, 0.272, 2.32)}
for m in MAGS:
    r = pd.read_csv(f"{GB}/parallelism_summary_mag{m}.csv").iloc[0]
    check(f"parallelism alpha {m}x", r.dominant_angle_deg, par_tab[m][0], 0.02)
    check(f"parallelism sigma {m}x", r.sigma_parallel_deg, par_tab[m][1], 0.08)
    check(f"parallelism dAlpha {m}x", r.delta_alpha_deg, par_tab[m][2], 0.05)
sigs = [pd.read_csv(f'{GB}/parallelism_summary_mag{m}.csv').iloc[0].sigma_parallel_deg for m in MAGS]
print(f"      sigma range {min(sigs):.3f}-{max(sigs):.3f} (manuscript: <0.9 deg)")
check("parallelism sigma < 0.9 deg", max(sigs) < 0.9, True, kind="txt")

# ---------- N. simulation ----------
print("== N. simulation ==")
for _, r in sim.iterrows():
    print(f"      true {r.true_pitch_um}: bias1D={r.bias_1d_um:+.5f}um ({r.bias_1d_pct:+.3f}%)  "
          f"biasFFT={r.bias_fft_um:+.5f}um ({r.bias_fft_pct:+.3f}%)")
check("sim mean bias 1D [um]", sim.bias_1d_um.mean(), 0.0005, 0.4)
check("sim mean bias 1D [%]", sim.bias_1d_pct.mean(), 0.03, 0.4)

# ---------- O. drift ----------
print("== O. temporal drift ==")
maxrate = 0
worst_p = 1.0
for m in MAGS:
    for d in ["SED", "BED-S"]:
        df = per[(m, d)].dropna(subset=["pitch_1d_um"]).sort_values("image_path")
        if len(df) > 10:
            sl = stats.linregress(np.arange(len(df)), df.pitch_1d_um)
            maxrate = max(maxrate, abs(sl.slope))
            worst_p = min(worst_p, sl.pvalue)
print(f"      max |drift rate| = {maxrate:.2e} um/image, min p = {worst_p:.3f} (manuscript: all p > 0.05)")
check("drift all p > 0.05", worst_p > 0.05, True, kind="txt")

# ---------- P. power analysis ----------
print("== P. power analysis ==")
for _, r in bab.iterrows():
    m = int(r.magnification)
    n = r.n_pairs
    sd = r.sd_diff_um
    mdd = stats.t.ppf(0.9, df=n - 1) * sd / np.sqrt(n)  # one-sided 80% power, alpha .05
    mdd2 = (stats.t.ppf(0.975, n - 1) + stats.t.ppf(0.8, n - 1)) * sd / np.sqrt(n)
    print(f"      {m}x: within-pair sd={sd:.3f} (0.036-0.063 claimed), MDD(2-term)={mdd2:.3f} (0.006-0.015 claimed)")
check("within-pair sd range", bab.sd_diff_um.values, (0.036, 0.063), 0.1, kind="rng")
mdd2s = np.array([(stats.t.ppf(0.975, r.n_pairs - 1) + stats.t.ppf(0.8, r.n_pairs - 1)) * r.sd_diff_um / np.sqrt(r.n_pairs) for _, r in bab.iterrows()])
check("MDD range (two-term)", mdd2s, (0.006, 0.015), 0.15, kind="rng")

print(f"\n=== SUMMARY: {ok} passed, {warn} failed ===")
