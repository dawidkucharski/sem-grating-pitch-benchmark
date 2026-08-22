#!/usr/bin/env python3
"""Audit part 2: FFT, agreement, WLI, parallelism, simulation, drift, power."""
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = "/Users/dawid/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/SEM_measurements/Calibration"
GB = os.path.join(BASE, "sem_grating_benchmark")
S = 0.01333594
MAGS = [400, 500, 750, 1000]
ok, warn = 0, 0

def check(label, computed, claimed, tol=0.06):
    global ok, warn
    passf = abs(computed - claimed) <= tol * max(abs(claimed), 1e-9)
    print(f"[{'PASS' if passf else 'FAIL'}] {label:60s} computed={computed} claimed={claimed}")
    ok += passf
    warn += (not passf)

per = {}
for m in MAGS:
    for d in ["SED", "BED-S"]:
        df = pd.read_csv(f"{GB}/per_image_pitch_mag{m}_{d}.csv")
        per[(m, d)] = df
allper = pd.concat(per.values(), ignore_index=True)
fft = pd.read_csv(f"{GB}/fft_subpixel_per_image.csv")
wli = pd.read_csv(f"{BASE}/WLI/wli_results.csv")
sim = pd.read_csv(f"{GB}/simulation_validation.csv")
bab = pd.read_csv(f"{GB}/bland_altman_bootstrap_summary.csv")

print("== J. FFT (cont.) ==")
snr = fft.groupby(["mag", "det"]).snr.median()
print(f"      SNR medians per mag/det: {snr.values.round(0).tolist()} -> range {snr.min():.0f}-{snr.max():.0f} (claimed 233-324)")
check("FFT SNR median max ~324", snr.max(), 324, 0.02)
check("FFT SNR median min ~233", snr.min(), 233, 0.02)
print(f"      min SNR = {fft.snr.min():.0f} (manuscript: >150)")
check("FFT min SNR >150", 1.0 if fft.snr.min() > 150 else 0.0, 1.0, 0.05)
k_int = (1280 * S / fft.quant).round()
print(f"      integer bin k distribution: {k_int.value_counts().to_dict()}")
check("FFT nearest-bin pitch (k=10 arithmetic)", 1280 * S / 10, 1.707, 0.005)
k_ref = (1280 * S / fft.pitch)
print(f"      refined k median={k_ref.median():.3f} (claimed ~10.4)")

print("== K. 1D-FFT agreement ==")
mrg = allper.dropna(subset=["pitch_1d_um", "pitch_fft_um"]).copy()
mrg["diff"] = mrg.pitch_1d_um - mrg.pitch_fft_um
check("agreement n", len(mrg), 1432, 0.01)
check("agreement bias", mrg["diff"].mean(), 0.0045, 0.15)
loa = mrg["diff"].mean() + np.array([-1.96, 1.96]) * mrg["diff"].std(ddof=1)
check("agreement LoA low", loa[0], -0.062, 0.15)
check("agreement LoA high", loa[1], 0.071, 0.15)

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
check("WLI means agree ~0.00001", abs(r30.pitch_fft.mean() - r36.pitch_fft.mean()), 0.00001, 0.5)
check("WLI sd ratio ~1.9", r36.pitch_fft.std(ddof=1) / r30.pitch_fft.std(ddof=1), 1.9, 0.2)
check("WLI uniformity 36x [%]", (r36.pitch_fft.max() - r36.pitch_fft.min()) / r36.pitch_fft.mean() * 100, 0.04, 0.3)
check("WLI coh pitch (example scan)", wli[wli.campaign == "30x"].coh_pitch.iloc[0], 1.6697, 0.01)

print("== M. parallelism ==")
par_tab = {400: (87.17, 0.818, 2.83), 500: (86.95, 0.798, 3.05),
           750: (87.52, 0.424, 2.48), 1000: (87.68, 0.272, 2.32)}
sigs = []
for m in MAGS:
    r = pd.read_csv(f"{GB}/parallelism_summary_mag{m}.csv").iloc[0]
    sigs.append(r.sigma_parallel_deg)
    check(f"parallelism alpha {m}x", r.dominant_angle_deg, par_tab[m][0], 0.02)
    check(f"parallelism sigma {m}x", r.sigma_parallel_deg, par_tab[m][1], 0.08)
    check(f"parallelism dAlpha {m}x", r.delta_alpha_deg, par_tab[m][2], 0.05)
print(f"      sigma range {min(sigs):.3f}-{max(sigs):.3f} (claimed <0.9 deg)")

print("== N. simulation ==")
for _, r in sim.iterrows():
    print(f"      true {r.true_pitch_um}: bias1D={r.bias_1d_um:+.5f}um ({r.bias_1d_pct:+.3f}%)  biasFFT={r.bias_fft_um:+.5f}um ({r.bias_fft_pct:+.3f}%)")
check("sim mean bias 1D [um]", sim.bias_1d_um.mean(), 0.0005, 0.4)
check("sim mean bias 1D [%]", sim.bias_1d_pct.mean(), 0.03, 0.4)
print(f"      FFT unrefined bias range: {sim.bias_fft_um.min():+.3f}..{sim.bias_fft_um.max():+.3f} um ({sim.bias_fft_pct.min():+.1f}%..{sim.bias_fft_pct.max():+.1f}%) (claimed -0.11..+0.08 um / -6%..+5%)")
check("sim FFT bias min", sim.bias_fft_um.min(), -0.107, 0.1)
check("sim FFT bias max", sim.bias_fft_um.max(), 0.079, 0.1)

print("== O. temporal drift ==")
maxrate, worst_p = 0, 1.0
for m in MAGS:
    for d in ["SED", "BED-S"]:
        df = per[(m, d)].dropna(subset=["pitch_1d_um"]).sort_values("image_path")
        if len(df) > 10:
            sl = stats.linregress(np.arange(len(df)), df.pitch_1d_um)
            maxrate = max(maxrate, abs(sl.slope))
            worst_p = min(worst_p, sl.pvalue)
print(f"      max |drift rate| = {maxrate:.2e} um/image; min p = {worst_p:.3f} (manuscript: all p>0.05)")
check("drift all p > 0.05", 1.0 if worst_p > 0.05 else 0.0, 1.0, 0.05)
print(f"      drift rate informational max={maxrate:.2e} um/image")

print("== P. power analysis ==")
sds, mdds = [], []
for _, r in bab.iterrows():
    n = r.n_pairs
    sd = r.sd_diff_um
    mdd = stats.t.ppf(0.9, df=n - 1) * sd / np.sqrt(n)
    sds.append(sd); mdds.append(mdd)
    print(f"      {int(r.magnification)}x: within-pair sd={sd:.3f} MDD={mdd:.3f}")
print(f"      sd range {min(sds):.3f}-{max(sds):.3f} (claimed 0.036-0.063); MDD range {min(mdds):.3f}-{max(mdds):.3f} (claimed 0.006-0.015)")

print(f"\n=== part2: {ok} passed, {warn} failed ===")
