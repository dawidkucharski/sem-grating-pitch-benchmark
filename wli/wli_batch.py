#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch analysis of the Optiv WLI campaign (multisensor CMM, Hexagon Optiv
DualZ 763) on the holographic grating (Thorlabs GH13-06U, 600 lines/mm).

Inputs:
  WLI/30x/1.txt ... 30.txt     -> 30 scans at ONE position (repeatability)
  WLI/36x_poz/1.txt ... 36.txt -> 36 scans at DIFFERENT positions (spatial)

Per-scan quantities (same methodology as wli_pitch.py):
  - grid reconstruction (1904x1200 px, ~0.4745 um/px), missing-point fill
  - least-squares plane removal
  - 2D FFT: dominant pitch, groove angle, peak SNR
  - 1D SEM-style procedure on the best-aligned 256x256 px block:
    B-spline profile, minima detection, regular-period filter
    (keep 0.7x-1.3x of the median -> scratches/defects rejected)

Outputs: WLI/wli_results.csv, WLI/WLI_repeatability_30x.pdf,
         WLI/WLI_spatial_map_36x.pdf (+ PNG copies)
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy.ndimage import rotate, convolve
from scipy import stats

BASE = os.path.dirname(os.path.abspath(__file__))
FOLDERS = {"30x": os.path.join(BASE, "30x"),
           "36x_poz": os.path.join(BASE, "36x_poz")}
CSV_OUT = os.path.join(BASE, "wli_results.csv")


def load_scan(path):
    """Load xyz file (mm), rebuild regular grid, fill missing points."""
    df = pd.read_csv(path, sep=r"\s+", header=None,
                     names=["x", "y", "z"], engine="c")
    x, y, z = df["x"].values, df["y"].values, df["z"].values
    xs, ys = np.unique(x), np.unique(y)
    nx, ny = xs.size, ys.size
    col = np.searchsorted(xs, x)
    row = np.searchsorted(ys, y)
    Z = np.full((ny, nx), np.nan)
    Z[row, col] = z
    n_nan = int(np.isnan(Z).sum())
    if n_nan:
        # wiersze sa postrzepione (braki do ~75 px): interpolacja liniowa
        # wzdluz kazdego wiersza; braki na koncach wiersza -> wartosc skrajna
        for i in range(ny):
            r = Z[i]
            nan = np.isnan(r)
            if not nan.any():
                continue
            cols = np.arange(nx)
            valid = ~nan
            nv = int(valid.sum())
            if nv == 0:
                Z[i, :] = np.nanmedian(Z)
                continue
            r[nan] = np.interp(cols[nan], cols[valid], r[valid])
            if nv == 1:
                Z[i, :] = r[valid][0]
                continue
            still = np.isnan(r)
            if still.any():
                first = int(np.argmax(valid))
                last = int(nx - 1 - np.argmax(valid[::-1]))
                Z[i, :first] = r[first]
                Z[i, last + 1:] = r[last]
    if np.isnan(Z).any():
        Z = np.where(np.isnan(Z), np.nanmedian(Z), Z)
    dx = (xs[-1] - xs[0]) / (nx - 1) * 1000.0     # um/px
    dy = (ys[-1] - ys[0]) / (ny - 1) * 1000.0
    return Z, dx, dy, n_nan, float(xs[0]), float(ys[0])


def find_peaks_manual(y, frac=0.02):
    rng = np.ptp(y)
    thr = frac * rng
    d = np.diff(y)
    p = np.where((d[:-1] > 0) & (d[1:] < 0))[0] + 1
    return p[y[p] >= thr]


def analyse_scan(path):
    Z, dx, dy, n_nan, x0, y0 = load_scan(path)
    ny, nx = Z.shape
    xx = np.arange(nx) * dx
    yy = np.arange(ny) * dy
    Xg, Yg = np.meshgrid(xx, yy)
    A = np.column_stack([Xg.ravel(), Yg.ravel(), np.ones(Xg.size)])
    coef, *_ = np.linalg.lstsq(A, Z.ravel(), rcond=None)
    Zl = Z - (A @ coef).reshape(ny, nx)
    ptv_um = (Zl.max() - Zl.min()) * 1000.0

    # ---- 2D FFT (nieprzesuniete widmo + fftfreq: jednoznaczna konwersja) ----
    S = Zl - Zl.mean()
    F2 = np.fft.fft2(S)
    Pa = np.abs(F2)
    fyv = np.fft.fftfreq(ny, d=dy)
    fxv = np.fft.fftfreq(nx, d=dx)
    fg = np.sqrt(fyv[:, None] ** 2 + fxv[None, :] ** 2)

    # pasmo siatki: 1.52-1.82 um (0.55-0.66 cykli/um)
    ann = (fg >= 0.55) & (fg <= 0.66)
    kp = np.unravel_index(np.argmax(Pa * ann), Pa.shape)
    coh_pitch = 1.0 / fg[kp] if fg[kp] > 0 else np.nan
    coh_amp_nm = 2.0 * Pa[kp] / (nx * ny) * 1000.0
    # pasmo artefaktu linii skanowania (1.00-1.10 cykli/um) - diagnostyka
    art = (fg >= 1.00) & (fg <= 1.10)
    art_amp_nm = 2.0 * Pa[art].max() / (nx * ny) * 1000.0

    # ---- estymacja blokowa: pitch z blokow 256x256 (stride 128) ----
    # faza rowkow wedruje w calym polu (dystorsja optyczna), przez co
    # koherentny pik calego pola jest rozmyty; lokalne bloki daja rzetelny
    # pitch -> mediana po blokach
    bs = 256
    fby = np.fft.fftfreq(bs, d=dy)
    fbx = np.fft.fftfreq(bs, d=dx)
    fbb = np.sqrt(fby[:, None] ** 2 + fbx[None, :] ** 2)
    mbb = (fbb >= 0.55) & (fbb <= 0.66)
    blk_pitch, blk_amp, blk_snr, blk_ang = [], [], [], []
    best = (-1.0, 0, 0, 0.0)
    for i0 in range(0, ny - bs + 1, 128):
        for j0 in range(0, nx - bs + 1, 128):
            blk = Zl[i0:i0 + bs, j0:j0 + bs]
            Fb = np.fft.fft2(blk - blk.mean())
            Pb = np.abs(Fb)
            kpb = np.unravel_index(np.argmax(Pb * mbb), Pb.shape)
            iyb, ixb = kpb
            # sub-binowe doprecyzowanie (interpolacja paraboliczna 3-pkt)
            def refine(P, a, b, axis):
                if axis == 0:
                    y = [P[max(a - 1, 0), b], P[a, b], P[min(a + 1, bs - 1), b]]
                else:
                    y = [P[a, max(b - 1, 0)], P[a, b], P[a, min(b + 1, bs - 1)]]
                den = y[0] - 2 * y[1] + y[2]
                if den == 0:
                    return 0.0
                return float(np.clip(0.5 * (y[0] - y[2]) / den, -0.5, 0.5))
            dyb = refine(Pb, iyb, ixb, 0)
            dxb = refine(Pb, iyb, ixb, 1)
            ky_r = iyb + dyb
            kx_r = ixb + dxb
            if ky_r > bs // 2:
                ky_r -= bs
            if kx_r > bs // 2:
                kx_r -= bs
            fy_r = ky_r / (bs * dy)
            fx_r = kx_r / (bs * dx)
            fm = np.hypot(fx_r, fy_r)
            amp = 2.0 * Pb[kpb] / (bs * bs) * 1000.0
            medb = 2.0 * np.median(Pb[mbb]) / (bs * bs) * 1000.0
            snrb = amp / medb if medb > 0 else 0.0
            if amp < 0.05 or fm <= 0:
                continue
            blk_pitch.append(1.0 / fm)
            blk_amp.append(amp)
            blk_snr.append(snrb)
            blk_ang.append(np.degrees(np.arctan2(fy_r, fx_r)))
            if amp > best[0]:
                best = (amp, i0, j0, np.degrees(np.arctan2(fy_r, fx_r)))
    blk_pitch = np.array(blk_pitch)
    blk_amp = np.array(blk_amp)
    blk_snr = np.array(blk_snr)
    blk_ang = np.array(blk_ang)
    if blk_pitch.size:
        pitch_fft = float(np.median(blk_pitch))
        pitch_sd = float(np.std(blk_pitch, ddof=1)) if blk_pitch.size > 1 else 0.0
        angle = float(np.median(blk_ang))
        snr = float(np.median(blk_snr))
        amp_best_nm = float(np.median(blk_amp))
        n_blk = int(blk_pitch.size)
    else:
        pitch_fft = pitch_sd = angle = snr = amp_best_nm = np.nan
        n_blk = 0

    # ---- 1D SEM-style na najlepszym bloku ----
    amp_best_nm_b, i0, j0, groove_angle_blk = best
    mean_1d = med_1d = reg_1d = np.nan
    n_drop = n_min = 0
    if amp_best_nm_b > 0:
        tile = Zl[i0:i0 + bs, j0:j0 + bs]
        tiler = rotate(tile, angle=-(groove_angle_blk + 90.0),
                       reshape=False, order=1)
        tiler = tiler[16:-16, 16:-16]
        prof = tiler.mean(axis=1)
        pos = np.arange(prof.size)
        spl = UnivariateSpline(pos, prof,
                               s=0.02 * np.sum((prof - prof.mean()) ** 2))
        sp = spl(pos)
        minima = find_peaks_manual(-sp)
        if minima.size > 1:
            n_min = int(minima.size)
            dist = np.diff(minima) * dy
            med = np.median(dist)
            keep = (dist >= 0.7 * med) & (dist <= 1.3 * med)
            reg = dist[keep]
            n_drop = int((~keep).sum())
            if reg.size:
                mean_1d = float(reg.mean())
                med_1d = float(np.median(reg))
            kidx = np.unique(np.concatenate([np.where(keep)[0],
                                            np.where(keep)[0] + 1]))
            if kidx.size > 5:
                km = minima[kidx]
                c = np.polyfit(np.arange(km.size), km, 1)
                reg_1d = float(c[0] * dy)
    return dict(file=os.path.basename(path), nx=nx, ny=ny, dx=dx, dy=dy,
                x0_mm=x0, y0_mm=y0, ptv_um=ptv_um, pitch_fft=pitch_fft,
                pitch_sd=pitch_sd, angle_deg=angle, snr=snr,
                amp_best_nm=amp_best_nm, n_blk=n_blk,
                coh_pitch=coh_pitch, coh_amp_nm=coh_amp_nm,
                art_amp_nm=art_amp_nm,
                mean_1d=mean_1d, med_1d=med_1d, reg_1d=reg_1d,
                n_drop=n_drop, n_min=n_min, n_nan=n_nan,
                status="accepted" if (n_blk >= 5 and snr >= 3)
                else "rejected_low_snr")


def main():
    rows = []
    for campaign, folder in FOLDERS.items():
        files = sorted([f for f in os.listdir(folder) if f.endswith(".txt")],
                       key=lambda s: int(s[:-4]))
        print(f"[{campaign}] {len(files)} files")
        for i, f in enumerate(files, 1):
            r = analyse_scan(os.path.join(folder, f))
            r["campaign"] = campaign
            r["file_no"] = i
            rows.append(r)
            print(f"  {f:>6}: FFT pitch {r['pitch_fft']:.4f} um | "
                  f"1D reg {r['reg_1d']:.4f} um | SNR {r['snr']:.0f} | "
                  f"angle {r['angle_deg']:.2f}°")
    df = pd.DataFrame(rows)
    df.to_csv(CSV_OUT, index=False)
    print(f"\nSaved: {CSV_OUT}")

    rep = df[df.campaign == "30x"]
    spa = df[df.campaign == "36x_poz"]
    rep_acc = rep[rep.status == "accepted"]
    spa_acc = spa[spa.status == "accepted"]

    print(f"\n=== classification (SNR >= 20) ===")
    print(f"30x:      accepted {len(rep_acc)}/{len(rep)}")
    print(f"36x_poz:  accepted {len(spa_acc)}/{len(spa)}")

    # ---- summary statistics ----
    def stats_summary(sub, label):
        if len(sub) == 0:
            print(f"\n=== {label} (n=0) — no accepted scans ===")
            return
        p = sub.pitch_fft
        print(f"\n=== {label} (n={len(p)}) — FFT pitch [um] ===")
        print(f"  mean  {p.mean():.5f}   sd {p.std(ddof=1):.5f}   "
              f"min {p.min():.5f}   max {p.max():.5f}")
        print(f"  median {p.median():.5f}   MAD {stats.median_abs_deviation(p):.5f}")
        print(f"  u_mean(sd/sqrt n) {p.std(ddof=1)/np.sqrt(len(p)):.5f}")
        r1 = sub.reg_1d.dropna()
        print(f"  1D reg-period: mean {r1.mean():.5f} sd {r1.std(ddof=1):.5f} (n={len(r1)})")
        r1m = sub.mean_1d.dropna()
        print(f"  1D regular mean: mean {r1m.mean():.5f} sd {r1m.std(ddof=1):.5f} (n={len(r1m)})")
        drift = stats.linregress(sub.file_no, p)
        print(f"  drift vs scan order: slope {drift.slope*1e6:.3f} nm/scan, "
              f"p = {drift.pvalue:.3f}")

    stats_summary(rep_acc, "30x repeatability (same position, accepted)")
    stats_summary(spa_acc, "36x_poz (different positions, accepted)")

    # ---- figures ----
    # 1) repeatability
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    ax = axes[0]
    p = rep_acc.pitch_fft.values
    ax.plot(rep_acc.file_no, p, "o-", color="tab:blue", label="accepted")
    rej = rep[rep.status != "accepted"]
    if len(rej):
        ax.plot(rej.file_no, np.full(len(rej), np.nan), "x", color="grey",
                ms=8, mew=1.5, label="rejected (low SNR)")
    ax.axhline(p.mean(), color="tab:red", lw=1.2, label=f"mean = {p.mean():.5f} µm")
    ax.axhspan(p.mean() - p.std(ddof=1), p.mean() + p.std(ddof=1),
               color="tab:red", alpha=0.12, label="±1 SD")
    ax.set_xlabel("scan number")
    ax.set_ylabel("pitch (2D FFT) [µm]")
    ax.set_title(f"WLI repeatability, same position ({len(rep_acc)}/{len(rep)} accepted)")
    ax.legend()
    ax.grid(ls=":", alpha=0.5)
    ax = axes[1]
    r1 = rep_acc.reg_1d.values
    ax.plot(rep_acc.file_no, r1, "s-", color="tab:green")
    ax.axhline(np.nanmean(r1), color="tab:red", lw=1.2,
               label=f"mean = {np.nanmean(r1):.5f} µm")
    ax.set_xlabel("scan number")
    ax.set_ylabel("pitch (1D regular periods) [µm]")
    ax.set_title("WLI repeatability, 1D procedure (accepted)")
    ax.legend()
    ax.grid(ls=":", alpha=0.5)
    fig.tight_layout()
    f1 = os.path.join(BASE, "WLI_repeatability_30x.pdf")
    fig.savefig(f1, dpi=200)
    fig.savefig(f1.replace(".pdf", ".png"), dpi=150)
    print(f"\nSaved: {f1}")

    # 2) spatial map
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    ax = axes[0]
    sc = ax.scatter(spa_acc.x0_mm, spa_acc.y0_mm, c=spa_acc.pitch_fft, s=55,
                    cmap="viridis", edgecolors="k", linewidths=0.4)
    rej = spa[spa.status != "accepted"]
    if len(rej):
        ax.scatter(rej.x0_mm, rej.y0_mm, marker="x", s=60, color="red",
                   label="rejected (low SNR)")
        ax.legend()
    cb = fig.colorbar(sc, ax=ax, label="pitch (2D FFT) [µm]")
    ax.set_xlabel("stage x [mm]")
    ax.set_ylabel("stage y [mm]")
    ax.set_title(f"WLI pitch across the grating ({len(spa_acc)}/{len(spa)} positions accepted)")
    ax.grid(ls=":", alpha=0.5)
    ax = axes[1]
    ax.hist(spa_acc.pitch_fft, bins=12, color="tab:blue", alpha=0.7,
            edgecolor="k")
    ax.axvline(spa_acc.pitch_fft.mean(), color="tab:red", lw=1.4,
               label=f"mean = {spa_acc.pitch_fft.mean():.5f} µm")
    ax.axvline(1.667, color="grey", ls="--", lw=1.2, label="nominal 1.667 µm")
    ax.set_xlabel("pitch (2D FFT) [µm]")
    ax.set_ylabel("count")
    ax.set_title("Distribution across positions (accepted)")
    ax.legend()
    fig.tight_layout()
    f2 = os.path.join(BASE, "WLI_spatial_map_36x.pdf")
    fig.savefig(f2, dpi=200)
    fig.savefig(f2.replace(".pdf", ".png"), dpi=150)
    print(f"Saved: {f2}")


if __name__ == "__main__":
    main()
