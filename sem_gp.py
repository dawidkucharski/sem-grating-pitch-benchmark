#!/usr/bin/env python3
"""Numerically robust GP analysis of tile-level linearity data.
Replaces gp_linearity.R (whose chol() solve silently failed for these
ill-conditioned kernels). Fits an SE-kernel GP on normalised tile
coordinates per magnification, computes the LRT against a noise-only
null, the physical lengthscale (stage extents) and the predictive-mean
range, writes gp_linearity_summary.csv and gp_linearity_combined.pdf."""
import os
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/dawid/Library/Mobile Documents/com~apple~CloudDocs/Dokumenty/SEM_measurements/Calibration"
RDIR = os.path.join(BASE, "R_output_dir")
GB = os.path.join(BASE, "sem_grating_benchmark")
MAGS = [400, 500, 750, 1000]
STAGE = {400: 305, 500: 237, 750: 153, 1000: 119}
PITCH = 1.644

def solve_psd(K, y):
    """Stable solve K alpha = y with log-det, PSD-safe."""
    n = len(y)
    K = K + 1e-12 * np.eye(n)
    try:
        L = np.linalg.cholesky(K)
        al = np.linalg.solve(L.T, np.linalg.solve(L, y))
        ld = 2 * np.log(np.diag(L)).sum()
    except np.linalg.LinAlgError:
        w, V = np.linalg.eigh(K)
        w = np.clip(w, 1e-14, None)
        al = V @ ((V.T @ y) / w)
        ld = np.log(w).sum()
    return al, ld

def fit(mag):
    df = pd.read_csv(os.path.join(RDIR, f"linearity_tiles_mag{mag}.csv"))
    med = df.local_pitch_um.median()
    madv = (df.local_pitch_um - med).abs().median() * 1.4826
    df = df[(df.local_pitch_um - med).abs() < 3 * madv].reset_index(drop=True)
    y = df.local_pitch_um.values - df.local_pitch_um.mean()
    X = np.c_[(df.tile_x - df.tile_x.min()) / max(1, df.tile_x.max() - df.tile_x.min()),
              (df.tile_y - df.tile_y.min()) / max(1, df.tile_y.max() - df.tile_y.min())]
    n = len(y)
    sq = ((X[:, None, :] - X[None, :, :]) ** 2).sum(-1)
    var0 = np.mean(y ** 2)
    nll0 = n / 2 * (np.log(2 * np.pi * var0) + 1)

    def nll(p):
        s2, ls, sn = np.exp(p)
        K = s2 * np.exp(-0.5 * sq / ls ** 2) + sn * np.eye(n)
        al, ld = solve_psd(K, y)
        return float(0.5 * y @ al + 0.5 * ld + 0.5 * n * np.log(2 * np.pi))

    rng = np.random.default_rng(7)
    inits = [[np.log(var0), np.log(0.3), np.log(var0 * 0.1)],
             [np.log(var0), np.log(0.2), np.log(var0 * 0.05)],
             [np.log(var0 * 0.5), np.log(0.5), np.log(var0 * 0.1)],
             [np.log(var0 * 0.2), np.log(0.7), np.log(var0 * 0.02)],
             [np.log(var0 * 0.1), np.log(1.0), np.log(var0 * 0.05)]]
    for _ in range(8):
        inits.append([np.log(var0 * rng.uniform(0.05, 2)),
                      np.log(rng.uniform(0.05, 1.9)),
                      np.log(var0 * rng.uniform(0.005, 0.5))])
    best = None
    for init in inits:
        r = minimize(nll, init, method="L-BFGS-B",
                     bounds=[(np.log(1e-6), 0), (np.log(0.01), np.log(2)),
                             (np.log(1e-8), 0)], options={"maxiter": 1000})
        if best is None or r.fun < best.fun:
            best = r
    s2, ls, sn = np.exp(best.x)
    K = s2 * np.exp(-0.5 * sq / ls ** 2) + sn * np.eye(n)
    al, ld = solve_psd(K, y)
    lrt = 2 * (nll0 - best.fun)
    pval = 1 - stats.chi2.cdf(lrt, 2)
    # predictions
    g = np.linspace(0, 1, 50)
    G = np.array([(a, b) for a in g for b in g])
    sqg = ((G[:, None, :] - X[None, :, :]) ** 2).sum(-1)
    fg = s2 * np.exp(-0.5 * sqg / ls ** 2) @ al
    ftr = s2 * np.exp(-0.5 * sq / ls ** 2) @ al
    frms = float(np.sqrt(np.mean(ftr ** 2)))
    maxdist = 100 * (fg.max() - fg.min()) / PITCH
    ls_phys = ls * STAGE[mag]
    print(f"{mag}x: n={n} sigma2={s2:.3e} ls_norm={ls:.4f} ls_um={ls_phys:.1f} "
          f"noise={sn:.3e} nlml={best.fun:.3f} LRT={lrt:.2f} p={pval:.4g} "
          f"pred_range={maxdist:.2f}% fRMS={frms:.5f} um")
    return dict(magnification=mag, n_tiles=n, sigma2=s2, lengthscale_norm=ls,
                lengthscale_um=ls_phys, noise_variance=sn,
                max_distortion_pct=maxdist, nlml=best.fun,
                lrt=lrt, lrt_p=pval, f_rms_um=frms,
                grid=fg.reshape(50, 50), y=y, X=X, df=df)

rows = [fit(m) for m in MAGS]
sum_df = pd.DataFrame([{k: r[k] for k in ["magnification", "n_tiles", "sigma2",
                                          "lengthscale_norm", "lengthscale_um",
                                          "noise_variance", "max_distortion_pct",
                                          "nlml", "lrt", "lrt_p", "f_rms_um"]}
                       for r in rows])
for out in [os.path.join(RDIR, "gp_linearity_summary.csv"),
            os.path.join(GB, "gp_linearity_summary.csv")]:
    sum_df.to_csv(out, index=False)
    print("wrote", out)

# figure (2x2 smooth-field maps, per-panel colour bars to avoid overlap)
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, axes = plt.subplots(2, 2, figsize=(10.0, 8.6))
for ax, r in zip(axes.ravel(), rows):
    z = r["grid"] + r["df"].local_pitch_um.mean()
    zm = z - z.mean()
    vmax = np.abs(zm).max()
    im = ax.imshow(zm, extent=(0, 1, 0, 1), origin="lower", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax)
    ax.contour(zm, levels=np.linspace(-vmax, vmax, 9)[1:-1],
               extent=(0, 1, 0, 1), origin="lower", colors="0.25", linewidths=0.5)
    m = r["magnification"]
    ax.set_title(f'{m}\u00d7  (lengthscale = {r["lengthscale_um"]:.0f} \u00b5m)\n'
                 f'peak deviation \u00b1{vmax * 1000:.1f} nm',
                 fontsize=9.5)
    ax.set_xlabel("Normalised X", fontsize=9)
    ax.set_ylabel("Normalised Y", fontsize=9)
    ax.tick_params(labelsize=8)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.06)
    fig.colorbar(im, cax=cax)
    cax.tick_params(labelsize=7)
fig.subplots_adjust(left=0.08, right=0.95, top=0.93, bottom=0.08,
                    wspace=0.35, hspace=0.55)
fig.savefig(os.path.join(BASE, "gp_linearity_combined.pdf"))
print("wrote gp_linearity_combined.pdf")
