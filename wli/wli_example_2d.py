#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2D example figure from the Optiv WLI campaign (scan 30x/1.txt):
 (a) levelled height map,
 (b) band-passed image (0.55-0.66 cycles/um) showing the grating grooves,
 (c) centred 2D FFT log-magnitude with the grating band highlighted and the
     detected peak marked.
"""
import os
import numpy as np
import importlib.util
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("wb", os.path.join(BASE, "wli_batch.py"))
wb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(wb)

SCAN = os.path.join(BASE, "30x", "1.txt")
OUT = os.path.join(BASE, "WLI_example_2d.pdf")

Z, dx, dy, n_nan, x0, y0 = wb.load_scan(SCAN)
ny, nx = Z.shape
xx = np.arange(nx) * dx
yy = np.arange(ny) * dy
Xg, Yg = np.meshgrid(xx, yy)
A = np.column_stack([Xg.ravel(), Yg.ravel(), np.ones(Xg.size)])
coef, *_ = np.linalg.lstsq(A, Z.ravel(), rcond=None)
Zl = Z - (A @ coef).reshape(ny, nx)

# band-pass for groove visualisation
Fs = np.fft.fft2(Zl - Zl.mean())
fyv = np.fft.fftfreq(ny, d=dy)
fxv = np.fft.fftfreq(nx, d=dx)
fg = np.sqrt(fyv[:, None] ** 2 + fxv[None, :] ** 2)
band = (fg >= 0.55) & (fg <= 0.66)
Fb = np.where(band, Fs, 0)
Zg = np.fft.ifft2(Fb).real

# centred log spectrum
P = np.abs(Fs)
Ps = np.fft.fftshift(P)
kx_s = np.fft.fftshift(fxv)
ky_s = np.fft.fftshift(fyv)
fg_s = np.sqrt(ky_s[:, None] ** 2 + kx_s[None, :] ** 2)
band_s = (fg_s >= 0.55) & (fg_s <= 0.66)
# peak inside band (unshifted coords -> shifted for display)
kp = np.unravel_index(np.argmax(P * band), P.shape)
pitch = 1.0 / fg[kp]
fy_p, fx_p = fyv[kp[0]], fxv[kp[1]]
# peak position in the shifted frame
kx_shift = fxv[kp[1]]
ky_shift = fyv[kp[0]]
# map to shifted-grid indices for plotting
ix_s = int(np.argmin(np.abs(kx_s - kx_shift)))
iy_s = int(np.argmin(np.abs(ky_s - ky_shift)))

crop = (slice(400, 1000), slice(600, 1400))
fig, axes = plt.subplots(1, 3, figsize=(13, 4.3))
ax = axes[0]
im = ax.imshow(Zl[crop] * 1000, cmap="viridis", aspect="equal")
cb = fig.colorbar(im, ax=ax, fraction=0.046, label="height [µm]")
ax.set_title("(a) Levelled height map")
ax.set_xlabel("pixel"); ax.set_ylabel("pixel")

ax = axes[1]
im = ax.imshow(Zg[crop] * 1000, cmap="gray", aspect="equal")
cb = fig.colorbar(im, ax=ax, fraction=0.046, label="height [nm]")
ax.set_title("(b) Grating band (0.55–0.66 cyc/µm)")
ax.set_xlabel("pixel"); ax.set_ylabel("pixel")

ax = axes[2]
im = ax.imshow(np.log10(Ps + 1e-12), cmap="inferno",
               extent=[kx_s[0], kx_s[-1], ky_s[0], ky_s[-1]], aspect="auto")
ax.plot(kx_shift, ky_shift, "o", mfc="none", mec="white", ms=12, mew=1.5,
        label=f"peak: {pitch:.4f} µm")
ax.set_title("(c) 2D FFT (log magnitude)")
ax.set_xlabel("$f_x$ [cyc/µm]"); ax.set_ylabel("$f_y$ [cyc/µm]")
ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.1, 1.1)
ax.legend(loc="lower right", fontsize=8)

fig.suptitle("WLI example (Optiv DualZ 763, scan 30x/1)", fontsize=11)
fig.tight_layout()
fig.savefig(OUT, dpi=200)
fig.savefig(OUT.replace(".pdf", ".png"), dpi=150)
print(f"Saved: {OUT}")
