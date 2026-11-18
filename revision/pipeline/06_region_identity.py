#!/usr/bin/env python3
"""Do two series image the same area of the grating?  (Referee 2, comment 5)

For every tile of the first series the BED-S frame is cross-correlated with
the BED-S frame of the second series recorded at the nearest stage
coordinates ("matched") and, as a control, with the frame recorded farthest
away.  The groove pattern is removed first (column and row means), so that
only the non-periodic surface texture contributes.  A common area gives a
correlation peak well above the control level at a consistent offset.

Output: revision/analysis/referee_region_identity.csv and a printed summary.
"""
import os

import numpy as np
import pandas as pd
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(HERE)
ROOT = os.path.dirname(REV)
OUT = os.path.join(REV, "analysis")
PX = 0.01333594
PAIRS = [(1000, 750), (750, 500), (500, 400), (1000, 400)]
MAX_TILES = 24


def texture(rel):
    a = np.asarray(Image.open(os.path.join(ROOT, rel)).convert("L"), dtype=np.float32)
    a = a - a.mean(axis=0, keepdims=True)
    a = a - a.mean(axis=1, keepdims=True)
    return a / a.std()


def peak(a, b):
    c = np.fft.irfft2(np.fft.rfft2(a) * np.conj(np.fft.rfft2(b)), s=a.shape) / a.size
    k = np.unravel_index(np.argmax(c), c.shape)
    dy = k[0] if k[0] < a.shape[0] / 2 else k[0] - a.shape[0]
    dx = k[1] if k[1] < a.shape[1] / 2 else k[1] - a.shape[1]
    return dx * PX, dy * PX, float((c.max() - c.mean()) / c.std())


meta = pd.read_csv(os.path.join(OUT, "sem_metadata.csv"))
meta = meta[meta.detector == "BED-S"]
rows = []
for first, second in PAIRS:
    A, B = meta[meta.group_label == first], meta[meta.group_label == second]
    A = A.sample(min(MAX_TILES, len(A)), random_state=1)
    for r in A.itertuples():
        dist = np.hypot(B.stage_x_um - r.stage_x_um, B.stage_y_um - r.stage_y_um)
        ta = texture(r.image_path)
        dx, dy, z = peak(ta, texture(B.loc[dist.idxmin()].image_path))
        _, _, zc = peak(ta, texture(B.loc[dist.idxmax()].image_path))
        rows.append(dict(first=first, second=second, stage_distance_um=dist.min(), shift_x_um=dx,
                         shift_y_um=dy, z_matched=z, z_control=zc))
res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "referee_region_identity.csv"), index=False)
for (first, second), g in res.groupby(["first", "second"], sort=False):
    hit = g[g.z_matched > g.z_control.max()]
    print(f"{first}x vs {second}x: tiles {len(g)}, matched z median {g.z_matched.median():.1f} "
          f"(min {g.z_matched.min():.1f}), control z median {g.z_control.median():.1f} (max {g.z_control.max():.1f}); "
          f"matched above all controls: {len(hit)}; their shift x median {hit.shift_x_um.median():.2f} "
          f"(range {hit.shift_x_um.min():.2f} to {hit.shift_x_um.max():.2f}) um, y median {hit.shift_y_um.median():.2f} um")
