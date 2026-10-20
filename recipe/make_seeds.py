#!/usr/bin/env python3
"""make_seeds.py - tracer seed rule used for the patches of this repository (standalone copy of the project's seeds268.py; MIT licence).  CPU, numpy + scipy + numcodecs.

INPUT: a local mirror of the public m7 surface prediction, level 0: MIRROR/pred.zarr/0/<cz>/<cy>/<cx> (uint8 0/255, 192^3 voxel blosc chunks; absent chunk = 0), and the mirror's voxel box.
RULE (fixed before any tracer run): cells = 2 x 2 x 2 chunks (384 voxels) aligned to the chunk grid.  In a cell take the mask at stride 2 (192^3).  Window = 11 stride-2 voxels (22 voxels), zero outside
the cell.  Candidate = mask voxel with window mask fraction M0 in [0.15, 0.50] (a thin band that fills a fifth to a half of the window, so not a fused mass), at least 12 stride-2 voxels from the cell
faces, whose distance to the window's mask centroid is <= 1.5 stride-2 voxels (near the centre of the band, not its edge).  From at most 6000 candidates (numpy rng seed 268 + chunk-block hash) compute
the window covariance of the mask voxels and r = lambda_min / lambda_mid (small = plate).  The cell's seed = the candidate with the smallest r, kept iff r <= 0.10, the seed lies >= 700 voxels inside
the mirror box, and the seed voxel is 255.  Cells are listed by r ascending (flattest first).

  python make_seeds.py MIRROR_DIR --box XLO XHI YLO YHI ZLO ZHI --out seed-candidates.json [--procs 4]

Check: on the project's PHerc0813 mirror this gives the same seeds as the project's stored seed-candidates.json (see checks/).
"""
import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

CH = 192
WALL_MARGIN_VOX = 700.0


class C:                       # constants namespace (the project's common268)
    CH = CH
    WALL_MARGIN_VOX = WALL_MARGIN_VOX


W = 11
M0_LO, M0_HI = 0.15, 0.50
OFF_MAX = 1.5
MARGIN_DS = 12
R_MAX = 0.10
MAX_CAND = 6000


def one_cell(args):
    import numcodecs
    from scipy.ndimage import uniform_filter
    mdir, cz, cy, cx = args
    codec = numcodecs.Blosc()
    a = np.zeros((2 * C.CH,) * 3, np.uint8)
    n = 0
    for dz in (0, 1):
        for dy in (0, 1):
            for dx in (0, 1):
                f = Path(mdir) / "pred.zarr" / "0" / str(cz * 2 + dz) / str(cy * 2 + dy) / str(cx * 2 + dx)
                if f.exists():
                    a[dz * C.CH:(dz + 1) * C.CH, dy * C.CH:(dy + 1) * C.CH, dx * C.CH:(dx + 1) * C.CH] = np.frombuffer(codec.decode(f.read_bytes()), np.uint8).reshape(C.CH, C.CH, C.CH)
                    n += 1
    if n == 0:
        return None
    m = (a[::2, ::2, ::2] > 0).astype(np.float32)
    if not m.any():
        return None
    s = m.shape[0]
    ax = (np.arange(s, dtype=np.float32) - s / 2)
    crd = [ax[:, None, None], ax[None, :, None], ax[None, None, :]]
    flt = lambda v: uniform_filter(v, W, mode="constant")
    M0 = flt(m)
    sel = (m > 0) & (M0 >= M0_LO) & (M0 <= M0_HI)
    sel[:MARGIN_DS] = sel[-MARGIN_DS:] = False
    sel[:, :MARGIN_DS] = sel[:, -MARGIN_DS:] = False
    sel[..., :MARGIN_DS] = sel[..., -MARGIN_DS:] = False
    idx = np.argwhere(sel)
    if len(idx) == 0:
        return {"cell": [cz, cy, cx], "chunks": n, "candidates": 0, "seed": None}
    rng = np.random.default_rng(268 + (cz * 1009 + cy) * 1009 + cx)
    if len(idx) > MAX_CAND:
        idx = idx[rng.choice(len(idx), MAX_CAND, replace=False)]
    M0s = np.maximum(M0, 1e-6)
    mu = [flt(m * c) / M0s for c in crd]
    cand = tuple(idx.T)
    off = np.sqrt(sum((i - mu_[cand]) ** 2 for i, mu_ in zip(idx.T.astype(np.float32) - s / 2, mu)))
    keep = off <= OFF_MAX
    idx, off = idx[keep], off[keep]
    if len(idx) == 0:
        return {"cell": [cz, cy, cx], "chunks": n, "candidates": 0, "seed": None}
    cand = tuple(idx.T)
    cov = np.empty((len(idx), 3, 3), np.float32)
    for i in range(3):
        for j in range(i, 3):
            v = flt(m * crd[i] * crd[j])[cand] / M0s[cand] - mu[i][cand] * mu[j][cand]
            cov[:, i, j] = v
            cov[:, j, i] = v
    ev = np.linalg.eigvalsh(cov)
    r = ev[:, 0] / np.maximum(ev[:, 1], 1e-6)
    k = int(np.argmin(r))
    zi, yi, xi = (int(v) for v in idx[k])
    gz, gy, gx = cz * 2 * C.CH + 2 * zi, cy * 2 * C.CH + 2 * yi, cx * 2 * C.CH + 2 * xi
    return {"cell": [cz, cy, cx], "chunks": n, "candidates": int(len(idx)), "seed": [gx, gy, gz], "r": float(r[k]), "M0": float(M0[zi, yi, xi]), "mask_at_seed": int(a[2 * zi, 2 * yi, 2 * xi]),
            "r_median": float(np.median(r))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mirror")
    ap.add_argument("--box", nargs=6, type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--procs", type=int, default=4)
    a = ap.parse_args()
    mdir = Path(a.mirror)
    box = a.box
    present = set()
    for z in (mdir / "pred.zarr" / "0").iterdir():
        if not z.name.isdigit():
            continue
        for y in z.iterdir():
            for x in y.iterdir():
                present.add((int(z.name) // 2, int(y.name) // 2, int(x.name) // 2))
    cells = sorted(present)
    print(a.mirror, "cells with data", len(cells), flush=True)
    t0 = time.time()
    out = []
    with ProcessPoolExecutor(a.procs) as ex:
        for i, res in enumerate(ex.map(one_cell, [(str(mdir), *c) for c in cells], chunksize=4)):
            if res is not None:
                out.append(res)
            if (i + 1) % 100 == 0:
                print(i + 1, "cells", round(time.time() - t0), "s", flush=True)
    xlo, xhi, ylo, yhi, zlo, zhi = box
    seeds = []
    for r in out:
        s = r.get("seed")
        if not s or r["r"] > R_MAX or r["mask_at_seed"] == 0:
            continue
        x, y, z = s
        if x < xlo + C.WALL_MARGIN_VOX or x > xhi - C.WALL_MARGIN_VOX or y < ylo + C.WALL_MARGIN_VOX or y > yhi - C.WALL_MARGIN_VOX or z < zlo + C.WALL_MARGIN_VOX or z > zhi - C.WALL_MARGIN_VOX:
            continue
        seeds.append(r)
    seeds.sort(key=lambda r: (r["r"], r["cell"]))
    for k, r in enumerate(seeds):
        r["no"] = k
    res = {"mirror": str(a.mirror), "rule": __doc__.split("RULE")[1].split("INPUT")[0].strip(), "box_xyz": box, "cells_with_data": len(cells), "cells_with_seed_any": sum(1 for r in out if r.get("seed")),
           "seeds": seeds, "n_seeds": len(seeds), "seconds": round(time.time() - t0, 1)}
    Path(a.out).write_text(json.dumps(res, indent=1))
    print("seeds kept", len(seeds), "of", len(cells), "cells;", "r pct", np.percentile([r["r"] for r in seeds], [0, 25, 50, 75, 100]).round(4).tolist() if seeds else None, "s", res["seconds"])


if __name__ == "__main__":
    main()
