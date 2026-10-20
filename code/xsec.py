"""xsec276.py - straightened CT cross-sections along a tifxyz surface, and a pooled structure-tensor measure of how the layer structure runs relative to the surface.
Standalone (numpy, scipy, tifffile, Pillow).  Registered as X-218 (Herculaneum project, T-276); MIT licence.

Conventions (the project renderer's, re-implemented): tifxyz arrays are (row, col) grids of (x, y, z) level-0 voxel coordinates; the grid is stored at `scale` (0.05 -> nodes 20 voxels apart) and a
'canvas pixel' is 1/scale grid steps (about 1 voxel); pixel (r, c) sits at grid coordinate (r*scale, c*scale) ('corner' convention); points and unit normals are bilinear in the grid; normals come
from central differences of the grid, normal = d/dcol x d/drow; CT sampling is trilinear, zero outside the volume.  CT = level-0 uint8 chunks of 128^3 voxels (zarr array order z, y, x), one raw
file per chunk at <root>/<cz>/<cy>/<cx> (a '.absent' sibling marks a chunk the store does not hold = fill value 0).  A chunk that is neither present nor marked absent is UNKNOWN: samples that
need it are flagged invalid (never read as zero).
"""
import itertools
import json
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi

CHUNK = 128


# ---------------------------------------------------------------------------------------------- mesh geometry
def load_tifxyz(directory):
    base = Path(directory)
    meta = json.loads((base / "meta.json").read_text())
    scale = meta["scale"]
    assert abs(scale[0] - scale[1]) < 1e-9, scale
    xyz = np.stack([tifffile.imread(base / f"{a}.tif").astype(np.float64) for a in "xyz"], axis=-1)
    valid = np.all(np.isfinite(xyz), axis=-1) & np.all(xyz > 0, axis=-1)
    points = np.where(valid[..., None], xyz, np.nan)
    return {"points": points, "normals": grid_normals(points), "scale": float(scale[0]), "shape": list(valid.shape), "meta": meta}


def grid_normals(points):
    def derivative(axis):
        n = points.shape[axis]

        def cut(a, b):
            return tuple(slice(a, b) if i == axis else slice(None) for i in range(3))
        central, forward, backward = (np.full_like(points, np.nan) for _ in range(3))
        central[cut(1, n - 1)] = points[cut(2, n)] - points[cut(0, n - 2)]
        forward[cut(0, n - 1)] = points[cut(1, n)] - points[cut(0, n - 1)]
        backward[cut(1, n)] = points[cut(1, n)] - points[cut(0, n - 1)]
        d = np.where(np.isnan(central), forward, central)
        return np.where(np.isnan(d), backward, d)
    normal = np.cross(derivative(1), derivative(0))
    length = np.linalg.norm(normal, axis=-1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        normal = normal / length
    normal[~(length[..., 0] > 0)] = np.nan
    return normal


def canvas_size(geom):
    h, w = geom["shape"]
    sc = geom["scale"]
    return int(np.floor((h - 1) / sc)) + 1, int(np.floor((w - 1) / sc)) + 1


def surface(geom, rows, cols):
    """points and unit normals (x, y, z) at canvas pixels rows x cols ('corner' convention); NaN where any of the four grid corners is invalid."""
    scale = geom["scale"]
    gr, gc = np.asarray(rows, np.float64) * scale, np.asarray(cols, np.float64) * scale
    height, width = geom["points"].shape[:2]
    i0, j0 = np.floor(gr).astype(int), np.floor(gc).astype(int)
    fr, fc = (gr - i0)[:, None, None], (gc - j0)[None, :, None]
    ok_r, ok_c = (i0 >= 0) & (i0 + 1 < height), (j0 >= 0) & (j0 + 1 < width)
    i0, j0 = np.clip(i0, 0, height - 2), np.clip(j0, 0, width - 2)

    def bilinear(values):
        top, bottom = values[i0], values[i0 + 1]
        out = (1 - fr) * ((1 - fc) * top[:, j0] + fc * top[:, j0 + 1]) + fr * ((1 - fc) * bottom[:, j0] + fc * bottom[:, j0 + 1])
        out[~ok_r] = np.nan
        out[:, ~ok_c] = np.nan
        return out
    points, normals = bilinear(geom["points"]), bilinear(geom["normals"])
    with np.errstate(invalid="ignore", divide="ignore"):
        normals = normals / np.linalg.norm(normals, axis=-1, keepdims=True)
    return points, normals


def frame_side(geom, shape_zyx):
    """+1 if the normal of the majority of nodes points away from the scan-volume XY centre (outward), -1 if toward it (the project's frame rule); also the minority share."""
    P, N = geom["points"], geom["normals"]
    ok = np.isfinite(P).all(-1) & np.isfinite(N).all(-1)
    c = np.array([(shape_zyx[2] - 1) / 2, (shape_zyx[1] - 1) / 2])
    dot = (N[..., :2] * (P[..., :2] - c)).sum(-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        cosr = dot / (np.linalg.norm(N[..., :2], axis=-1) * np.linalg.norm(P[..., :2] - c, axis=-1))
    det = ok & (np.abs(cosr) >= 0.1)
    out = (dot > 0) & det
    inn = (dot < 0) & det
    n_out, n_in = int(out.sum()), int(inn.sum())
    side = 1 if n_out >= n_in else -1
    minority = min(n_out, n_in) / max(1, n_out + n_in)
    return side, minority


# ---------------------------------------------------------------------------------------------- CT volume (raw chunk files)
class ChunkVolume:
    def __init__(self, roots, shape_zyx, max_chunks=96):
        self.roots = [Path(r) for r in (roots if isinstance(roots, (list, tuple)) else [roots])]
        self.shape = np.array(shape_zyx, np.int64)
        self.cache = {}
        self.max_chunks = max_chunks
        self.status = {}

    def _find(self, key):
        rel = "/".join(str(k) for k in key)
        for r in self.roots:
            p = r / rel
            if p.exists():
                return p
        for r in self.roots:
            if (r / (rel + ".absent")).exists():
                return "absent"
        return None

    def chunk(self, key):
        """array (128,128,128) uint8, 'absent' (zeros) or None (unknown)"""
        if key in self.cache:
            return self.cache[key]
        f = self._find(key)
        if f is None:
            val = None
        elif f == "absent":
            val = "absent"
        else:
            val = np.fromfile(f, np.uint8)
            if val.size != CHUNK ** 3:
                raise OSError(f"{f}: {val.size} bytes")
            val = val.reshape(CHUNK, CHUNK, CHUNK)
        if len(self.cache) >= self.max_chunks:
            self.cache.pop(next(iter(self.cache)))
        self.cache[key] = val
        return val

    def gather(self, idx):
        """voxel values at integer (z, y, x) rows of idx; outside the volume 0; returns values float64 and an `unknown` mask."""
        n = len(idx)
        vals = np.zeros(n, np.float64)
        unknown = np.zeros(n, bool)
        inside = ((idx >= 0) & (idx < self.shape)).all(1)
        ii = np.where(inside)[0]
        if len(ii) == 0:
            return vals, unknown
        ck = idx[ii] // CHUNK
        code = (ck[:, 0] * 4096 + ck[:, 1]) * 4096 + ck[:, 2]
        order = np.argsort(code, kind="stable")
        code_s = code[order]
        bounds = np.flatnonzero(np.r_[True, code_s[1:] != code_s[:-1], True])
        for a, b in zip(bounds[:-1], bounds[1:]):
            sel = ii[order[a:b]]
            key = tuple(int(v) for v in (idx[sel[0]] // CHUNK))
            arr = self.chunk(key)
            if arr is None:
                unknown[sel] = True
            elif isinstance(arr, str):
                pass
            else:
                loc = idx[sel] - np.array(key) * CHUNK
                vals[sel] = arr[loc[:, 0], loc[:, 1], loc[:, 2]]
        return vals, unknown

    def sample(self, zyx):
        """trilinear (zero outside the volume).  returns float64 values and the mask of samples that touch an unknown chunk."""
        f = np.floor(zyx).astype(np.int64)
        t = zyx - f
        out = np.zeros(len(zyx))
        bad = np.zeros(len(zyx), bool)
        for corner in itertools.product((0, 1), repeat=3):
            cz = np.array(corner)
            w = np.prod(np.where(cz == 1, t, 1 - t), axis=1)
            v, u = self.gather(f + cz)
            out += w * v
            bad |= u & (w > 0)
        return out, bad


# ---------------------------------------------------------------------------------------------- straightened sections
def section(vol, geom, kind, index, span, depths, side=1):
    """straightened section.  kind 'row': canvas row `index`, canvas columns span=(c0, c1); kind 'col': canvas column `index`, canvas rows span.
    Samples p + side*d*n for d in depths (so +d is away from the scroll axis when side is the patch frame side).  Returns S (n_depth, n_s) float32, valid (n_s) bool, arclength in voxels, canvas coordinates."""
    s_idx = np.arange(span[0], span[1])
    if kind == "row":
        P, N = surface(geom, [index], s_idx)
        P, N = P[0], N[0]
    else:
        P, N = surface(geom, s_idx, [index])
        P, N = P[:, 0], N[:, 0]
    good = np.isfinite(P).all(-1) & np.isfinite(N).all(-1)
    d = np.asarray(depths, np.float64)
    S = np.zeros((len(d), len(s_idx)), np.float32)
    valid = good.copy()
    if good.any():
        base, direction = P[good][:, ::-1], N[good][:, ::-1]          # z, y, x
        samples = base[None] + (side * d)[:, None, None] * direction[None]
        vals, bad = vol.sample(samples.reshape(-1, 3))
        vals = vals.reshape(len(d), -1)
        bad = bad.reshape(len(d), -1)
        S[:, good] = vals
        colbad = bad.any(0)
        vi = np.where(good)[0]
        valid[vi[colbad]] = False
    step = np.linalg.norm(np.diff(np.where(good[:, None], P, np.nan), axis=0), axis=-1)
    arclen = float(np.nansum(step)) if step.size else 0.0
    return S, valid, arclen, s_idx


# ---------------------------------------------------------------------------------------------- structure-tensor measure
SIGMA = 2.0
WINDOW = 24            # |d| <= 24 voxels enters the sums


def _tensor(S, valid, depths, sigma=SIGMA, window=WINDOW):
    S = S.astype(np.float64)
    if (~valid).any() and valid.any():
        S = S.copy()
        S[:, ~valid] = S[:, valid].mean()
    gs = ndi.gaussian_filter(S, sigma=sigma, order=(0, 1), mode="nearest")
    gd = ndi.gaussian_filter(S, sigma=sigma, order=(1, 0), mode="nearest")
    rows = np.abs(np.asarray(depths)) <= window
    # drop columns within 4 px of an invalid column
    ok = valid.copy()
    if (~valid).any():
        ok = ~ndi.binary_dilation(~valid, iterations=4)
    ok[:4] = False
    ok[-4:] = False
    a, b = gs[rows][:, ok], gd[rows][:, ok]
    if a.size == 0:
        return None
    return float((a * a).sum()), float((a * b).sum()), float((b * b).sum()), int(ok.sum())


def tilt_coh(J):
    if J is None:
        return float("nan"), float("nan")
    Jss, Jsd, Jdd, _ = J
    tot = Jss + Jdd
    if tot <= 0:
        return float("nan"), float("nan")
    tg = 0.5 * np.degrees(np.arctan2(2 * Jsd, Jss - Jdd))
    tilt = (tg % 180.0) - 90.0
    coh = float(np.sqrt((Jss - Jdd) ** 2 + 4 * Jsd ** 2) / tot)
    return float(tilt), coh


def measure(S, valid, depths, rng, n_null=20):
    """tilt (degrees; 0 = layers parallel to the surface), coherence, noise threshold (99th pct of coherence over n_null fully shuffled copies), informative flag."""
    J = _tensor(S, valid, depths)
    tilt, coh = tilt_coh(J)
    if not np.isfinite(coh):
        return {"tilt": np.nan, "coherence": np.nan, "threshold": np.nan, "informative": False, "n_valid_cols": int(valid.sum())}
    vals = S[:, valid].ravel()
    cn = []
    for _ in range(n_null):
        Sh = S.copy()
        Sh[:, valid] = rng.permutation(vals).reshape(S[:, valid].shape)
        cn.append(tilt_coh(_tensor(Sh, valid, depths))[1])
    thr = float(np.percentile(cn, 99))
    return {"tilt": tilt, "coherence": coh, "threshold": thr, "informative": bool(coh > thr), "n_valid_cols": int(valid.sum())}


def local_tilts(S, valid, depths, win=384, stride=192):
    out = []
    for a in range(0, max(1, S.shape[1] - win + 1), stride):
        t, c = tilt_coh(_tensor(S[:, a:a + win], valid[a:a + win], depths))
        out.append((a + win // 2, t, c))
    return out


def shear(S, depths, a_deg):
    """S'(d, s) = S(d - (s - s_mid) tan a, s), linear interpolation in d, edge rows clamped."""
    d = np.asarray(depths, np.float64)
    ns = S.shape[1]
    sm = (ns - 1) / 2.0
    ta = np.tan(np.radians(a_deg))
    out = np.empty_like(S)
    for j in range(ns):
        src = d - (j - sm) * ta
        out[:, j] = np.interp(src, d, S[:, j].astype(np.float64))
    return out
