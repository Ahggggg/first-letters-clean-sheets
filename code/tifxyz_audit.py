#!/usr/bin/env python3
"""tifxyz_audit.py - mesh-level checks of a tifxyz patch (the "clean sheet" rules used to cut the patches in this repository).  MIT licence.

  python tifxyz_audit.py PATCH.tifxyz [PATCH2.tifxyz ...] [--shape Z Y X] [--json]

A tifxyz patch is a grid of (x, y, z) level-0 voxel coordinates (x.tif, y.tif, z.tif, float32, rows x cols) plus meta.json.  Rules (all from the project that produced the patches):

  frame side   per node: +1 'outward' if the grid normal (normal = d/dcol x d/drow, central differences) points away from the scan-volume XY centre, -1 'inward' if toward it; a node whose normal
               is nearly tangent to the radial direction (|cos| < 0.1) is undetermined.  frame_minority = share of the minority sign among determined nodes.  (The scroll axis is assumed to
               pass through the XY centre of the scan volume: a crude rule, stated as such; it needs the volume shape --shape Z Y X, or the meta.json keys volume_shape / the lookup in index.csv.)
  jump         per node: the largest angle between its normal and its finite 4-neighbours' normals; a 'fold' node has a jump > 30 degrees.
  stretch      an edge longer than 1.5 x the median edge length is stretched.
  clean node   finite, jump <= 30 degrees, no incident stretched edge, and on the patch's majority frame side.
  area         clean-rectangle area = nodes x (median edge length in voxels x voxel size in um x 1e-4)^2 cm2   (the area quoted everywhere in this repository);
               also the triangulated surface area of the grid, for comparison.

Every patch in data/ is a rectangle of a larger trace that was cut so that ALL its nodes are clean by these rules (largest all-clean axis-aligned node rectangle, with ground that earlier
patches already cover removed).  Re-running this script on the cropped files gives the 'rect_*' columns of index.csv (border nodes see one-sided neighbours after the crop, so a few border
nodes can differ from the values computed on the uncropped trace).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import tifffile

NORMAL_JUMP_DEG = 30.0
STRETCH = 1.5
COS_MARGIN = 0.1


def load(directory):
    base = Path(directory)
    meta = json.loads((base / "meta.json").read_text())
    xyz = np.stack([tifffile.imread(base / f"{a}.tif").astype(np.float64) for a in "xyz"], axis=-1)
    valid = np.all(np.isfinite(xyz), axis=-1) & np.all(xyz > 0, axis=-1)
    P = np.where(valid[..., None], xyz, np.nan)
    return P, grid_normals(P), meta


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


def angle_jump(N):
    """per-node max angle (degrees) to its finite 4-neighbours; NaN where no neighbour pair is finite."""
    h, w, _ = N.shape
    best = np.full((h, w), np.nan)
    for da, db in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        sh = np.full(N.shape, np.nan)
        sa = slice(0, h - 1) if da == 1 else slice(1, h) if da == -1 else slice(0, h)
        ta = slice(1, h) if da == 1 else slice(0, h - 1) if da == -1 else slice(0, h)
        sb = slice(0, w - 1) if db == 1 else slice(1, w) if db == -1 else slice(0, w)
        tb = slice(1, w) if db == 1 else slice(0, w - 1) if db == -1 else slice(0, w)
        sh[sa, sb] = N[ta, tb]
        with np.errstate(invalid="ignore"):
            ang = np.degrees(np.arccos(np.clip((N * sh).sum(-1), -1, 1)))
        best = np.where(np.isfinite(ang) & (~np.isfinite(best) | (ang > best)), ang, best)
    return best


def edge_lengths(P):
    e0 = np.linalg.norm(P[1:, :] - P[:-1, :], axis=-1)
    e1 = np.linalg.norm(P[:, 1:] - P[:, :-1], axis=-1)
    return e0, e1


def incident_max_edge(P):
    h, w, _ = P.shape
    e0, e1 = edge_lengths(P)
    inc = np.zeros((h, w))
    a = np.zeros((h, w)); b = np.zeros((h, w))
    a[:-1] = np.nan_to_num(e0, nan=0.0); b[1:] = np.nan_to_num(e0, nan=0.0)
    inc = np.maximum(inc, np.maximum(a, b))
    a = np.zeros((h, w)); b = np.zeros((h, w))
    a[:, :-1] = np.nan_to_num(e1, nan=0.0); b[:, 1:] = np.nan_to_num(e1, nan=0.0)
    return np.maximum(inc, np.maximum(a, b))


def frame_flags(P, N, shape_zyx):
    """+1 outward / -1 inward / NaN undetermined per node (shape_zyx = (slices, height, width) of the scan volume)."""
    ok = np.isfinite(P).all(-1) & np.isfinite(N).all(-1)
    c = np.array([(shape_zyx[2] - 1) / 2, (shape_zyx[1] - 1) / 2])
    dot = (N[..., :2] * (P[..., :2] - c)).sum(-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        cosr = dot / (np.linalg.norm(N[..., :2], axis=-1) * np.linalg.norm(P[..., :2] - c, axis=-1))
    s = np.where(dot > 0, 1.0, -1.0)
    s[~ok] = np.nan
    s[ok & ~(np.abs(cosr) >= COS_MARGIN)] = np.nan
    return s


def triangulated_area_vox2(P):
    a, b, c, d = P[:-1, :-1], P[:-1, 1:], P[1:, 1:], P[1:, :-1]
    t1 = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=-1)
    t2 = 0.5 * np.linalg.norm(np.cross(c - a, d - a), axis=-1)
    return float(np.nansum(t1) + np.nansum(t2))


def audit(mesh_dir, shape_zyx=None, voxel_um=None):
    P, N, meta = load(mesh_dir)
    voxel_um = voxel_um or meta.get("voxel_size_um") or meta.get("vc_gsfs_params", {}).get("voxelsize") or 9.362
    shape_zyx = shape_zyx or meta.get("volume_shape_zyx")
    fin = np.isfinite(P).all(-1) & np.isfinite(N).all(-1)
    ang = angle_jump(N)
    e0, e1 = edge_lengths(P)
    e = np.concatenate([e0.ravel(), e1.ravel()])
    e = e[np.isfinite(e)]
    med = float(np.median(e))
    inc = incident_max_edge(P)
    stretch_ok = inc <= STRETCH * med
    jump_ok = np.isfinite(ang) & (ang <= NORMAL_JUMP_DEG)
    out = {"grid": list(P.shape[:2]), "finite_nodes": int(fin.sum()), "median_edge_vox": med,
           "fold_share": float((ang[fin & np.isfinite(ang)] > NORMAL_JUMP_DEG).mean()) if (fin & np.isfinite(ang)).any() else None,
           "max_jump_deg": float(np.nanmax(ang)) if np.isfinite(ang).any() else None,
           "stretch_share": float((~stretch_ok[fin]).mean()) if fin.any() else None}
    # nodes in the outer two rings of a cropped patch see one-sided neighbours (their normals differ from the uncropped trace), so the interior is reported separately
    inner = np.zeros(fin.shape, bool)
    if min(fin.shape) > 4:
        inner[2:-2, 2:-2] = True
    ok_i = fin & np.isfinite(ang) & inner
    out["n_fold_nodes"] = int((fin & np.isfinite(ang) & (ang > NORMAL_JUMP_DEG)).sum())
    out["fold_nodes_in_outer_2_rings"] = int((fin & np.isfinite(ang) & (ang > NORMAL_JUMP_DEG) & ~inner).sum())
    out["fold_share_interior"] = float((ang[ok_i] > NORMAL_JUMP_DEG).mean()) if ok_i.any() else None
    out["max_jump_interior_deg"] = float(np.nanmax(ang[ok_i])) if ok_i.any() else None
    cell_cm2 = (med * voxel_um * 1e-4) ** 2
    out["cell_cm2"] = cell_cm2
    out["area_cm2_nodes_x_cell"] = float(fin.sum() * cell_cm2)
    out["area_cm2_triangulated"] = triangulated_area_vox2(P) * (voxel_um * 1e-4) ** 2
    if shape_zyx is not None:
        S = frame_flags(P, N, shape_zyx)
        o, i = int((S > 0).sum()), int((S < 0).sum())
        out["frame"] = {"outward": o, "inward": i, "undetermined": int(fin.sum()) - o - i, "minority": min(o, i) / max(1, o + i), "side": 1 if o >= i else -1}
        side = out["frame"]["side"]
        base = fin & jump_ok & stretch_ok & (S == side)
        out["clean_nodes_on_majority_side"] = int(base.sum())
        out["all_nodes_clean"] = bool(base.sum() == fin.sum())
    return out


def max_rectangle(mask):
    """largest all-True axis-aligned rectangle of a 2-D bool array: (area, r0, r1, c0, c1) half-open (histogram method)."""
    h, w = mask.shape
    heights = np.zeros(w, int)
    best = (0, 0, 0, 0, 0)
    for r in range(h):
        heights = np.where(mask[r], heights + 1, 0)
        stack = []
        for c in range(w + 1):
            cur = heights[c] if c < w else 0
            start = c
            while stack and stack[-1][1] >= cur:
                s, hh = stack.pop()
                area = hh * (c - s)
                if area > best[0]:
                    best = (int(area), int(r - hh + 1), int(r + 1), int(s), int(c))
                start = s
            stack.append((start, cur))
    return best


def clean_mask(P, N, shape_zyx, side):
    fin = np.isfinite(P).all(-1) & np.isfinite(N).all(-1)
    ang = angle_jump(N)
    e0, e1 = edge_lengths(P)
    e = np.concatenate([e0.ravel(), e1.ravel()])
    med = float(np.median(e[np.isfinite(e)]))
    ok = fin & np.isfinite(ang) & (ang <= NORMAL_JUMP_DEG) & (incident_max_edge(P) <= STRETCH * med)
    return ok & (frame_flags(P, N, shape_zyx) == side)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("meshes", nargs="+")
    ap.add_argument("--shape", nargs=3, type=int, metavar=("Z", "Y", "X"))
    ap.add_argument("--voxel-um", type=float)
    a = ap.parse_args()
    for m in a.meshes:
        r = audit(m, tuple(a.shape) if a.shape else None, a.voxel_um)
        print(json.dumps({"mesh": str(m), **r}, indent=1))
