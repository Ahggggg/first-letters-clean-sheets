#!/usr/bin/env python3
"""cut_patches.py - cut the largest clean single-frame, fold-free rectangle out of a traced tifxyz (the step that turns a whole trace into a patch of this repository).  MIT licence.

  python cut_patches.py TRACE_SEG_DIR --shape Z Y X [--voxel-um 9.362] [--min-cm2 1.0] [--exclude PATCH.tifxyz ...] [--box XLO XHI YLO YHI ZLO ZHI] [--out DIR --name NAME]

Rules (identical to the project's audit, see code/tifxyz_audit.py):  a node is clean if it is finite, has a neighbouring-normal jump <= 30 degrees, no edge longer than 1.5 x the median edge, and lies on the
chosen frame side; nodes within 12 voxels of an --exclude patch (ground already kept) and nodes outside --box (the edge of the local data) are not clean.  The patch is the larger of the two sides' maximal
all-clean axis-aligned node rectangles; it is kept if its area (nodes x (median edge x voxel um)^2) is >= --min-cm2.  Prints the rectangle; with --out writes the cropped tifxyz.
Check: on the stored trace PHerc0813_s000 (seed 5028 4098 9036) this reproduces rows 69-117, cols 58-133 (3,600 nodes, 1.30 cm2) of data/PHerc0813/PHerc0813_t269_s000.tifxyz (see checks/).
"""
import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import tifffile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
import tifxyz_audit as A  # noqa: E402

EXCL_VOX = 12.0


def covered_mask(P, others):
    from scipy.spatial import cKDTree
    pts = []
    for m in others:
        Q, _, _ = A.load(m)
        pts.append(Q[np.isfinite(Q).all(-1)])
    cov = np.zeros(P.shape[:2], bool)
    if not pts:
        return cov
    tree = cKDTree(np.concatenate(pts))
    fin = np.isfinite(P).all(-1)
    d, _ = tree.query(P[fin])
    cov[fin] = d <= EXCL_VOX
    return cov


def cut(mesh_dir, shape_zyx, voxel_um=9.362, min_cm2=1.0, exclude=(), box=None):
    P, N, meta = A.load(mesh_dir)
    fin = np.isfinite(P).all(-1) & np.isfinite(N).all(-1)
    ang = A.angle_jump(N)
    e0, e1 = A.edge_lengths(P)
    e = np.concatenate([e0.ravel(), e1.ravel()])
    med = float(np.median(e[np.isfinite(e)]))
    base = fin & np.isfinite(ang) & (ang <= A.NORMAL_JUMP_DEG) & (A.incident_max_edge(P) <= A.STRETCH * med)
    if box is not None:
        lo = np.array([box[0], box[2], box[4]], float)
        hi = np.array([box[1], box[3], box[5]], float)
        base &= np.isfinite(P).all(-1) & (P >= lo).all(-1) & (P <= hi).all(-1)
    cov = covered_mask(P, exclude)
    S = A.frame_flags(P, N, shape_zyx)
    best = None
    for side in (1.0, -1.0):
        a, r0, r1, c0, c1 = A.max_rectangle(base & (S == side) & ~cov)
        if best is None or a > best["nodes"]:
            best = {"nodes": int(a), "side": int(side), "rows": [int(r0), int(r1)], "cols": [int(c0), int(c1)]}
    best["cm2"] = best["nodes"] * (med * voxel_um * 1e-4) ** 2
    best["median_edge_vox"] = med
    best["keep"] = bool(best["cm2"] >= min_cm2)
    best["frame_correct_order"] = "F" if best["side"] == 1 else "R"
    return best, P, meta


def write_crop(mesh_dir, rect, out_dir, name, meta):
    r0, r1 = rect["rows"]
    c0, c1 = rect["cols"]
    out = Path(out_dir) / f"{name}.tifxyz"
    out.mkdir(parents=True, exist_ok=True)
    for a in "xyz":
        arr = tifffile.imread(Path(mesh_dir) / f"{a}.tif")[r0:r1, c0:c1]
        tifffile.imwrite(out / f"{a}.tif", arr)
    m = dict(meta)
    m.update({"uuid": name, "crop_rows": [r0, r1], "crop_cols": [c0, c1], "area_cm2": rect["cm2"], "note": "crop of the largest clean rectangle (cut_patches.py)"})
    (out / "meta.json").write_text(json.dumps(m, indent=1))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mesh")
    ap.add_argument("--shape", nargs=3, type=int, required=True, metavar=("Z", "Y", "X"))
    ap.add_argument("--voxel-um", type=float, default=9.362)
    ap.add_argument("--min-cm2", type=float, default=1.0)
    ap.add_argument("--exclude", nargs="*", default=[])
    ap.add_argument("--box", nargs=6, type=float)
    ap.add_argument("--out")
    ap.add_argument("--name")
    a = ap.parse_args()
    rect, P, meta = cut(a.mesh, tuple(a.shape), a.voxel_um, a.min_cm2, a.exclude, a.box)
    print(json.dumps(rect, indent=1))
    if a.out and rect["keep"]:
        print("wrote", write_crop(a.mesh, rect, a.out, a.name or Path(a.mesh).parent.name, meta))
