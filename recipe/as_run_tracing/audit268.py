"""T-268 (X-212) per-trace audit.  Same thresholds and code logic as work/t267-trace/common267.audit (X-205 rules: frame flag by framed254.node_flag, neighbouring-normal jump > 30 degrees,
edge > 1.5 x median = stretch), re-expressed so that (a) a mirror-box wall mask and (b) an exclusion of nodes covered by earlier kept patches can be applied before the largest rectangle is searched.
Largest clean rectangle = larger over the two frame sides of the maximal all-clean axis-aligned node rectangle (common267.max_rectangle, histogram method).
  python audit268.py selftest      B1 of X-210 must give 23,865 nodes, rows 2-217, cols 2-113 (the registered equivalence check)
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common268 as C  # noqa: E402
import common267 as K  # noqa: E402

_F248 = None


def f248():
    global _F248
    if _F248 is None:
        _F248 = K.framed254.load_frame248()
    return _F248


def clean_state(mesh_dir, shape, box_xyz=None):
    S, gm, nrec = K.framed254.node_flag(f248(), mesh_dir, shape)
    P, Nn = gm["points"], gm["normals"]
    fin = np.isfinite(P).all(-1) & np.isfinite(Nn).all(-1)
    ang = K.audit261.angle_jump(Nn)
    med_edge = K.audit261.edge_stretch(P)[1]
    h, w, _ = P.shape
    inc = np.zeros((h, w))
    for axis in (0, 1):
        e = np.linalg.norm(np.diff(P, axis=axis), axis=-1)
        a = np.zeros((h, w))
        b = np.zeros((h, w))
        if axis == 0:
            a[:-1] = np.nan_to_num(e, nan=0.0)
            b[1:] = np.nan_to_num(e, nan=0.0)
        else:
            a[:, :-1] = np.nan_to_num(e, nan=0.0)
            b[:, 1:] = np.nan_to_num(e, nan=0.0)
        inc = np.maximum(inc, np.maximum(a, b))
    stretch_ok = inc <= K.STRETCH * med_edge
    jump_ok = np.isfinite(ang) & (ang <= K.NORMAL_JUMP_DEG)
    base = fin & jump_ok & stretch_ok
    inside = np.ones_like(base)
    if box_xyz is not None:
        lo = np.array([box_xyz[0], box_xyz[2], box_xyz[4]], float)
        hi = np.array([box_xyz[1], box_xyz[3], box_xyz[5]], float)
        inside = np.isfinite(P).all(-1) & (P >= lo).all(-1) & (P <= hi).all(-1)
        base = base & inside
    return {"S": S, "P": P, "N": Nn, "fin": fin, "ang": ang, "base": base, "med_edge": med_edge, "gm": gm, "nrec": nrec, "stretch_ok": stretch_ok, "inside": inside, "jump_ok": jump_ok}


def best_rect(st, covered=None):
    """largest clean rectangle over both sides, with the optional boolean `covered` node mask removed.  returns dict or None"""
    best = None
    for side in (1.0, -1.0):
        m = st["base"] & (st["S"] == side)
        if covered is not None:
            m = m & ~covered
        a, r0, r1, c0, c1 = K.max_rectangle(m)
        if best is None or a > best["nodes"]:
            best = {"nodes": int(a), "side": int(side), "rows": [int(r0), int(r1)], "cols": [int(c0), int(c1)]}
    return best


def covered_mask(P, tree, vox=C.EXCL_VOX):
    """nodes whose nearest point of `tree` (cKDTree of earlier kept nodes) is within vox"""
    cov = np.zeros(P.shape[:2], bool)
    if tree is None:
        return cov
    fin = np.isfinite(P).all(-1)
    d, _ = tree.query(P[fin])
    cov[fin] = d <= vox
    return cov


def summarize(st, voxel_um):
    P, S, ang, fin = st["P"], st["S"], st["ang"], st["fin"]
    cell_cm2 = (st["med_edge"] * voxel_um * 1e-4) ** 2
    fl = K.audit261.flag_stats(S)
    return {"grid": list(P.shape[:2]), "finite_nodes": int(fin.sum()), "frame": fl, "frame_minority": fl.get("minority"),
            "jump_nodes": int((fin & np.isfinite(ang) & (ang > K.NORMAL_JUMP_DEG)).sum()),
            "fold_share": float((ang[fin & np.isfinite(ang)] > K.NORMAL_JUMP_DEG).mean()) if (fin & np.isfinite(ang)).any() else None,
            "stretch_share": float((~st["stretch_ok"][fin]).mean()) if fin.any() else None, "median_edge_vox": float(st["med_edge"]), "cell_cm2": cell_cm2,
            "outside_box_nodes": int((fin & ~st["inside"]).sum()), "clean_nodes": int(st["base"].sum()), "median_node_angle_deg": float(np.nanmedian(ang)) if np.isfinite(ang).any() else None}


def selftest():
    mesh = C.DATA267 / "runs/B1_seed75_75_G110/out/seg"
    st = clean_state(mesh, C.vol_info("PHerc0826")["shape"], None)
    r = best_rect(st)
    print(json.dumps({"rect": r, "summary": summarize(st, 9.362)}, indent=1))
    ok = r["nodes"] == 23865 and r["rows"] == [2, 217] and r["cols"] == [2, 113] and r["side"] == -1
    # equivalence with the X-210 audit json
    ref = json.loads((C.DATA267 / "audit-B1.json").read_text())
    s = summarize(st, 9.362)
    ok = ok and s["finite_nodes"] == ref["finite_nodes"] and s["jump_nodes"] == ref["jump_nodes"] and abs(s["frame"]["minority"] - ref["frame"]["minority"]) < 1e-12
    print("SELFTEST", "PASS" if ok else "FAIL")
    # exclusion check: cover the left half of the rectangle -> the new rectangle must lie entirely off the covered nodes
    from scipy.spatial import cKDTree
    P = st["P"]
    cov_pts = P[2:217, 2:60].reshape(-1, 3)
    cov = covered_mask(P, cKDTree(cov_pts))
    r2 = best_rect(st, cov)
    r0, r1 = r2["rows"]
    c0, c1 = r2["cols"]
    ok2 = not cov[r0:r1, c0:c1].any() and r2["nodes"] > 0 and r2["nodes"] < r["nodes"]
    print("EXCLUSION", r2, "PASS" if ok2 else "FAIL")
    return ok and ok2


if __name__ == "__main__":
    if sys.argv[1] == "selftest":
        sys.exit(0 if selftest() else 1)
