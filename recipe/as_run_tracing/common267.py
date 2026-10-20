"""T-267 shared code: comparison of a tracer-grown tifxyz surface with a reference mesh, and the frame / fold audit with the largest clean rectangle.
CPU only; mesh files only (no CT read).  Reuses UNCHANGED: work/t119-pherc1203/render.py load_mesh (normal = d_col x d_row, central differences),
work/t254-framed/framed254.node_flag (outward iff n_xy . (P_xy - c) > 0, undetermined if |cos| < 0.1), work/t261-single-sheet/audit261 angle_jump / edge_stretch / flag_stats.
Conventions: tifxyz arrays (row, col, [x, y, z]); node spacing 20 voxels at scale 0.05; distances in voxels of the scan.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
PROB = HERE.parents[1]
sys.path.insert(0, str(PROB / "work/t261-single-sheet"))
sys.path.insert(0, str(PROB / "work/t254-framed"))
import common261 as C  # noqa: E402
import audit261  # noqa: E402
import framed254  # noqa: E402

DATA = PROB / "data/t267-trace"
FINE = 0.25                    # reference lattice step in node units
NEAR_VOX = 8.0                 # "within the sheet" distance (about a quarter of a layer spacing)
NORMAL_JUMP_DEG = 30.0         # as S1 (X-205)
STRETCH = 1.5                  # edge longer than 1.5 x median = stretch
MIN_AREA_NODES_CM2 = None


def load(mesh_dir):
    r119 = C.load_r119()
    gm = r119.load_mesh(str(mesh_dir))
    return gm


def fine_lattice(gm, step=FINE):
    """bilinear points and normals of the mesh on a lattice `step` node units apart, only inside quads whose four corners are finite.
    Returns (points (N,3), normals (N,3), interior (N,) bool = quad has all 8 surrounding quads' corners finite, rc (N,2) node coordinates)."""
    P, Nn = gm["points"], gm["normals"]
    h, w, _ = P.shape
    ok = np.isfinite(P).all(-1) & np.isfinite(Nn).all(-1)
    rr = np.arange(0, h - 1 + 1e-9, step)
    cc = np.arange(0, w - 1 + 1e-9, step)
    i0 = np.minimum(np.floor(rr).astype(int), h - 2)
    j0 = np.minimum(np.floor(cc).astype(int), w - 2)
    fr = (rr - i0)[:, None, None]
    fc = (cc - j0)[None, :, None]
    I, J = i0[:, None], j0[None, :]
    corners_ok = ok[I, J] & ok[I + 1, J] & ok[I, J + 1] & ok[I + 1, J + 1]

    def bil(A):
        a = A[I, J] * (1 - fr) * (1 - fc) + A[I + 1, J] * fr * (1 - fc) + A[I, J + 1] * (1 - fr) * fc + A[I + 1, J + 1] * fr * fc
        return a
    Pz = np.where(ok[..., None], P, 0.0)
    Nz = np.where(ok[..., None], Nn, 0.0)
    pts, nrm = bil(Pz), bil(Nz)
    ln = np.linalg.norm(nrm, axis=-1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        nrm = nrm / ln
    # interior: eroded valid-node mask by 2 nodes (boundary effects of the reference are not scored)
    from scipy.ndimage import binary_erosion
    core = binary_erosion(ok, iterations=2, border_value=0)
    inter = core[I, J] & core[I + 1, J + 1]
    good = corners_ok & np.isfinite(nrm).all(-1) & (ln[..., 0] > 0)
    rcgrid = np.stack(np.broadcast_arrays(rr[:, None], cc[None, :]), -1)
    return pts[good], nrm[good], inter[good], rcgrid[good]


def compare(grown, ref, seed=267, gen=None, gen_max=None, ref_nodebox=None):
    """Distances from each valid node of `grown` (dict from load_mesh) to the dense reference surface, with displaced-reference nulls.
    gen: optional (h,w) generation array of the grown mesh; nodes with gen > gen_max are dropped."""
    rp, rn, rint, rrc = fine_lattice(ref)
    tree = cKDTree(rp)
    P, Nn = grown["points"], grown["normals"]
    ok = np.isfinite(P).all(-1)
    if gen is not None and gen_max is not None:
        ok &= (gen <= gen_max) & (gen > 0)
    gp = P[ok]
    gn = Nn[ok]
    d, idx = tree.query(gp)
    comparable = rint[idx]
    out = {"grown_nodes": int(ok.sum()), "comparable_nodes": int(comparable.sum()), "edge_nodes_excluded": int((~comparable).sum())}

    def stats(dd, mask):
        x = dd[mask]
        if x.size == 0:
            return None
        return {"n": int(x.size), "median": float(np.median(x)), "p90": float(np.percentile(x, 90)), "p99": float(np.percentile(x, 99)),
                "within_3": float((x <= 3).mean()), "within_8": float((x <= NEAR_VOX).mean()), "within_20": float((x <= 20).mean())}
    out["distance"] = stats(d, comparable)
    if ref_nodebox is not None:                      # optional extra (T-267 step 2): nodes whose nearest reference point lies in a reference node box (r0, r1, c0, c1)
        r0, r1, c0, c1 = ref_nodebox
        rc = rrc[idx]
        inbox = comparable & (rc[:, 0] >= r0) & (rc[:, 0] <= r1 - 1) & (rc[:, 1] >= c0) & (rc[:, 1] <= c1 - 1)
        out["distance_in_refbox"] = stats(d, inbox)
        out["distance_outside_refbox"] = stats(d, comparable & ~inbox)
    # orientation: cosine of grown render normal and the reference normal at the nearest reference point (only nodes within NEAR_VOX)
    near = comparable & (d <= NEAR_VOX)
    gn_ok = np.isfinite(gn).all(-1)
    cos = (gn * rn[idx]).sum(-1)
    m = near & gn_ok
    out["orientation"] = {"n": int(m.sum()), "share_cos_pos": float((cos[m] > 0).mean()) if m.any() else None,
                          "share_abs_cos_ge_0p9": float((np.abs(cos[m]) >= 0.9).mean()) if m.any() else None,
                          "median_abs_cos": float(np.median(np.abs(cos[m]))) if m.any() else None}
    # nulls: reference displaced along its own normal by +-40, +-60 voxels; and by 60 voxels in 12 random directions
    rng = np.random.default_rng(seed)
    nulls = {}
    for dz in (-60, -40, 40, 60):
        t2 = cKDTree(rp + dz * rn)
        dd, ii = t2.query(gp)
        nulls[f"normal_shift_{dz:+d}"] = stats(dd, rint[ii])
    rnd = []
    for _ in range(12):
        v = rng.normal(size=3)
        v = 60 * v / np.linalg.norm(v)
        t2 = cKDTree(rp + v)
        dd, ii = t2.query(gp)
        s = stats(dd, rint[ii])
        if s:
            rnd.append(s["within_8"])
    nulls["random_shift_60_within_8"] = {"values": rnd, "max": max(rnd) if rnd else None, "mean": float(np.mean(rnd)) if rnd else None}
    out["null"] = nulls
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


def audit(mesh_dir, shape, label, cell_area_cm2=None, box_xyz=None, box_margin=0.0):
    """Frame side, neighbouring-normal jump, edge stretch and the largest clean rectangle of a whole mesh.
    clean node (for side s in {+1 outward, -1 inward}) = finite point and normal, frame flag == s, max neighbouring-normal angle <= 30 deg (all finite neighbours),
    all incident edges <= 1.5 x median edge.  Largest clean rectangle = larger of the two sides' maximal all-clean axis-aligned node rectangles."""
    F248 = framed254.load_frame248()
    S, gm, nrec = framed254.node_flag(F248, mesh_dir, shape)
    P, Nn = gm["points"], gm["normals"]
    fin = np.isfinite(P).all(-1) & np.isfinite(Nn).all(-1)
    ang = audit261.angle_jump(Nn)
    med_edge = audit261.edge_stretch(P)[1]
    h, w, _ = P.shape
    # per-node max incident edge length
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
    stretch_ok = inc <= STRETCH * med_edge
    jump_ok = np.isfinite(ang) & (ang <= NORMAL_JUMP_DEG)
    base = fin & jump_ok & stretch_ok
    if box_xyz is not None:                     # optional (T-267 wall check): only nodes inside the local-mirror voxel box shrunk by box_margin count as clean
        lo = np.array([box_xyz[0], box_xyz[2], box_xyz[4]], float) + box_margin
        hi = np.array([box_xyz[1], box_xyz[3], box_xyz[5]], float) - box_margin
        inside = np.isfinite(P).all(-1) & (P >= lo).all(-1) & (P <= hi).all(-1)
        base = base & inside
    # area per node: median edge (vox) squared x (voxel um)^2 -> cm2; voxel size passed by the caller through cell_area_cm2 (cm2 per node cell)
    res = {"label": label, "grid": [h, w], "finite_nodes": int(fin.sum()), "frame": audit261.flag_stats(S), "native_nrec": nrec,
           "jump_nodes": int((fin & np.isfinite(ang) & (ang > NORMAL_JUMP_DEG)).sum()),
           "jump_share": float((ang[fin & np.isfinite(ang)] > NORMAL_JUMP_DEG).mean()) if (fin & np.isfinite(ang)).any() else None,
           "median_edge_vox": med_edge, "stretch_share_nodes": float((~stretch_ok[fin]).mean()) if fin.any() else None,
           "median_node_angle_deg": float(np.nanmedian(ang)) if np.isfinite(ang).any() else None}
    best = None
    for side in (1.0, -1.0):
        m = base & (S == side)
        a, r0, r1, c0, c1 = max_rectangle(m)
        res[f"clean_side_{int(side):+d}"] = {"nodes": a, "rows": [r0, r1], "cols": [c0, c1], "clean_node_count": int(m.sum())}
        if best is None or a > best[0]:
            best = (a, side, r0, r1, c0, c1)
    a, side, r0, r1, c0, c1 = best
    res["largest_clean_rectangle"] = {"nodes": int(a), "side": int(side), "rows": [int(r0), int(r1)], "cols": [int(c0), int(c1)]}
    if cell_area_cm2 is not None:
        res["largest_clean_rectangle"]["area_cm2"] = a * cell_area_cm2
        res["cell_area_cm2"] = cell_area_cm2
    return res, S, ang
