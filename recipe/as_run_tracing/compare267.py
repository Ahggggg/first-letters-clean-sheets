"""T-267 driver (CPU, mesh files only, fixed seed 267).
  compare --grown DIR --ref DIR --label L --out JSON [--gen-max N]   distances of grown nodes to a dense reference surface, orientation, displaced-reference nulls
  audit   --mesh DIR --scroll S --label L --out JSON                 frame side, folds, stretch, largest clean rectangle of one mesh
  cleanref --grown DIR --ref DIR --audit JSON --refbox r0 r1 c0 c1 --label L --out JSON   share of the largest clean rectangle lying on / in a reference node box
  refs    --out JSON                                                 the same audit on the 24 official PHerc0826 factory meshes (the null for "how big is a clean rectangle")
python compare267.py ..."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common267 as K  # noqa: E402

SHAPES = {"PHerc0826": (16920, 8169, 8169), "PHerc0139": (20974, 6621, 6621)}
VOX_UM = 9.362


def cell_cm2(med_edge_vox):
    return (med_edge_vox * VOX_UM * 1e-4) ** 2


def read_gen(d):
    import tifffile
    p = Path(d) / "generations.tif"
    return tifffile.imread(p).astype(float) if p.exists() else None


def cmd_compare(a):
    g = K.load(a.grown)
    r = K.load(a.ref)
    gen = read_gen(a.grown)
    res = {"label": a.label, "grown": str(a.grown), "ref": str(a.ref), "gen_max": a.gen_max, "grown_grid": g["shape"], "ref_grid": r["shape"]}
    box = tuple(a.refbox) if a.refbox else None
    res.update(K.compare(g, r, gen=gen, gen_max=a.gen_max, ref_nodebox=box))
    Path(a.out).write_text(json.dumps(res, indent=1))
    print(json.dumps({k: res[k] for k in ("label", "grown_nodes", "comparable_nodes", "distance", "orientation")}, indent=1))
    print("null:", json.dumps(res["null"]["random_shift_60_within_8"]), {k: v["within_8"] for k, v in res["null"].items() if k.startswith("normal") and v})


def cmd_cleanref(a):
    """how much of the grown mesh's largest clean rectangle (from an audit json) lies on the reference surface and inside the reference node box"""
    from scipy.spatial import cKDTree
    au = json.loads(Path(a.audit).read_text())
    lr = au["largest_clean_rectangle"]
    g = K.load(a.grown)
    r = K.load(a.ref)
    rp, rn, rint, rrc = K.fine_lattice(r)
    (r0, r1), (c0, c1) = lr["rows"], lr["cols"]
    P = g["points"][r0:r1, c0:c1].reshape(-1, 3)
    d, idx = cKDTree(rp).query(P)
    rc = rrc[idx]
    b = a.refbox
    inbox = (rc[:, 0] >= b[0]) & (rc[:, 0] <= b[1] - 1) & (rc[:, 1] >= b[2]) & (rc[:, 1] <= b[3] - 1)
    res = {"label": a.label, "clean_rectangle": lr, "nodes": int(len(P)), "share_within_8_of_ref": float((d <= K.NEAR_VOX).mean()), "share_in_refbox": float(inbox.mean()),
           "share_in_refbox_and_within_8": float((inbox & (d <= K.NEAR_VOX)).mean()), "median_dist": float(np.median(d))}
    Path(a.out).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


def cmd_audit(a):
    gm = K.load(a.mesh)
    med = K.audit261.edge_stretch(gm["points"])[1]
    box = [float(v) for v in a.box_xyz] if a.box_xyz else None
    res, _, _ = K.audit(a.mesh, SHAPES[a.scroll], a.label, cell_cm2(med), box_xyz=box, box_margin=a.box_margin)
    res["mesh"] = str(a.mesh)
    res["box_xyz"], res["box_margin"] = box, a.box_margin
    Path(a.out).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


def cmd_refs(a):
    idx = json.loads((K.C.FAC / "index.json").read_text())["surfaces"]
    names = sorted(s["surface"] for s in (idx if isinstance(idx, list) else idx.values()) if s["surface"].startswith("PHerc0826"))
    out = []
    for nm in names:
        d = K.C.FAC / "meshes" / (nm + ".tifxyz")
        gm = K.load(d)
        med = K.audit261.edge_stretch(gm["points"])[1]
        res, _, _ = K.audit(d, SHAPES["PHerc0826"], nm, cell_cm2(med))
        out.append(res)
        print(nm, res["largest_clean_rectangle"], "jump", round(res["jump_share"], 4), flush=True)
    Path(a.out).write_text(json.dumps(out, indent=1))
    areas = [r["largest_clean_rectangle"]["area_cm2"] for r in out]
    print("largest clean rectangle cm2 over", len(out), "official PHerc0826 meshes: median", float(np.median(areas)), "max", max(areas), "n>=4", int(sum(x >= 4 for x in areas)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compare")
    c.add_argument("--grown", required=True)
    c.add_argument("--ref", required=True)
    c.add_argument("--label", required=True)
    c.add_argument("--out", required=True)
    c.add_argument("--gen-max", type=float, default=None)
    c.add_argument("--refbox", type=int, nargs=4, default=None)
    c = sub.add_parser("cleanref")
    c.add_argument("--grown", required=True)
    c.add_argument("--ref", required=True)
    c.add_argument("--audit", required=True)
    c.add_argument("--refbox", type=int, nargs=4, required=True)
    c.add_argument("--label", required=True)
    c.add_argument("--out", required=True)
    c = sub.add_parser("audit")
    c.add_argument("--mesh", required=True)
    c.add_argument("--scroll", required=True)
    c.add_argument("--box-xyz", type=float, nargs=6, default=None)
    c.add_argument("--box-margin", type=float, default=0.0)
    c.add_argument("--label", required=True)
    c.add_argument("--out", required=True)
    c = sub.add_parser("refs")
    c.add_argument("--out", required=True)
    args = ap.parse_args()
    {"compare": cmd_compare, "audit": cmd_audit, "refs": cmd_refs, "cleanref": cmd_cleanref}[args.cmd](args)
