"""T-269 (X-213) mirror-block choice by X-202's on-sheet ranking (RULE fixed in X-213 before any tracer run; geometry of stored meshes only, no tracer output is used).
ANCHOR = the stored mesh(es) of a scroll that X-202 (data/picture-factory/onsheet_d/ranking.json, rank_area = on-sheet cm2) ranks highest.  BLOCK = anchor bbox (finite positive node coordinates, x y z) widened by 700 voxels on
every side (the seed wall margin); a side shorter than 3456 is widened symmetrically to 3456 (the X-212 cube side); a side longer than 5760 is cut to 5760 around the anchor centre; then rounded outward to the 192-voxel
chunk grid and clipped to the volume.
 PHerc0813 : the highest-ranked community (auto_grown) mesh;   PHerc0358 : same;   PHerc1545 : the highest-ranked mesh (it has no auto_grown mesh in the factory queue).
 PHerc0826 (box 1, existing X-210 mirror, no fetch): box [2894, 8099, 3539, 6217, 10985, 16275] unchanged, continued from X-212's state.
 PHerc0813D (control, independent of PHerc0813): the X-212 densest-cube rule, written by fetch269.py block.
 PHerc0826B (second box): PHerc0826 surfaces with < 50 pct of their nodes inside box 1, taken in rank order; the first is the anchor, each next one is added iff the union bbox widened by 700 stays <= 5760 on every side.
python sites269.py   -> data/t269-trace/sites.json"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import tifffile

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common269 as C  # noqa: E402

MARGIN, MINSIDE, MAXSIDE, CH = 700, 3456, 5760, 192
BOX1 = [2894, 8099, 3539, 6217, 10985, 16275]


def bbox(mesh_dir):
    P = np.stack([tifffile.imread(Path(mesh_dir) / f"{k}.tif").astype(np.float64) for k in "xyz"], -1).reshape(-1, 3)
    P = P[np.isfinite(P).all(1) & (P > 0).all(1)]
    return P, P.min(0), P.max(0)


def block(lo, hi, shape_zyx):
    out = []
    for i in range(3):
        a, b = lo[i] - MARGIN, hi[i] + MARGIN
        c = (lo[i] + hi[i]) / 2
        if b - a < MINSIDE:
            a, b = c - MINSIDE / 2, c + MINSIDE / 2
        if b - a > MAXSIDE:
            a, b = c - MAXSIDE / 2, c + MAXSIDE / 2
        lim = shape_zyx[2 - i]                        # x -> shape[2], y -> shape[1], z -> shape[0]
        a, b = int(np.floor(a / CH) * CH), int(np.ceil(b / CH) * CH)
        out += [max(0, a), min(lim, b)]
    return out


def main():
    rk = json.loads((C.FAC / "onsheet_d/ranking.json").read_text())["rows"]
    q = json.loads((C.FAC / "queue.json").read_text())
    surf = {s["name"]: s for s in q["surfaces"]}
    rows = sorted(rk, key=lambda r: r["rank_area"])
    cache = {}

    def get(name):
        if name not in cache:
            cache[name] = bbox(surf[name]["mesh_dir"])
        return cache[name]

    sites = {}
    sites["PHerc0826"] = {"scroll": "PHerc0826", "mirror": str(C.DATA267 / "PHerc0826"), "docker_volume": "t267-0826", "box_xyz": BOX1, "prior_state_dirs": [str(C.DATA268 / "PHerc0826")],
                          "anchor": "PHerc0826_auto_grown_20260902140841364 (X-210 mirror; box = its bbox + 700, unchanged)", "fetch": False, "target_traces": 64, "order": 0}
    vol = C.queue_volumes()
    inside = {}
    for r in rows:
        if r["scroll"] != "PHerc0826":
            continue
        P, lo, hi = get(r["surface"])
        inside[r["surface"]] = float(((P >= [BOX1[0], BOX1[2], BOX1[4]]) & (P <= [BOX1[1], BOX1[3], BOX1[5]])).all(1).mean())
    cand = [r for r in rows if r["scroll"] == "PHerc0826" and inside[r["surface"]] < 0.5]
    anchors = []
    ulo = uhi = None
    for r in cand:
        P, lo, hi = get(r["surface"])
        nlo, nhi = (lo, hi) if ulo is None else (np.minimum(ulo, lo), np.maximum(uhi, hi))
        if (nhi - nlo + 2 * MARGIN).max() <= MAXSIDE:
            anchors.append({"surface": r["surface"], "rank_area": r["rank_area"], "onsheet_cm2": r["area_cm2"], "share": r["share"], "inside_box1": inside[r["surface"]],
                            "bbox_lo": lo.round().tolist(), "bbox_hi": hi.round().tolist()})
            ulo, uhi = nlo, nhi
    shape = tuple(vol["PHerc0826"]["shape"])
    sites["PHerc0826B"] = {"scroll": "PHerc0826", "mirror": str(C.DATA / "PHerc0826B"), "docker_volume": "t269-pherc0826b", "box_xyz": block(ulo, uhi, shape), "prior_state_dirs": [], "anchors": anchors,
                           "fetch": True, "target_traces": 32, "order": 3}
    for sc, pref, order in (("PHerc0813", "auto_grown", 1), ("PHerc0358", "auto_grown", 2), ("PHerc1545", "", 4)):
        top = [r for r in rows if r["scroll"] == sc and pref in r["surface"]][0]
        P, lo, hi = get(top["surface"])
        shape = tuple(vol[sc]["shape"])
        sites[sc] = {"scroll": sc, "mirror": str(C.DATA / sc), "docker_volume": f"t269-{sc.lower()}", "box_xyz": block(lo, hi, shape), "prior_state_dirs": [],
                     "anchors": [{"surface": top["surface"], "rank_area": top["rank_area"], "onsheet_cm2": top["area_cm2"], "share": top["share"], "bbox_lo": lo.round().tolist(), "bbox_hi": hi.round().tolist()}],
                     "fetch": True, "target_traces": 16, "order": order}
    sites["PHerc0813D"] = {"scroll": "PHerc0813", "mirror": str(C.DATA / "PHerc0813D"), "docker_volume": "t269-pherc0813d", "box_xyz": None, "prior_state_dirs": [], "independent": True, "control": True,
                           "anchor": "none: the X-212 densest level-3 cube rule (fetch269.py block), box_xyz written there", "fetch": True, "target_traces": 16, "order": 5}
    for s in sites.values():
        b = s["box_xyz"]
        s["box_sides"] = [b[1] - b[0], b[3] - b[2], b[5] - b[4]] if b else None
    C.jwrite(C.SITES_JSON, {"made": time.strftime("%Y-%m-%d %H:%M:%S"), "rule": __doc__, "margin": MARGIN, "min_side": MINSIDE, "max_side": MAXSIDE, "sites": sites, "box1_inside_share": inside})
    for k, s in sites.items():
        print(k, s["scroll"], s["box_xyz"], s["box_sides"], [a["surface"] for a in s.get("anchors", [])])


if __name__ == "__main__":
    main()
