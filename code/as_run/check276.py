"""check276.py (T-276): checks on the staged patches, written to the staging repository's checks/ folder.  (1) cut_patches.py reproduces the stored rectangle of the first patch of each T-269 site;
(2) ground counted once: nodes of each patch within 12 voxels of a node of another patch of the same scroll; (3) every tifxyz loads with finite nodes; sizes.  CPU only."""
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROB = HERE.parents[1]
STAGE = Path("<path of this repository>")
sys.path.insert(0, str(STAGE / "code"))
sys.path.insert(0, str(STAGE / "recipe"))
import tifxyz_audit as A  # noqa: E402
import cut_patches as CP  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

(STAGE / "checks").mkdir(exist_ok=True)
recs = json.loads((PROB / "data/t276-xsec/records.json").read_text())
vols = json.loads((PROB / "data/picture-factory/queue.json").read_text())["volumes"]

# (1) reproduction of the cut on stored traces (first patch of round 1 of each T-269 site; no exclusion needed, no box used)
d = json.loads((PROB / "data/t269-trace/index.json").read_text())
rows = []
seen = set()
for p in d["patches"]:
    site = p["site"]
    if site in seen or site == "PHerc0826":
        continue
    seen.add(site)
    q = min((x for x in d["patches"] if x["site"] == site and x["round"] == 1), key=lambda x: x["tid"])
    seg = PROB / f"data/t269-trace/{site}/runs/r001/out/{q['tid']}/seg"
    r, _, _ = CP.cut(seg, tuple(vols[q["scroll"]]["shape"]), 9.362, 1.0)
    rows.append({"patch": q["name"], "stored_rows": q["rows"], "stored_cols": q["cols"], "stored_nodes": q["nodes"], "recut_rows": r["rows"], "recut_cols": r["cols"], "recut_nodes": r["nodes"],
                 "identical": bool(q["rows"] == r["rows"] and q["cols"] == r["cols"] and q["nodes"] == r["nodes"] and q["side"] == r["side"])})
(STAGE / "checks/cut_patches_reproduction.json").write_text(json.dumps(rows, indent=1))
print("REPRO", json.dumps(rows))

# (2) ground counted once
by = {}
for r in recs:
    P, _, _ = A.load(Path(r["stage_mesh"]))
    by.setdefault(r["scroll"], []).append((r["name"], P[np.isfinite(P).all(-1)]))
out = []
for sc, lst in by.items():
    for i, (name, P) in enumerate(lst):
        others = [Q for j, (n, Q) in enumerate(lst) if j != i]
        if not others:
            out.append((sc, name, len(P), 0, 0.0))
            continue
        tree = cKDTree(np.concatenate(others))
        dd, _ = tree.query(P)
        out.append((sc, name, len(P), int((dd <= 12.0).sum()), float((dd <= 12.0).mean())))
with open(STAGE / "checks/overlap_between_patches.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["scroll", "patch", "nodes", "nodes_within_12_vox_of_another_patch", "share"])
    for o in out:
        w.writerow([o[0], o[1], o[2], o[3], round(o[4], 4)])
sh = np.array([o[4] for o in out])
print("OVERLAP share >0:", int((sh > 0).sum()), "max", sh.max().round(4), "mean", sh.mean().round(5), "total nodes within 12:", sum(o[3] for o in out), "of", sum(o[2] for o in out))
# (3) integrity
bad = []
tot = 0
for r in recs:
    P, N, m = A.load(Path(r["stage_mesh"]))
    if not np.isfinite(P).all():
        bad.append(r["name"])
    tot += int(np.isfinite(P).all(-1).sum())
print("INTEGRITY non-finite nodes in", bad, "total nodes", tot)
(STAGE / "checks/summary.json").write_text(json.dumps({"patches": len(recs), "patches_with_nonfinite_nodes": bad, "total_nodes": tot, "reproduction_identical": sum(1 for r in rows if r["identical"]), "reproduction_n": len(rows),
                                                      "overlap_patches_with_any_overlap": int((sh > 0).sum()), "overlap_max_share": float(sh.max()), "overlap_nodes": sum(o[3] for o in out)}, indent=1))
