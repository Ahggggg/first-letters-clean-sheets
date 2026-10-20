"""build_index.py (T-276): index.csv + xsec/ tables of the staging repository from records.json (build_data.py), X-218 summary and the post-hoc support numbers.  CPU only."""
import csv
import json
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROB = HERE.parents[1]
STAGE = Path("<path of this repository>")
OUT = PROB / "data/t276-xsec"

recs = json.loads((OUT / "records.json").read_text())
summ = json.loads((OUT / "summary.json").read_text())
pp = {p["name"]: p for p in summ["per_patch"]}
X212_OFF = {"PHerc0826_t268_s004", "PHerc0826_t268_s014", "PHerc0826_t268_s033"}

# B1 ink-map shares from the stored maps (X-211 stored numbers: R 0.032 / F 0.076)
for r in recs:
    if r["read_summary"] is None:
        z = np.load(r["maps"])
        v = z["valid"].astype(bool)
        sm = {}
        for arm in ("F", "R"):
            m = np.minimum(z["rv2__" + arm].astype(np.float32), z["d9v4C__" + arm].astype(np.float32))[v]
            sm[arm] = {"share_gt_0.5": float((m > 0.5).mean()), "share_gt_0.8": float((m > 0.8).mean())}
        r["read_summary"] = sm

cols = ["scroll", "patch", "volume_id", "voxel_um", "seed_x", "seed_y", "seed_z", "area_cm2", "nodes", "nodes_rows", "nodes_cols", "frame_side", "frame_correct_order", "rect_frame_minority",
        "rect_fold_share_interior", "rect_stretch_share", "trace_frame_minority", "trace_fold_share", "trace_stretch_share", "source_experiment", "source_task", "trace_id", "round", "tracer_generations",
        "tracer_seconds", "on_prediction_share", "xsec_median_abs_tilt_deg", "xsec_median_abs_drift_vox", "xsec_share_sections_drift_ge_20vox", "ink_share_gt_0.5_frame_correct", "ink_share_gt_0.8_frame_correct",
        "ink_nonfinite_px", "trace_wall_outside_share", "flags", "tifxyz", "qa_picture", "ink_maps_both_orders"]
rows = []
for r in recs:
    ra = r["rect_audit"]
    p = pp[r["name"]]
    fc = r["frame_correct"]
    flags = []
    if r["name"] in X212_OFF or (r["on_prediction_share"] is not None and r["on_prediction_share"] < 0.30):
        flags.append("off_public_prediction")      # X-212's three named patches (0.17-0.41) or, post-hoc, any patch with < 0.30 of its nodes on the public m7 prediction (displaced-node null about 0.24)
    if sum(r["nonfinite"].values()) > 0:
        flags.append("ink_map_nonfinite_px")
    if (r.get("wall_outside_share") or 0) > 0.2:
        flags.append("trace_touches_mirror_box_wall")
    if p["median_abs_tilt"] is not None and p["median_abs_tilt"] >= 2.0:
        flags.append("xsec_tilt_ge_2deg_post_hoc")
    rows.append({"scroll": r["scroll"], "patch": r["name"], "volume_id": r["volume_id"], "voxel_um": r["voxel_um"], "seed_x": int(r["seed"][0]), "seed_y": int(r["seed"][1]), "seed_z": int(r["seed"][2]),
                 "area_cm2": round(r["cm2"], 3), "nodes": r["nodes"], "nodes_rows": ra["grid"][0], "nodes_cols": ra["grid"][1], "frame_side": r["side"], "frame_correct_order": fc,
                 "rect_frame_minority": round(ra["frame"]["minority"], 5), "rect_fold_share_interior": round(ra["fold_share_interior"], 5), "rect_stretch_share": round(ra["stretch_share"], 5),
                 "trace_frame_minority": round(r["trace_frame_minority"], 4), "trace_fold_share": round(r["trace_fold_share"], 4), "trace_stretch_share": "" if r["trace_stretch_share"] is None else round(r["trace_stretch_share"], 4),
                 "source_experiment": r["experiment"], "source_task": r["task"], "trace_id": r["tid"], "round": "" if r["round"] is None else r["round"], "tracer_generations": r["generations"],
                 "tracer_seconds": "" if r["tracer_seconds"] is None else round(r["tracer_seconds"]), "on_prediction_share": "" if r["on_prediction_share"] is None else round(r["on_prediction_share"], 3),
                 "xsec_median_abs_tilt_deg": round(p["median_abs_tilt"], 2), "xsec_median_abs_drift_vox": round(p["median_abs_drift"], 1), "xsec_share_sections_drift_ge_20vox": round(p["share_drift_ge_20"], 2),
                 "ink_share_gt_0.5_frame_correct": round(r["read_summary"][fc]["share_gt_0.5"], 4), "ink_share_gt_0.8_frame_correct": round(r["read_summary"][fc]["share_gt_0.8"], 4),
                 "ink_nonfinite_px": sum(r["nonfinite"].values()) // 4 if r["nonfinite"].get("rv2__F") is not None else 0,
                 "trace_wall_outside_share": "" if r.get("wall_outside_share") is None else round(r["wall_outside_share"], 3), "flags": ";".join(flags),
                 "tifxyz": f"data/{r['scroll']}/{r['name']}.tifxyz", "qa_picture": f"qa/{r['scroll']}/{r['name']}_qa.png", "ink_maps_both_orders": f"qa/{r['scroll']}/{r['name']}_both_orders.png"})
rows.sort(key=lambda x: (x["scroll"], x["patch"]))
with open(STAGE / "index.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    w.writerows(rows)
print(len(rows), "rows;", round(sum(x["area_cm2"] for x in rows), 2), "cm2")
by = {}
for x in rows:
    b = by.setdefault(x["scroll"], [0, 0.0])
    b[0] += 1
    b[1] += x["area_cm2"]
print({k: (v[0], round(v[1], 2)) for k, v in by.items()})
# xsec tables
(STAGE / "xsec").mkdir(exist_ok=True)
for n in ("summary.json", "g1.json", "identity.json"):
    shutil.copyfile(OUT / n, STAGE / "xsec" / n)
with open(STAGE / "xsec/patch_sections.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["patch", "scroll", "kind", "canvas_index", "frac", "valid_cols", "arclength_vox", "tilt_deg", "coherence", "noise_threshold", "informative", "implied_drift_vox", "uncached_cols"])
    for r in recs:
        pj = json.loads((OUT / "patches" / f"{r['name']}.json").read_text())
        for s in pj["sections"]:
            w.writerow([r["name"], r["scroll"], s["kind"], s["index"], s["frac"], s["n_valid_cols"], round(s["arclen_vox"], 1), round(s["tilt"], 3), round(s["coherence"], 4), round(s["threshold"], 4),
                        int(s["informative"]), round(s["drift_vox"], 1), s["unknown_cols"]])
with open(STAGE / "xsec/control_sections.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["surface", "group", "kind", "canvas_index", "span_start", "span_end", "tilt_deg", "coherence", "noise_threshold", "informative", "implied_drift_vox"])
    for fn in sorted((OUT / "controls").glob("*.json")):
        c = json.loads(fn.read_text())
        for s in c["sections"]:
            if "tilt" in s:
                w.writerow([c["name"], c["group"], s["kind"], s["index"], s["span"][0], s["span"][1], round(s["tilt"], 3), round(s["coherence"], 4), round(s["threshold"], 4), int(s["informative"]), round(s["drift_vox"], 1)])
