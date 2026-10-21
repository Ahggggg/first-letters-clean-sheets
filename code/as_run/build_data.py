"""build_data.py (T-276): gather every kept clean patch of X-211 (t267), X-212 (t268), X-213 (t269) and any finished T-275 patch, copy the cropped tifxyz into the staging repository with a clean
meta.json (no local paths), recompute the rectangle-level checks on the cropped files, and write data/t276-xsec/records.json (the single table every later step reads).
CPU only, no network.  python build_data.py [--stage DIR]"""
import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import tifffile

HERE = Path(__file__).resolve().parent
PROB = HERE.parents[1]
STAGE = Path("<path of this repository>")
sys.path.insert(0, str(STAGE / "code"))
import tifxyz_audit as A  # noqa: E402

FAC = PROB / "data/picture-factory"
VOLS = json.loads((FAC / "queue.json").read_text())["volumes"]
IMAGE_DIGEST = "sha256:bad516f66001abca759454cc43e4fd11e5b19aa55d36bdc2043817291c8083c4"


def jl(p):
    return json.loads(Path(p).read_text())


def all_traces(tn):
    out = {}
    for f in sorted((PROB / f"data/{tn}-trace").glob("*/state.json")):
        out.update(jl(f)["traces"])
    return out


def collect():
    recs = []
    # ---- t267 B1 (X-210 / X-211)
    idx = jl(PROB / "data/t267-trace/render/index.json")["surfaces"][0]
    aud = jl(PROB / "data/t267-trace/audit-B1.json")
    sup = jl(PROB / "data/t267-trace/support-B1-clean.json")
    mesh_meta = jl(PROB / "data/t267-trace/render/meshes/PHerc0826_t267_B1_clean.tifxyz/meta.json")
    summ = None
    recs.append({"name": idx["surface"], "scroll": "PHerc0826", "run": "t267", "experiment": "X-210 (trace) / X-211 (read)", "task": "T-267", "mesh": PROB / "data/t267-trace/render/meshes/PHerc0826_t267_B1_clean.tifxyz",
                 "round": None, "tid": "B1_seed75_75_G110", "seed": mesh_meta["seed"], "rows": [2, 217], "cols": [2, 113], "nodes": 23865, "cm2": idx["cm2"], "side": -1, "frame_correct": "R",
                 "median_edge_vox": aud["median_edge_vox"], "trace_frame_minority": aud["frame"]["minority"], "trace_fold_share": aud["jump_share"], "trace_stretch_share": aud["stretch_share_nodes"],
                 "generations": 110, "tracer_seconds": mesh_meta["elapsed_time_s"], "wall_outside_share": None, "on_prediction_share": sup["on_prediction_share"],
                 "nonfinite": idx["nonfinite_valid_px"], "picture": PROB / "data/t267-trace/render/png/PHerc0826_t267_B1_clean.png", "maps": PROB / "data/t267-trace/render/maps/PHerc0826_t267_B1_clean.npz",
                 "stack": PROB / "data/t267-trace/render/stacks/PHerc0826_t267_B1_clean.npy", "um_per_canvas_px": idx["um_per_canvas_px"], "site": "PHerc0826", "read_summary": None})
    # ---- t268, t269
    for tn, exp, task in (("t268", "X-212", "T-268"), ("t269", "X-213", "T-269")):
        d = jl(PROB / f"data/{tn}-trace/index.json")
        traces = all_traces(tn)
        sup = {}
        for f in (PROB / f"data/{tn}-trace/support.json", PROB / "data/t276-xsec/support_t269.json"):
            if f.exists() and (tn == "t268" or "t269" in f.name):
                sup.update(jl(f))
        for p in d["patches"]:
            t = traces.get(p["tid"], {})
            recs.append({"name": p["name"], "scroll": p["scroll"], "run": tn, "experiment": exp, "task": task, "mesh": PROB / f"data/{tn}-trace/render/meshes/{p['name']}.tifxyz", "round": p["round"], "tid": p["tid"],
                         "seed": p["seed"], "rows": p["rows"], "cols": p["cols"], "nodes": p["nodes"], "cm2": p["cm2"], "side": p["side"], "frame_correct": p["frame_correct"], "median_edge_vox": p["median_edge_vox"],
                         "trace_frame_minority": p["trace_frame_minority"], "trace_fold_share": p["trace_fold_share"], "trace_stretch_share": p["trace_stretch_share"], "generations": t.get("max_gen", 75),
                         "tracer_seconds": t.get("tracer_seconds"), "wall_outside_share": t.get("wall_outside_share"), "on_prediction_share": sup.get(p["name"], {}).get("on_prediction_share"),
                         "nonfinite": p["nonfinite_valid_px"], "picture": PROB / f"outputs/{tn}-trace/{p['name']}.png", "maps": Path(p["maps_npz"]), "stack": PROB / f"data/{tn}-trace/render/stacks/{p['name']}.npy",
                         "um_per_canvas_px": p["um_per_canvas_px"], "site": p.get("site", p["scroll"]), "read_summary": p["summary"], "share_within_8_of_read_ground": p.get("share_within_8_of_read_ground")})
    # ---- T-275 patches, if any finished (index at data/t275-trace/index.json with the same layout)
    f275 = PROB / "data/t275-trace/index.json"
    if f275.exists():
        d = jl(f275)
        traces = all_traces("t275")
        sup275 = jl(PROB / "data/t275-trace/support.json") if (PROB / "data/t275-trace/support.json").exists() else {}
        for p in d["patches"]:
            if not (PROB / f"data/t275-trace/render/maps/{p['name']}.npz").exists():
                continue
            t = traces.get(p["tid"], {})
            recs.append({"name": p["name"], "scroll": p["scroll"], "run": "t275", "experiment": "X-217", "task": "T-275", "mesh": PROB / f"data/t275-trace/render/meshes/{p['name']}.tifxyz", "round": p["round"],
                         "tid": p["tid"], "seed": p["seed"], "rows": p["rows"], "cols": p["cols"], "nodes": p["nodes"], "cm2": p["cm2"], "side": p["side"], "frame_correct": p["frame_correct"],
                         "median_edge_vox": p["median_edge_vox"], "trace_frame_minority": p["trace_frame_minority"], "trace_fold_share": p["trace_fold_share"], "trace_stretch_share": p["trace_stretch_share"],
                         "generations": t.get("max_gen", 75), "tracer_seconds": t.get("tracer_seconds"), "wall_outside_share": t.get("wall_outside_share"), "on_prediction_share": sup275.get(p["name"], {}).get("on_prediction_share"),
                         "nonfinite": p["nonfinite_valid_px"], "picture": PROB / f"outputs/t275-trace/{p['name']}.png", "maps": Path(p["maps_npz"]), "stack": PROB / f"data/t275-trace/render/stacks/{p['name']}.npy",
                         "um_per_canvas_px": p["um_per_canvas_px"], "site": p.get("site", p["scroll"]), "read_summary": p["summary"]})
    return recs


def write_patch(r, stage):
    v = VOLS[r["scroll"]]
    src = Path(r["mesh"])
    dst = stage / "data" / r["scroll"] / f"{r['name']}.tifxyz"
    dst.mkdir(parents=True, exist_ok=True)
    for a in "xyz":
        shutil.copyfile(src / f"{a}.tif", dst / f"{a}.tif")
    shutil.copyfile(src / "meta.json", dst / "meta.json")          # provisional (overwritten below)
    om = jl(src / "meta.json")
    xyz = np.stack([tifffile.imread(dst / f"{a}.tif").astype(np.float64) for a in "xyz"], -1)
    ok = np.isfinite(xyz).all(-1) & (xyz > 0).all(-1)
    lo, hi = xyz[ok].min(0), xyz[ok].max(0)
    voxel_um = v["voxel_um"]
    a = A.audit(dst, tuple(v["shape"]), voxel_um)
    cm2 = float(r["cm2"])
    meta = {"format": "tifxyz", "type": "seg", "uuid": r["name"], "scale": om["scale"],
            "bbox": [[float(lo[0]), float(lo[1]), float(lo[2])], [float(hi[0]), float(hi[1]), float(hi[2])]],
            "area_cm2": cm2, "area_vx2": cm2 / (voxel_um * 1e-4) ** 2, "voxel_size_um": voxel_um, "voxelsize": voxel_um,
            "scroll": r["scroll"], "volume_id": v["volume"], "volume_shape_zyx": list(v["shape"]),
            "source": "vc_grow_seg_from_seed", "vc_gsfs_version": om.get("vc_gsfs_version"), "seed": om["seed"], "target_volume": om.get("target_volume"), "max_gen": om.get("max_gen"),
            "vc_gsfs_params": {k: om["vc_gsfs_params"][k] for k in ("generations", "min_area_cm", "mode", "thread_limit", "use_cuda", "voxelsize")},
            "elapsed_time_s_whole_trace": om.get("elapsed_time_s"),
            "clean_patch": {"crop_rows_of_trace_grid_half_open": om["crop_rows"], "crop_cols_of_trace_grid_half_open": om["crop_cols"], "nodes": r["nodes"], "frame_side": r["side"],
                            "frame_correct_order": r["frame_correct"], "node_spacing_vox_median": r["median_edge_vox"] or a["median_edge_vox"], "area_cm2_note": "nodes x (median edge x voxel size)^2, the area quoted in index.csv",
                            "source_experiment": r["experiment"], "trace_id": r["tid"], "round": r["round"]},
            "note": "Rectangle cut from a larger trace made with villa's vc_grow_seg_from_seed on the public m7 surface prediction; all nodes pass the mesh-level clean rules of code/tifxyz_audit.py.  Mesh only: no check beyond the mesh and cross-sections in qa/.",
            "license": "CC BY-NC 4.0 (derived from Vesuvius Challenge scan data and prediction, which are CC BY-NC 4.0)"}
    (dst / "meta.json").write_text(json.dumps(meta, indent=1))
    return dst, a, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default=str(STAGE))
    a = ap.parse_args()
    stage = Path(a.stage)
    recs = collect()
    out = []
    for r in recs:
        dst, aud, meta = write_patch(r, stage)
        rr = {k: (str(v) if isinstance(v, Path) else v) for k, v in r.items()}
        rr["stage_mesh"] = str(dst)
        rr["rect_audit"] = aud
        rr["volume_id"] = VOLS[r["scroll"]]["volume"]
        rr["voxel_um"] = VOLS[r["scroll"]]["voxel_um"]
        out.append(rr)
        print(r["name"], r["cm2"], "rect minority", round(aud["frame"]["minority"], 4), "fold", None if aud["fold_share"] is None else round(aud["fold_share"], 4), "max jump", round(aud["max_jump_deg"], 1),
              "all clean", aud["all_nodes_clean"], "tri/area", round(aud["area_cm2_triangulated"] / r["cm2"], 3))
    (PROB / "data/t276-xsec").mkdir(parents=True, exist_ok=True)
    (PROB / "data/t276-xsec/records.json").write_text(json.dumps(out, indent=1, default=str))
    print(len(out), "patches", round(sum(x["cm2"] for x in out), 2), "cm2")


if __name__ == "__main__":
    main()
