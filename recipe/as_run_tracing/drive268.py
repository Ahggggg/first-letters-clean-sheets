"""T-268 (X-212) orchestrator.  ONE process (run it under tools/run.py --experiment X-212, hidden through tools/bg.py).  Per scroll, in rounds:
  1. choose up to 4 seeds from data/t268-trace/<scroll>/seed-candidates.json (flattest first), skipping seeds within 16 voxels of a node of an earlier whole trace (same sheet already grown),
     and keeping the seeds of a round >= 1000 voxels apart;
  2. ONE detached Docker container (image volume-cartographer:main, entry point vc_grow_seg_from_seed, params = the official recipe with generations 75, thread_limit 1, cache 3e9, no network, mirror in a Docker
     named volume) runs the 4 tracers side by side (1 thread each);
  3. audit of each trace (audit268.py): frame minority, fold share, stretch share, wall; the largest clean rectangle (a) of the trace alone and (b) with nodes within 12 voxels of an earlier kept patch (and of the
     X-210 B1 clean rectangle on PHerc0826) removed;  (b) >= 1.0 cm2 -> kept, cut out as a tifxyz crop and queued;
  4. kept patches rendered (project renderer, 28 planes) and read (Reader v2 + d9v4C, both orders, Arc FP32, arc.lock by exclusive create) by read268.py;
  5. pictures (40 um/px, 1 cm bar, frame-correct order marked), contact sheets, reading summaries, index data/t268-trace/index.json.
Resumable: state per scroll in data/t268-trace/<scroll>/state.json.  Stop between rounds by creating data/t268-trace/STOP.
python drive268.py --scroll PHerc0826 [--max-traces 32] [--rounds N]"""
import argparse
import json
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common268 as C  # noqa: E402

os.environ["PATH"] = C.DOCKER_BIN.replace("/", "\\") + os.pathsep + os.environ.get("PATH", "")
PY = sys.executable
INDEX = C.DATA / "index.json"
STOP = C.DATA / "STOP"


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def sh(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


# ------------------------------------------------------------------------------------------------ state
def load_state(sc):
    p = C.scroll_dir(sc) / "state.json"
    return json.loads(p.read_text()) if p.exists() else {"scroll": sc, "traces": {}, "rounds": [], "patches": []}


def save_state(sc, st):
    C.jwrite(C.scroll_dir(sc) / "state.json", st)


def load_index():
    return json.loads(INDEX.read_text()) if INDEX.exists() else {"experiment": "X-212", "task": "T-268", "patches": [], "traces": [], "format": __doc__.split("Resumable")[0]}


# ------------------------------------------------------------------------------------------------ seeds
def traced_nodes(sc, st):
    pts = []
    if sc == "PHerc0826":
        import tifffile
        for r in ("B1_seed75_75_G110", "B2_seed40_100_G75"):
            d = C.DATA267 / "runs" / r / "out" / "seg"
            a = [tifffile.imread(d / f"{k}.tif").astype(np.float32) for k in "xyz"]
            p = np.stack(a, -1).reshape(-1, 3)
            pts.append(p[np.isfinite(p).all(-1) & (p > 0).all(-1)])
    for tid in st["traces"]:
        f = C.scroll_dir(sc) / "traces" / f"{tid}.nodes.npy"
        if f.exists():
            pts.append(np.load(f))
    return np.concatenate(pts) if pts else np.zeros((0, 3), np.float32)


def pick_seeds(sc, st, cands, n):
    from scipy.spatial import cKDTree
    tn = traced_nodes(sc, st)
    tree = cKDTree(tn) if len(tn) else None
    chosen = []
    for c in cands:
        tid = f"{sc}_s{c['no']:03d}"
        if tid in st["traces"]:
            continue
        s = np.array(c["seed"], float)
        if tree is not None and tree.query(s)[0] <= C.SKIP_TRACED_VOX:
            st["traces"].setdefault(tid, {"status": "skipped_on_traced_sheet", "seed": c["seed"], "no": c["no"]})
            continue
        if any(np.linalg.norm(s - np.array(o["seed"], float)) < C.ROUND_SEP_VOX for o in chosen):
            continue
        chosen.append(c)
        if len(chosen) == n:
            break
    return chosen


# ------------------------------------------------------------------------------------------------ tracer container
def run_round_container(sc, rnd, chosen):
    vi = C.vol_info(sc)
    rdir = C.scroll_dir(sc) / "runs" / f"r{rnd:03d}"
    (rdir / "out").mkdir(parents=True, exist_ok=True)
    params = {"cache_root": "/tmp/cache", "cache_size": C.CACHE_SIZE, "generations": C.GENERATIONS, "min_area_cm": 0.3, "mode": "seed", "normal_grid_path": "/data/ngrid", "thread_limit": 1,
              "use_cuda": False, "voxelsize": vi["voxel_um"]}
    (rdir / "params.json").write_bytes(json.dumps(params).encode())
    lines = []
    for c in chosen:
        tid = f"{sc}_s{c['no']:03d}"
        (rdir / "out" / tid).mkdir(exist_ok=True)
        lines.append(f"{tid} {c['seed'][0]} {c['seed'][1]} {c['seed'][2]}")
    (rdir / "seeds.txt").write_bytes(("\n".join(lines) + "\n").encode())
    script = ("#!/bin/sh\nwhile read tid x y z; do\n"
              "  ( vc_grow_seg_from_seed -v /data/pred.zarr -t /out/$tid/seg -p /params.json --seed $x $y $z --segment-name seg > /out/$tid.log 2>&1 ; echo exit $? >> /out/$tid.log ) &\n"
              "done < /seeds.txt\nwait\n")
    (rdir / "run.sh").write_bytes(script.encode())
    name = f"t268-{sc.lower()}-r{rnd:03d}"
    sh(["docker", "rm", "-f", name])
    cmd = ["docker", "run", "-d", "--name", name, "--cpus", "5", "--memory", "14g", "-v", f"{C.docker_volume(sc)}:/data:ro", "-v", f"{rdir / 'params.json'}:/params.json:ro",
          "-v", f"{rdir / 'seeds.txt'}:/seeds.txt:ro", "-v", f"{rdir / 'run.sh'}:/run.sh:ro", "-v", f"{rdir / 'out'}:/out", "--entrypoint", "sh", C.IMAGE, "/run.sh"]
    r = sh(cmd)
    if r.returncode != 0:
        raise RuntimeError("docker run failed: " + r.stderr[-500:])
    t0 = time.time()
    log("container", name, "started with", len(chosen), "tracers", [f"s{c['no']}" for c in chosen])
    w = sh(["docker", "wait", name])
    secs = time.time() - t0
    insp = sh(["docker", "inspect", name, "--format", "{{.State.ExitCode}} {{.State.OOMKilled}}"]).stdout.strip()
    (rdir / "container.log").write_text(sh(["docker", "logs", name]).stdout[-4000:])
    sh(["docker", "rm", name])
    log("container", name, "exited", insp, "after", round(secs), "s")
    return rdir, secs, insp


# ------------------------------------------------------------------------------------------------ official coverage (already-read ground)
_REF = {}


def read_tree(sc):
    if sc in _REF:
        return _REF[sc]
    import audit268 as A
    from scipy.spatial import cKDTree
    q = json.loads((C.FAC / "queue.json").read_text())
    pts = []
    for s in q["surfaces"]:
        if s["scroll"] != sc or s["voxel_um"] != C.vol_info(sc)["voxel_um"]:
            continue
        try:
            gm = A.K.load(s["mesh_dir"])
            p, _, _, _ = A.K.fine_lattice(gm, 0.25)
            pts.append(p.astype(np.float32))
        except Exception as e:  # noqa: BLE001
            log("ref mesh failed", s["name"], type(e).__name__, e)
    t = cKDTree(np.concatenate(pts)) if pts else None
    _REF[sc] = (t, len(pts))
    return _REF[sc]


# ------------------------------------------------------------------------------------------------ audit of one round
def audit_round(sc, rnd, rdir, chosen, st):
    import audit268 as A
    import tifffile
    from scipy.spatial import cKDTree
    sys.path.insert(0, str(C.PROB / "work/t267-trace"))
    import read267
    vi = C.vol_info(sc)
    box = C.box_of(sc)
    kept_new = []
    for c in chosen:
        tid = f"{sc}_s{c['no']:03d}"
        rec = {"tid": tid, "no": c["no"], "seed": c["seed"], "round": rnd, "seed_r": c["r"], "seed_M0": c["M0"], "status": None}
        seg = rdir / "out" / tid / "seg"
        logtxt = (rdir / "out" / f"{tid}.log").read_text(errors="replace") if (rdir / "out" / f"{tid}.log").exists() else ""
        rec["tracer_exit_line"] = logtxt.strip().splitlines()[-1] if logtxt.strip() else None
        if not (seg / "x.tif").exists():
            rec["status"] = "no_output"
            rec["log_tail"] = logtxt[-600:]
            st["traces"][tid] = rec
            log(tid, "no output", rec["log_tail"][-200:].replace("\n", " | "))
            continue
        meta = json.loads((seg / "meta.json").read_text())
        rec.update({"tracer_area_cm2": meta.get("area_cm2"), "max_gen": meta.get("max_gen"), "tracer_seconds": meta.get("elapsed_time_s"), "bbox": meta.get("bbox")})
        gen = tifffile.imread(seg / "generations.tif") if (seg / "generations.tif").exists() else None
        rec["generations_reached"] = float(np.nanmax(gen)) if gen is not None else None
        stt = A.clean_state(seg, vi["shape"], box)
        rec["audit"] = A.summarize(stt, vi["voxel_um"])
        P = stt["P"]
        fin = np.isfinite(P).all(-1)
        (C.scroll_dir(sc) / "traces").mkdir(exist_ok=True)
        np.save(C.scroll_dir(sc) / "traces" / f"{tid}.nodes.npy", P[fin].astype(np.float32))
        alone = A.best_rect(stt)
        cell = rec["audit"]["cell_cm2"]
        rec["largest_clean_alone"] = {**alone, "cm2": alone["nodes"] * cell}
        kt = kept_tree(sc, st, kept_new)
        cov = A.covered_mask(P, kt)
        rec["covered_nodes"] = int((cov & fin).sum())
        excl = A.best_rect(stt, cov)
        rec["largest_clean_new"] = {**excl, "cm2": excl["nodes"] * cell}
        # wall: share of nodes outside the mirror box, and whether the trace reached its generation limit
        rec["wall_outside_share"] = rec["audit"]["outside_box_nodes"] / max(1, rec["audit"]["finite_nodes"])
        if excl["nodes"] * cell >= C.MIN_KEEP_CM2:
            r0, r1 = excl["rows"]
            c0, c1 = excl["cols"]
            name = f"{sc}_t268_{tid.split('_')[-1]}"
            read267.OUT = C.RENDER
            read267.crop(str(seg), r0, r1, c0, c1, name)
            rect = P[r0:r1, c0:c1].reshape(-1, 3).astype(np.float32)
            np.save(C.scroll_dir(sc) / "traces" / f"{tid}.kept.npy", rect)
            rt, nref = read_tree(sc)
            d_ref = rt.query(rect)[0] if rt is not None else None
            patch = {"name": name, "scroll": sc, "tid": tid, "round": rnd, "seed": c["seed"], "rows": excl["rows"], "cols": excl["cols"], "nodes": excl["nodes"], "cm2": excl["nodes"] * cell,
                     "side": excl["side"], "frame_correct": "F" if excl["side"] == 1 else "R", "median_edge_vox": rec["audit"]["median_edge_vox"], "voxel_um": vi["voxel_um"],
                     "trace_frame_minority": rec["audit"]["frame_minority"], "trace_fold_share": rec["audit"]["fold_share"], "trace_stretch_share": rec["audit"]["stretch_share"],
                     "share_within_8_of_read_ground": float((d_ref <= 8).mean()) if d_ref is not None else None, "n_read_meshes_ref": nref, "mesh": str(C.RENDER / "meshes" / f"{name}.tifxyz")}
            rec["status"] = "kept"
            rec["patch"] = name
            kept_new.append(patch)
            log(tid, "KEPT", f"{patch['cm2']:.2f} cm2", "side", excl["side"], "trace minority", round(rec["audit"]["frame_minority"], 3), "fold", round(rec["audit"]["fold_share"], 3))
        else:
            rec["status"] = "no_new_clean_rectangle"
            log(tid, f"not kept: new clean {excl['nodes'] * cell:.2f} cm2 (alone {alone['nodes'] * cell:.2f})")
        st["traces"][tid] = rec
    return kept_new


def kept_tree(sc, st, kept_new):
    """KD-tree of all nodes of earlier kept patches (+ the T-267 B1 clean rectangle on PHerc0826) plus those kept earlier in this round."""
    from scipy.spatial import cKDTree
    import tifffile
    pts = []
    if sc == "PHerc0826":
        d = C.DATA267 / "runs/B1_seed75_75_G110/out/seg"
        a = np.stack([tifffile.imread(d / f"{k}.tif").astype(np.float32) for k in "xyz"], -1)[2:217, 2:113].reshape(-1, 3)
        pts.append(a)
    for tid in st["traces"]:
        f = C.scroll_dir(sc) / "traces" / f"{tid}.kept.npy"
        if f.exists():
            pts.append(np.load(f))
    return cKDTree(np.concatenate(pts)) if pts else None


# ------------------------------------------------------------------------------------------------ render + read + pictures
def queue_patches(kept):
    qp = C.RENDER / "queue.json"
    q = json.loads(qp.read_text()) if qp.exists() else {"surfaces": []}
    have = {s["name"] for s in q["surfaces"]}
    for p in kept:
        if p["name"] in have:
            continue
        vi = C.vol_info(p["scroll"])
        q["surfaces"].append({"name": p["name"], "scroll": p["scroll"], "volume": vi["volume"], "voxel_um": vi["voxel_um"], "area_cm2": round(p["cm2"], 3), "order": 1, "onsheet_score": None,
                              "surface_text": f"T-268 trace {p['tid']} (vc_grow_seg_from_seed, public m7 surface prediction), largest clean rectangle with earlier patches excluded: {p['rows'][1]-p['rows'][0]} x {p['cols'][1]-p['cols'][0]} nodes, side {'outward' if p['side']==1 else 'inward'} (X-212)",
                              "ink_text": "none (new trace)", "source": f"T-268 X-212 {p['tid']} rows {p['rows']} cols {p['cols']}", "mesh_dir": p["mesh"], "tier": "T268"})
    C.jwrite(qp, q)


def read_patches(names):
    out = {}
    for cmd in ("produce", "consume"):
        lf = C.RENDER / "logs" / f"{cmd}-{time.strftime('%Y%m%d-%H%M%S')}.log"
        lf.parent.mkdir(parents=True, exist_ok=True)
        with open(lf, "w") as f:
            r = subprocess.run([PY, str(HERE / "read268.py"), cmd, ",".join(names)], stdout=f, stderr=subprocess.STDOUT)
        out[cmd] = r.returncode
        log(cmd, "exit", r.returncode, str(lf))
    return out


def finish_patches(kept):
    import picture268 as P
    idx_f = C.RENDER / "index.json"
    fidx = {r["surface"]: r for r in json.loads(idx_f.read_text())["surfaces"]} if idx_f.exists() else {}
    done = []
    for p in kept:
        row = fidx.get(p["name"])
        if row is None:
            p["read"] = False
            continue
        p["um_per_canvas_px"] = row["um_per_canvas_px"]
        p["canvas_px"] = row["canvas_px"]
        p["valid_mpx"] = row["valid_mpx"]
        p["nonfinite_valid_px"] = row["nonfinite_valid_px"]
        p["arc_read_minutes"] = row["arc_read_minutes"]
        p["render_identity"] = row["render_identity"]
        p["ct_GB_fetched"] = row["ct_GB_fetched"]
        p["maps_npz"] = row["maps_npz"]
        try:
            p["summary"] = P.summary(p["name"])
            pic = P.picture(p["name"], p, C.PICS / f"{p['name']}.png")
        except Exception as e:  # noqa: BLE001
            traceback.print_exc()
            p["read"] = False
            p["error"] = f"{type(e).__name__}: {e}"
            continue
        p.update({"picture": pic["png"], "picture_size_mm": pic["size_mm"], "picture_um_per_px": 40, "read": True})
        done.append(p)
        log("read", p["name"], f"{p['cm2']:.2f} cm2", "frame-correct", p["frame_correct"], "share>0.5 fc", round(p["summary"][p["frame_correct"]]["share_gt_0.5"], 3))
    return done


def contact_sheets(sc, patches):
    import picture268 as P
    rd = [p for p in patches if p["scroll"] == sc and p.get("read")]
    pages = []
    for i in range(0, len(rd), 12):
        items = []
        for p in rd[i:i + 12]:
            s = p["summary"][p["frame_correct"]]
            items.append({"name": p["name"], "patch": p, "label": f"{p['name'].split('_t268_')[1]}  {p['cm2']:.2f} cm2  {p['picture_size_mm'][0]:g}x{p['picture_size_mm'][1]:g} mm",
                          "label2": f"{'inward' if p['side']==-1 else 'outward'}, frame-correct {p['frame_correct']}; >0.5: {s['share_gt_0.5']:.3f}"})
        out = C.PICS / f"contact_{sc}_{i // 12 + 1:02d}.png"
        P.contact(items, out, f"{sc}  patches {i + 1}-{i + len(items)}")
        pages.append(str(out))
    return pages


def write_index(sc_list):
    idx = load_index()
    patches, traces = [], []
    totals = {}
    for sc in sc_list:
        st = load_state(sc)
        traces += [{k: v for k, v in t.items() if k not in ("bbox", "patch_obj")} for t in st["traces"].values()]
        patches += st["patches"]
        pk = [p for p in st["patches"] if p.get("read")]
        totals[sc] = {"traces": sum(1 for t in st["traces"].values() if t.get("status") not in ("skipped_on_traced_sheet",)), "skipped_on_traced_sheet": sum(1 for t in st["traces"].values() if t.get("status") == "skipped_on_traced_sheet"),
                      "patches_kept": len(st["patches"]), "patches_read": len(pk), "clean_cm2_traced_new": round(sum(p["cm2"] for p in st["patches"]), 3), "clean_cm2_read": round(sum(p["cm2"] for p in pk), 3),
                      "patches_ge_4cm2": sum(1 for p in st["patches"] if p["cm2"] >= 4.0), "rounds": len(st["rounds"])}
    old = {p["name"] for p in patches}
    for p in idx["patches"]:                     # keep other scrolls' entries not handled in this call
        if p["scroll"] not in sc_list and p["name"] not in old:
            patches.append(p)
    idx["patches"], idx["traces"] = patches, traces
    idx["totals"] = {**idx.get("totals", {}), **totals}
    idx["updated"] = time.strftime("%Y-%m-%d %H:%M:%S")
    C.jwrite(INDEX, idx)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scroll", required=True)
    ap.add_argument("--max-traces", type=int, default=32)
    ap.add_argument("--rounds", type=int, default=1000)
    ap.add_argument("--parallel", type=int, default=C.PARALLEL)
    a = ap.parse_args()
    sc = a.scroll
    cands = json.loads((C.scroll_dir(sc) / "seed-candidates.json").read_text())["seeds"]
    (C.scroll_dir(sc) / "traces").mkdir(parents=True, exist_ok=True)
    st = load_state(sc)
    fails = 0
    for _ in range(a.rounds):
        if STOP.exists():
            log("STOP file present; stopping before the next round")
            break
        ntr = sum(1 for t in st["traces"].values() if t.get("status") != "skipped_on_traced_sheet")
        if ntr >= a.max_traces:
            log("max traces reached", ntr)
            break
        chosen = pick_seeds(sc, st, cands, min(a.parallel, a.max_traces - ntr))
        save_state(sc, st)
        if not chosen:
            log("no seeds left")
            break
        rnd = len(st["rounds"]) + 1
        t0 = time.time()
        try:
            rdir, secs, insp = run_round_container(sc, rnd, chosen)
            kept = audit_round(sc, rnd, rdir, chosen, st)
            rinfo = {"round": rnd, "seeds": [c["no"] for c in chosen], "container_seconds": round(secs), "container_exit": insp, "kept": [p["name"] for p in kept]}
            st["rounds"].append(rinfo)
            save_state(sc, st)
            if kept:
                st["patches"] += kept
                save_state(sc, st)
                queue_patches(kept)
                rr = read_patches([p["name"] for p in kept])
                rinfo["read_exit"] = rr
                done = finish_patches(kept)
                save_state(sc, st)
                contact_sheets(sc, st["patches"])
            rinfo["round_seconds"] = round(time.time() - t0)
            save_state(sc, st)
            write_index([sc])
            fails = 0
        except Exception as e:  # noqa: BLE001
            traceback.print_exc()
            fails += 1
            st["rounds"].append({"round": rnd, "seeds": [c["no"] for c in chosen], "error": f"{type(e).__name__}: {e}"})
            for c in chosen:
                st["traces"].setdefault(f"{sc}_s{c['no']:03d}", {"status": "round_error", "seed": c["seed"], "no": c["no"]})
            save_state(sc, st)
            if fails >= 3:
                log("three failed rounds in a row: stopping")
                break
    write_index([sc])
    log("done", json.dumps(load_index().get("totals", {}).get(sc)))


if __name__ == "__main__":
    main()
