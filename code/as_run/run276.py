"""run276.py - X-218 driver (T-276).  Run only through  python tools/run.py --problem herculaneum-ink --experiment X-218 -- <python> run276.py all
Stages (all as registered in X-218; xsec276.py holds the measure):
  identity   my sampler vs the stored 28-plane stack of PHerc0813_t269_s000 (256 x 256 block at canvas (384, 1000)); if it fails nothing else runs
  patches    the 10 straightened sections per kept patch (X-211, X-212, X-213 meshes), measure -> data/t276-xsec/patches/<name>.json, sections -> data/t276-xsec/sections/<name>.npz
  g1         synthetic-shear gate on 192-px windows of the informative sections of 20 patches drawn with numpy seed 276
  g2         public control surfaces (flat / swirled) with the same lines
  summary    all numbers -> data/t276-xsec/summary.json (printed)
CPU only, one process, no Arc, no Docker, no download.  numpy seed 276."""
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
PROB = HERE.parents[1]
sys.path.insert(0, str(HERE))
import xsec276 as X  # noqa: E402

OUT = PROB / "data/t276-xsec"
FAC = PROB / "data/picture-factory"
DEPTHS = np.arange(-32, 33)
FRACS = (0.1, 0.3, 0.5, 0.7, 0.9)
CT_ROOTS_EXTRA = [Path("<archived CT cache>/picture-factory/ctcache"), PROB / "data/t237-new-candidates/ctcache"]
QUEUE_VOLS = json.loads((FAC / "queue.json").read_text())["volumes"]


def vol_for(scroll):
    v = QUEUE_VOLS[scroll]
    tag = f"{scroll}-{v['volume']}"
    roots = [FAC / "ctcache" / tag] + [r / tag for r in CT_ROOTS_EXTRA]
    return X.ChunkVolume([r for r in roots if r.exists()], v["shape"]), tuple(v["shape"])


def patch_list():
    out = []
    for tn, idx in (("t267", "render/index.json"), ("t268", "index.json"), ("t269", "index.json")):
        d = json.loads((PROB / f"data/{tn}-trace/{idx}").read_text())
        for p in (d["surfaces"] if "surfaces" in d else d["patches"]):
            name = p.get("name") or p["surface"]
            mesh = PROB / f"data/{tn}-trace/render/meshes/{name}.tifxyz"
            out.append({"name": name, "scroll": p["scroll"], "run": tn, "mesh": mesh, "stack": PROB / f"data/{tn}-trace/render/stacks/{name}.npy",
                        "side": p.get("side"), "stack_valid": PROB / f"data/{tn}-trace/render/stacks/{name}.valid.npy"})
    return out


def identity():
    p = [q for q in patch_list() if q["name"] == "PHerc0813_t269_s000"][0]
    geom = X.load_tifxyz(p["mesh"])
    vol, shape = vol_for(p["scroll"])
    H, W = X.canvas_size(geom)
    stack = np.load(p["stack"], mmap_mode="r")
    valid = np.load(p["stack_valid"], mmap_mode="r")
    r0, c0 = 384, 1000
    rows, cols = np.arange(r0, r0 + 256), np.arange(c0, c0 + 256)
    P, N = X.surface(geom, rows, cols)
    good = np.isfinite(P).all(-1) & np.isfinite(N).all(-1)
    offs = -1.0 * (np.arange(28) - 13.5)
    base, direction = P[good][:, ::-1], N[good][:, ::-1]
    tot, ok1, maxd, unknown = 0, 0, 0, 0
    for k in range(28):
        samples = base + offs[k] * direction
        v, bad = vol.sample(samples)
        mine = np.rint(np.clip(v, 0, 255)).astype(np.int16)
        ref = np.asarray(stack[k, r0:r0 + 256, c0:c0 + 256])[good].astype(np.int16)
        m = ~bad
        unknown += int(bad.sum())
        d = np.abs(mine[m] - ref[m])
        tot += int(m.sum())
        ok1 += int((d <= 1).sum())
        maxd = max(maxd, int(d.max()))
    res = {"stack": p["name"], "block": [r0, c0, 256, 256], "canvas": [H, W], "valid_px_equal": bool(np.array_equal(good, np.asarray(valid[r0:r0 + 256, c0:c0 + 256]))),
           "compared": tot, "share_within_1": ok1 / tot, "max_abs_diff": maxd, "unknown_samples": unknown, "pass": bool(ok1 / tot >= 0.999 and unknown == 0)}
    (OUT).mkdir(parents=True, exist_ok=True)
    (OUT / "identity.json").write_text(json.dumps(res, indent=1))
    print("IDENTITY", json.dumps(res))
    return res["pass"]


def run_sections(geom, vol, side, rng, cap=None, name=""):
    H, W = X.canvas_size(geom)
    secs, arrays = [], {}
    for kind in ("row", "col"):
        n_idx, n_s = (H, W) if kind == "row" else (W, H)
        for f in FRACS:
            idx = int(round(f * (n_idx - 1)))
            span = (0, n_s)
            if cap:       # public controls: longest geometrically valid run, capped, centred
                s_idx = np.arange(0, n_s)
                P, N = X.surface(geom, [idx], s_idx) if kind == "row" else X.surface(geom, s_idx, [idx])
                P = P[0] if kind == "row" else P[:, 0]
                g = np.isfinite(P).all(-1)
                best, cur, start = (0, 0), 0, 0
                for i, v in enumerate(np.r_[g, False]):
                    if v and cur == 0:
                        start = i
                    if v:
                        cur += 1
                    elif cur:
                        if cur > best[1] - best[0]:
                            best = (start, start + cur)
                        cur = 0
                a, b = best
                if b - a > cap:
                    mid = (a + b) // 2
                    a, b = mid - cap // 2, mid - cap // 2 + cap
                span = (int(a), int(b))
                if b - a < 200:
                    secs.append({"kind": kind, "index": idx, "skipped": "no valid run >= 200 px"})
                    continue
            S, valid, L, s_idx = X.section(vol, geom, kind, idx, span, DEPTHS, side)
            m = X.measure(S, valid, DEPTHS, rng)
            rec = {"kind": kind, "index": idx, "span": list(span), "frac": f, "arclen_vox": L, **m}
            if np.isfinite(m["tilt"]):
                rec["drift_vox"] = float(L * np.tan(np.radians(m["tilt"])))
                rec["local"] = [[int(a), float(t), float(c)] for a, t, c in X.local_tilts(S, valid, DEPTHS)]
            rec["unknown_cols"] = int((~valid).sum())
            secs.append(rec)
            arrays[f"{kind}_{idx}"] = S
            arrays[f"{kind}_{idx}_valid"] = valid
    return secs, arrays


def stage_patches():
    rng = np.random.default_rng(276)
    (OUT / "patches").mkdir(parents=True, exist_ok=True)
    (OUT / "sections").mkdir(parents=True, exist_ok=True)
    for i, p in enumerate(patch_list()):
        f = OUT / "patches" / f"{p['name']}.json"
        if f.exists():
            continue
        t0 = time.time()
        geom = X.load_tifxyz(p["mesh"])
        vol, shape = vol_for(p["scroll"])
        side, minority = X.frame_side(geom, shape)
        secs, arrays = run_sections(geom, vol, side, rng)
        np.savez_compressed(OUT / "sections" / f"{p['name']}.npz", **{k: v.astype(np.float16) if v.dtype != bool else v for k, v in arrays.items()})
        rec = {"name": p["name"], "scroll": p["scroll"], "run": p["run"], "side_measured": side, "side_index": p["side"], "frame_minority_rect": minority, "sections": secs, "seconds": round(time.time() - t0, 1)}
        f.write_text(json.dumps(rec, indent=1))
        inf = [s for s in secs if s.get("informative")]
        print(i, p["name"], "informative", len(inf), "/", len(secs), "med|tilt|", round(float(np.median([abs(s["tilt"]) for s in inf])), 2) if inf else None, "unknown cols", sum(s.get("unknown_cols", 0) for s in secs), rec["seconds"], "s", flush=True)


def stage_g1():
    rng = np.random.default_rng(276)
    names = sorted(p.stem for p in (OUT / "patches").glob("*.json"))
    pick = list(rng.choice(len(names), 20, replace=False))
    res = []
    for ix in pick:
        name = names[ix]
        rec = json.loads((OUT / "patches" / f"{name}.json").read_text())
        z = np.load(OUT / "sections" / f"{name}.npz")
        for s in rec["sections"]:
            if not s.get("informative"):
                continue
            key = f"{s['kind']}_{s['index']}"
            S = z[key].astype(np.float32)
            valid = z[key + "_valid"]
            for a0 in range(0, S.shape[1] - 192 + 1, 192):
                Sw, vw = S[:, a0:a0 + 192], valid[a0:a0 + 192]
                if vw.sum() < 150:
                    continue
                m0 = X.measure(Sw, vw, DEPTHS, rng)
                if not m0["informative"]:
                    continue
                for a in (-4, -2, 2, 4):
                    m1 = X.measure(X.shear(Sw, DEPTHS, a), vw, DEPTHS, rng)
                    if m1["informative"]:
                        res.append({"name": name, "kind": s["kind"], "idx": s["index"], "a0": a0, "a": a, "t0": m0["tilt"], "t1": m1["tilt"], "err": m1["tilt"] - m0["tilt"] - a})
    err = np.array([r["err"] for r in res])
    # wrap errors to (-90, 90]
    err = (err + 90) % 180 - 90
    out = {"patches_drawn": [names[i] for i in pick], "n_windows_sheared": len(res), "median_err": float(np.median(err)) if len(err) else None,
           "share_within_1": float((np.abs(err) <= 1).mean()) if len(err) else None, "share_within_2": float((np.abs(err) <= 2).mean()) if len(err) else None,
           "by_a": {str(a): {"n": int(sum(1 for r in res if r["a"] == a)), "median_err": float(np.median([(r['err'] + 90) % 180 - 90 for r in res if r['a'] == a]))} for a in (-4, -2, 2, 4) if any(r["a"] == a for r in res)}}
    out["pass"] = bool(len(err) and abs(out["median_err"]) <= 1.0 and out["share_within_2"] >= 0.70)
    (OUT / "g1.json").write_text(json.dumps(out, indent=1))
    print("G1", json.dumps(out))


CONTROLS = {"flat": ["PHerc0826_auto_grown_20260902140841364", "PHerc0826_auto_grown_20260902140900360", "PHerc0826_auto_grown_20260902140849408"],
            "swirled": ["PHerc0826_v3_s0", "PHerc0358_v3_s7", "PHerc0257_v2_s0"]}


def stage_g2():
    rng = np.random.default_rng(276)
    (OUT / "controls").mkdir(parents=True, exist_ok=True)
    (OUT / "control_sections").mkdir(parents=True, exist_ok=True)
    for grp, names in CONTROLS.items():
        for nm in names:
            f = OUT / "controls" / f"{nm}.json"
            if f.exists():
                continue
            scroll = nm.split("_")[0]
            geom = X.load_tifxyz(FAC / "meshes" / f"{nm}.tifxyz")
            vol, shape = vol_for(scroll)
            side, minority = X.frame_side(geom, shape)
            secs, arrays = run_sections(geom, vol, side, rng, cap=2500)
            np.savez_compressed(OUT / "control_sections" / f"{nm}.npz", **{k: v.astype(np.float16) if v.dtype != bool else v for k, v in arrays.items()})
            f.write_text(json.dumps({"name": nm, "group": grp, "scroll": scroll, "side_measured": side, "frame_minority": minority, "sections": secs}, indent=1))
            inf = [s for s in secs if s.get("informative")]
            print("G2", grp, nm, "informative", len(inf), "/", len(secs), flush=True)


def stage_summary():
    P = [json.loads(f.read_text()) for f in sorted((OUT / "patches").glob("*.json"))]
    C = [json.loads(f.read_text()) for f in sorted((OUT / "controls").glob("*.json"))]
    g1 = json.loads((OUT / "g1.json").read_text()) if (OUT / "g1.json").exists() else None

    def inf(secs):
        return [s for s in secs if s.get("informative")]
    allp = [s for r in P for s in r["sections"] if "tilt" in s]
    infp = [s for s in allp if s["informative"]]
    per_patch = []
    for r in P:
        i = inf(r["sections"])
        per_patch.append({"name": r["name"], "scroll": r["scroll"], "run": r["run"], "n_sections": len(r["sections"]), "n_informative": len(i),
                          "median_abs_tilt": float(np.median([abs(s["tilt"]) for s in i])) if i else None,
                          "median_abs_drift": float(np.median([abs(s["drift_vox"]) for s in i])) if i else None,
                          "share_drift_ge_20": float(np.mean([abs(s["drift_vox"]) >= 20 for s in i])) if i else None})
    ctl = {}
    for grp in ("flat", "swirled"):
        s = [x for r in C if r["group"] == grp for x in inf(r["sections"])]
        allc = [x for r in C if r["group"] == grp for x in r["sections"] if "tilt" in x]
        ctl[grp] = {"n_sections": len(allc), "n_informative": len(s), "median_abs_tilt": float(np.median([abs(x["tilt"]) for x in s])) if s else None,
                    "median_abs_drift": float(np.median([abs(x["drift_vox"]) for x in s])) if s else None,
                    "share_drift_ge_20": float(np.mean([abs(x["drift_vox"]) >= 20 for x in s])) if s else None,
                    "per_surface": {r["name"]: {"n_informative": len(inf(r["sections"])), "median_abs_tilt": float(np.median([abs(x["tilt"]) for x in inf(r["sections"])])) if inf(r["sections"]) else None} for r in C if r["group"] == grp}}
    g2 = None
    if ctl["flat"]["median_abs_tilt"] and ctl["swirled"]["median_abs_tilt"] is not None:
        g2 = {"ratio_swirled_over_flat": ctl["swirled"]["median_abs_tilt"] / ctl["flat"]["median_abs_tilt"], "pass": bool(ctl["swirled"]["median_abs_tilt"] >= 1.5 * ctl["flat"]["median_abs_tilt"])}
    by_scroll = {}
    for pp in per_patch:
        by_scroll.setdefault(pp["scroll"], []).append(pp)
    sc = {k: {"patches": len(v), "median_of_patch_median_abs_tilt": float(np.median([x["median_abs_tilt"] for x in v if x["median_abs_tilt"] is not None]))} for k, v in by_scroll.items()}
    out = {"n_patches": len(P), "n_sections": len(allp), "n_informative": len(infp), "share_informative": len(infp) / max(1, len(allp)),
           "patches_with_ge3_informative": sum(1 for pp in per_patch if pp["n_informative"] >= 3),
           "uncached_cols_share": float(sum(s.get("unknown_cols", 0) for r in P for s in r["sections"]) / max(1, sum(s.get("n_valid_cols", 0) + s.get("unknown_cols", 0) for r in P for s in r["sections"]))),
           "patch_median_abs_tilt_all_informative": float(np.median([abs(s["tilt"]) for s in infp])) if infp else None,
           "patch_median_abs_drift_all_informative": float(np.median([abs(s["drift_vox"]) for s in infp])) if infp else None,
           "share_informative_sections_drift_ge_20": float(np.mean([abs(s["drift_vox"]) >= 20 for s in infp])) if infp else None,
           "share_patches_median_drift_lt_20": float(np.mean([pp["median_abs_drift"] < 20 for pp in per_patch if pp["median_abs_drift"] is not None])),
           "controls": ctl, "G1": g1, "G2": g2, "by_scroll": sc, "per_patch": per_patch}
    off = [pp for pp in per_patch if pp["name"] in ("PHerc0826_t268_s004", "PHerc0826_t268_s014", "PHerc0826_t268_s033")]
    rest = [pp for pp in per_patch if pp["scroll"] == "PHerc0826" and pp not in off and pp["median_abs_tilt"] is not None]
    out["off_prediction_t268"] = {"off": [(x["name"], x["median_abs_tilt"]) for x in off], "median_off": float(np.median([x["median_abs_tilt"] for x in off if x["median_abs_tilt"] is not None])) if off else None,
                                  "median_other_0826": float(np.median([x["median_abs_tilt"] for x in rest])) if rest else None}
    (OUT / "summary.json").write_text(json.dumps(out, indent=1))
    print("SUMMARY", json.dumps({k: v for k, v in out.items() if k != "per_patch"}, indent=1))


if __name__ == "__main__":
    st = sys.argv[1]
    OUT.mkdir(parents=True, exist_ok=True)
    if st in ("all", "identity"):
        if not identity():
            print("IDENTITY FAILED: nothing else is run")
            sys.exit(2)
    if st in ("all", "patches"):
        stage_patches()
    if st in ("all", "g1"):
        stage_g1()
    if st in ("all", "g2"):
        stage_g2()
    if st in ("all", "summary"):
        stage_summary()
