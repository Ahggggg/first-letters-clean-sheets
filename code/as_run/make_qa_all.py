"""make_qa_all.py (T-276): one QA picture per patch into the staging repository qa/<scroll>/<name>_qa.png (+ the project's original both-orders picture as <name>_both_orders.png).
Reads data/t276-xsec/records.json (build_data.py), the X-218 sections (data/t276-xsec/sections, patches/*.json), the stored 28-plane stack (plane 14) and the stored Arc maps.
CPU only.  python make_qa_all.py [--only NAME ...] [--stage DIR]"""
import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROB = HERE.parents[1]
STAGE = Path("<path of this repository>")
sys.path.insert(0, str(STAGE / "code"))
import make_qa as Q  # noqa: E402

OUT = PROB / "data/t276-xsec"
DEPTHS = np.arange(-32, 33)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--stage", default=str(STAGE))
    a = ap.parse_args()
    stage = Path(a.stage)
    recs = json.loads((OUT / "records.json").read_text())
    for r in recs:
        if a.only and r["name"] not in a.only:
            continue
        pj = json.loads((OUT / "patches" / f"{r['name']}.json").read_text())
        z = np.load(OUT / "sections" / f"{r['name']}.npz")
        mid = {}
        for s in pj["sections"]:
            if s.get("frac") == 0.5:
                mid[s["kind"]] = s
        stack = np.load(r["stack"], mmap_mode="r")
        valid = np.load(Path(str(r["stack"]).replace(".npy", ".valid.npy")), mmap_mode="r")
        scan = np.asarray(stack[14]).astype(np.float32)
        valid = np.asarray(valid).astype(bool)
        zm = np.load(r["maps"])
        mv = zm["valid"].astype(bool)
        F = np.minimum(zm["rv2__F"].astype(np.float32), zm["d9v4C__F"].astype(np.float32))
        R = np.minimum(zm["rv2__R"].astype(np.float32), zm["d9v4C__R"].astype(np.float32))
        ink = Q.ink_panel(F, R, mv, r["frame_correct"], r["um_per_canvas_px"])
        secs = []
        for kind in ("row", "col"):
            s = mid[kind]
            key = f"{kind}_{s['index']}"
            S = z[key].astype(np.float32)
            v = z[key + "_valid"]
            if s.get("informative"):
                txt = f"tilt {s['tilt']:+.2f} deg, coherence {s['coherence']:.2f}, implied drift {s['drift_vox']:+.0f} voxels over {s['arclen_vox']:.0f} voxels"
            else:
                txt = "not informative by the noise-null rule" if "tilt" in s else "no valid samples"
            secs.append({"S": S, "valid": v, "label_text": txt})
        lines = [("row", mid["row"]["index"]), ("col", mid["col"]["index"])]
        ra = r["rect_audit"]
        side_txt = "outward" if r["side"] == 1 else "inward"
        title = f"{r['name']}   {r['scroll']}   {r['cm2']:.2f} cm2   ({ra['grid'][0]} x {ra['grid'][1]} nodes)"
        sub = (f"scan {r['voxel_um']} um, volume {r['volume_id']}   seed (x,y,z) {tuple(int(v) for v in r['seed'])}   {r['experiment']}, trace {r['tid']}   "
               f"frame side {side_txt} -> frame-correct order {r['frame_correct']}")
        foot = (f"mesh checks on this crop: frame minority {ra['frame']['minority']:.4f}, fold share (interior) {ra['fold_share_interior']:.4f}, stretch share {ra['stretch_share']:.4f}.  "
                f"Cross-sections and maps: see README (limits: small patch, mesh-level checks plus two sections, no letters found).")
        vox = r["um_per_canvas_px"]
        out = stage / "qa" / r["scroll"]
        out.mkdir(parents=True, exist_ok=True)
        size = Q.render_qa(out / f"{r['name']}_qa.png", title, sub, scan, valid, lines, secs, ink, DEPTHS, vox, foot)
        shutil.copyfile(r["picture"], out / f"{r['name']}_both_orders.png")
        print(r["name"], size, flush=True)


if __name__ == "__main__":
    main()
