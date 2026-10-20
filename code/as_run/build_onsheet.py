"""build_onsheet.py (T-276): onsheet/ folder of the staging repository: the X-202 on-sheet score code as run, thresholds, the 120-surface ranking table (+ source group), group statistics, example pictures.  CPU only."""
import csv
import json
import shutil
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROB = HERE.parents[1]
STAGE = Path("<path of this repository>")
FAC = PROB / "data/picture-factory"
D = FAC / "onsheet_d"
out = STAGE / "onsheet"
(out / "code_as_run").mkdir(parents=True, exist_ok=True)
(out / "examples").mkdir(exist_ok=True)
for src in (PROB / "work/t257-onsheet/onsheet257.py", PROB / "work/t257-onsheet/onsheet257c.py", PROB / "work/t259-onsheet-d/onsheet259.py", PROB / "work/t245-picture-factory/pic245.py"):
    shutil.copyfile(src, out / "code_as_run" / src.name)
shutil.copyfile(D / "thresholds259.json", out / "thresholds259.json")
shutil.copyfile(D / "thresholds259.sha256", out / "thresholds259.sha256")
rk = json.loads((D / "ranking.json").read_text())
q = {s["name"]: s for s in json.loads((FAC / "queue.json").read_text())["surfaces"]}


def group(src):
    if "flummoxjr" in src:
        return "flummoxjr/gp13-ink-detectability"
    if "pscamillo" in src:
        return "pscamillo/vesuvius-eligible-meshes"
    if "kadenpool" in src:
        return "kadenpool/scroll-lineup"
    if "Bullo27" in src:
        return "Bullo27/first-letters-survey"
    return "other"


rows = []
for r in sorted(rk["rows"], key=lambda r: r["rank_area"]):
    s = q[r["surface"]]
    rows.append({"rank_by_on_sheet_area": r["rank_area"], "surface": r["surface"], "scroll": r["scroll"], "source_group": group(s["source"]), "auto_grown_name": int("auto_grown" in r["surface"]),
                 "on_sheet_area_cm2": round(r["area_cm2"], 2), "valid_area_cm2": round(r["valid_area_cm2"], 2), "on_sheet_share": round(r["share"], 3), "tiles_valid": r["tiles_valid"], "tiles_on_sheet": r["tiles_on_sheet"],
                 "rank_by_share": r["rank_share"], "depth_contrast_C_median_pct": round(r["C_median"], 1), "texture_Tx_median": round(r["Tx_median"], 4), "um_per_px": round(r["um_per_px"], 3)})
with open(out / "ranking_120_surfaces.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
sh = np.array([r["on_sheet_share"] for r in rows])
stats = {"n": len(rows), "median_share": float(np.median(sh)), "share_quartiles_p10_p25_p50_p75_p90": np.percentile(sh, [10, 25, 50, 75, 90]).round(3).tolist(), "n_share_ge_0.5": int((sh >= 0.5).sum()),
         "on_sheet_cm2": round(sum(r["on_sheet_area_cm2"] for r in rows), 1), "valid_cm2": round(sum(r["valid_area_cm2"] for r in rows), 1), "by_source": {}, "by_name_kind": {}}
for g in sorted({r["source_group"] for r in rows}):
    x = np.array([r["on_sheet_share"] for r in rows if r["source_group"] == g])
    stats["by_source"][g] = {"n": int(len(x)), "median_share": round(float(np.median(x)), 3), "n_share_ge_0.5": int((x >= 0.5).sum())}
for k, lab in ((1, "name contains auto_grown"), (0, "other names")):
    x = np.array([r["on_sheet_share"] for r in rows if r["auto_grown_name"] == k])
    stats["by_name_kind"][lab] = {"n": int(len(x)), "median_share": round(float(np.median(x)), 3), "n_share_ge_0.5": int((x >= 0.5).sum())}
stats["validation_from_X-202"] = {"V1_heldout_known_writing": "0.74 / 0.74 / 0.96", "nulls": "0.000 (depth shuffle) / 0.005 (spatial shuffle)", "V2": rk["V2"]}
(out / "ranking_stats.json").write_text(json.dumps(stats, indent=1))
print(json.dumps(stats, indent=1))
ex = [(D / "top20/02_PHerc0826_auto_grown_20260902140841364.png", "ex1_rank2_regular_fibre_grid_95pct_on_sheet.png"),
      (D / "top20/03_PHerc0826_v3_s0.png", "ex2_rank3_scores_high_but_swirled_false_positive.png"),
      (FAC / "framed/PHerc0826_auto_grown_20260902140900360.png", "ex3_flat_reference_PHerc0826_auto_grown_900360.png"),
      (FAC / "framed/PHerc0800_z15264_w100.png", "ex4_flat_reference_PHerc0800_w100.png"),
      (FAC / "framed/PHerc0358_v3_s7.png", "ex5_swirled_PHerc0358_v3_s7.png"),
      (FAC / "framed/PHerc0257_v2_s0.png", "ex6_swirled_PHerc0257_v2_s0.png")]
for s, n in ex:
    shutil.copyfile(s, out / "examples" / n)
shutil.copyfile(D / "onsheet_map_all.png", out / "examples/onsheet_map_all_120.png")
