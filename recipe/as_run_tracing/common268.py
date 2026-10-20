"""T-268 shared constants and small helpers (X-212).  CPU only.  Everything registered in X-212 lives in these files:
  seeds268.py   seed cells (grid of 384-voxel cells -> one planar-band seed per cell) from a local level-0 mirror of the public m7 surface prediction
  audit268.py   per-trace audit (frame side, folds, stretch, mirror-box wall) and the largest clean rectangle with earlier patches excluded
  drive268.py   the orchestrator: rounds of 4 tracer runs in ONE Docker container, audit, render, Arc read, pictures, index
  read268.py    thin wrapper of the UNCHANGED T-245 picture factory (like work/t267-trace/read267.py) for a queue of many patches
  picture268.py 40 um/px picture with the frame-correct order marked, contact sheets, reading summaries
Conventions: tifxyz arrays (row, col, [x, y, z]); node spacing 20 voxels at scale 0.05; level-0 voxel coordinates (x, y, z); zarr array order (z, y, x).
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROB = HERE.parents[1]
DATA = PROB / "data/t268-trace"
DATA267 = PROB / "data/t267-trace"
RENDER = DATA / "render"
PICS = PROB / "outputs/t268-trace"
FAC = PROB / "data/picture-factory"
sys.path.insert(0, str(PROB / "work/t267-trace"))

IMAGE = "ghcr.io/scrollprize/villa/volume-cartographer:main"
DOCKER_BIN = "<local path>"
CH = 192
CELL_CH = 2                  # a seed cell = 2 x 2 x 2 chunks = 384 voxels
PARALLEL = 4                 # tracer processes in the one container
GENERATIONS = 75             # the official recipe (X-210: 152 x 152 nodes)
CACHE_SIZE = 3_000_000_000
MIN_KEEP_CM2 = 1.0           # a patch is kept (rendered and read) iff its largest clean rectangle, with earlier patches excluded, is at least this
EXCL_VOX = 12.0              # a node within this many voxels of an earlier kept patch's node counts as covered
SKIP_TRACED_VOX = 16.0       # a seed within this many voxels of a node of an earlier whole trace is skipped (same sheet already grown)
ROUND_SEP_VOX = 1000.0       # seeds of one round are at least this far apart
WALL_MARGIN_VOX = 700.0      # seeds are at least this far from the mirror-box faces


def queue_volumes():
    return json.loads((FAC / "queue.json").read_text())["volumes"]


def vol_info(sc):
    v = queue_volumes()[sc]
    return {"scroll": sc, "volume": v["volume"], "voxel_um": v["voxel_um"], "shape": tuple(v["shape"])}      # shape = (slices z, height y, width x)


def scroll_dir(sc):
    return DATA / sc


def mirror_dir(sc):
    """directory holding pred.zarr/ and ngrid/ for the scroll (PHerc0826 reuses the T-267 mirror)."""
    return DATA267 / sc if sc == "PHerc0826" else DATA / sc


def box_of(sc):
    if sc == "PHerc0826":
        return json.loads((DATA267 / "plan-B_0826wide.json").read_text())["box_xyz"]
    return json.loads((DATA / sc / "box.json").read_text())["box_xyz"]


def docker_volume(sc):
    return "t267-0826" if sc == "PHerc0826" else f"t268-{sc.lower()}"


def jwrite(path, obj):
    import os
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(path) + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1, default=str))
    os.replace(tmp, path)
