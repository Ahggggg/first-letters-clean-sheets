"""T-269 (X-213) shared constants and site registry.  CPU only.  The T-268 code (seeds268.py, audit268.py, picture268.py, fetch268.py helpers, summarize268.py) is reused UNCHANGED:
this module is installed under the name `common268` (sys.modules alias) so that those files, which do `import common268 as C`, see the T-269 site registry and the data/t269-trace/ folders.
It first imports the REAL common268 and re-exports every constant (CH, CELL_CH, PARALLEL, GENERATIONS, CACHE_SIZE, MIN_KEEP_CM2, EXCL_VOX, SKIP_TRACED_VOX, ROUND_SEP_VOX, WALL_MARGIN_VOX, IMAGE, DOCKER_BIN) unchanged.
SITE = one mirror box.  A site names a scroll, a level-0 mirror (pred.zarr + ngrid), a Docker volume and a box [xlo, xhi, ylo, yhi, zlo, zhi] (level-0 voxels).  Sites: data/t269-trace/sites.json (written by sites269.py).
Conventions as in common268: tifxyz arrays (row, col, [x, y, z]); node spacing 20 voxels at scale 0.05; zarr array order (z, y, x).
"""
import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "t268-trace"))
sys.path.insert(0, str(HERE.parent / "t267-trace"))
import common268 as _C268  # noqa: E402  (the real module; the alias is installed at the bottom)

for _k, _v in vars(_C268).items():
    if _k.isupper() and not _k.startswith("_"):
        globals()[_k] = _v

PROB = _C268.PROB
FAC = _C268.FAC
DATA267 = _C268.DATA267
DATA268 = PROB / "data/t268-trace"
DATA = PROB / "data/t269-trace"
RENDER = DATA / "render"
PICS = PROB / "outputs/t269-trace"
SITES_JSON = DATA / "sites.json"
GUARD_LOCK = DATA / "guard.lock"
queue_volumes = _C268.queue_volumes
jwrite = _C268.jwrite


def sites():
    return json.loads(SITES_JSON.read_text())["sites"]


def site(s):
    return sites()[s]


def scroll_of(s):
    return site(s)["scroll"]


def vol_info(s):
    sc = scroll_of(s)
    v = queue_volumes()[sc]
    return {"scroll": sc, "site": s, "volume": v["volume"], "voxel_um": v["voxel_um"], "shape": tuple(v["shape"])}      # shape = (slices z, height y, width x)


def scroll_dir(s):
    return DATA / s


def mirror_dir(s):
    return Path(site(s)["mirror"])


def box_of(s):
    return site(s)["box_xyz"]


def docker_volume(s):
    return site(s)["docker_volume"]


def prior_dirs(s):
    """folders of EARLIER (other-experiment) state of the same site, read-only: their traced nodes / kept patches must not be traced or kept again."""
    return [Path(p) for p in site(s).get("prior_state_dirs", [])]


def related_dirs(s):
    """every folder holding traces (state.json + traces/*.npy) of the same scroll that count as already-traced ground for site s: its own, its priors, and those of the other sites of the scroll."""
    out = []
    for name, d in sites().items():
        if d["scroll"] != scroll_of(s):
            continue
        out.append(scroll_dir(name))
        out += [Path(p) for p in d.get("prior_state_dirs", [])]
    return out


class guard_lock:
    """Serialises every use of the D-005 download guard by this experiment (a fetch and the CT fetch inside a render run): exclusive create of data/t269-trace/guard.lock; a lock older than 8 h is stale."""
    def __init__(self, who):
        self.who = who

    def __enter__(self):
        GUARD_LOCK.parent.mkdir(parents=True, exist_ok=True)
        while True:
            try:
                fd = os.open(str(GUARD_LOCK), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                os.write(fd, f"{self.who} {time.strftime('%Y-%m-%d %H:%M:%S')} pid {os.getpid()}\n".encode())
                os.close(fd)
                return self
            except FileExistsError:
                try:
                    if time.time() - GUARD_LOCK.stat().st_mtime > 8 * 3600:
                        GUARD_LOCK.unlink()
                        continue
                except FileNotFoundError:
                    continue
                time.sleep(10)

    def __exit__(self, *a):
        try:
            GUARD_LOCK.unlink()
        except FileNotFoundError:
            pass


sys.modules["common268"] = sys.modules[__name__]
