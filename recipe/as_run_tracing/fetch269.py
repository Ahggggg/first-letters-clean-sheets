"""T-269 (X-213) data for the new sites: listing, metadata, (control only) level-3 pyramid and densest-cube block, plan, guarded fetch, Docker volume.  The T-268 / T-267 fetch code is reused UNCHANGED; only
(a) the S3 prefix is taken from the site's SCROLL while local folders are per SITE (fetch267's base() wrapped), and (b) every guarded pack is taken under the experiment's guard lock (common269.guard_lock) so that a
block fetch and the CT fetch of a render never run the D-005 guard at the same time.  Every byte still goes through work/download_guard_d005.py via fetch267.guarded_pack.
  python fetch269.py list SITE | meta SITE | level3 SITE | block SITE | plan SITE | fetch SITE | volume SITE
 'block' is for the control site only (X-212's rule); the anchored sites' boxes come from sites269.py."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common269 as C  # noqa: E402
import fetch268 as F8  # noqa: E402

F = F8.F
_base = F.base
_gp = F.guarded_pack
F.base = lambda site: _base(C.scroll_of(site))                       # S3 prefix of the scroll; dest_of() keeps the site as the local folder name


def locked_pack(name, objs, tag):
    with C.guard_lock(f"experimenter-39 fetch269 {name}"):
        return _gp(name, objs, tag)


F.guarded_pack = locked_pack


def cmd_block(sc):
    F8.cmd_block(sc)                                                  # writes data/t269-trace/<site>/box.json by the X-212 rule
    box = json.loads((C.DATA / sc / "box.json").read_text())["box_xyz"]
    s = json.loads(C.SITES_JSON.read_text())
    s["sites"][sc]["box_xyz"] = box
    s["sites"][sc]["box_sides"] = [box[1] - box[0], box[3] - box[2], box[5] - box[4]]
    C.jwrite(C.SITES_JSON, s)


def cmd_plan(sc):
    assert C.box_of(sc), "box not set"
    F8.cmd_plan(sc)


if __name__ == "__main__":
    cmd, sc = sys.argv[1], sys.argv[2]
    {"list": F8.cmd_list, "meta": F8.cmd_meta, "level3": F8.cmd_level3, "block": cmd_block, "plan": cmd_plan, "fetch": F8.cmd_fetch, "volume": F8.cmd_volume}[cmd](sc)
