"""T-268 (X-212): render traced patches with the project renderer and read them with Reader v2 and d9v4C in both orders.  The T-245 picture-factory code (work/t245-picture-factory/factory245.py) is used
UNCHANGED, exactly as work/t267-trace/read267.py does (same re-pointing of its module globals to a private folder, here data/t268-trace/render/, the same shared guarded CT chunk cache, the same pinned weights and
Arc FP32 agreement files); the only differences are that the queue is data/t268-trace/render/queue.json (many patches, appended by drive268.py) and LOOKAHEAD is lifted so that produce can run before consume.
  python read268.py produce NAME[,NAME...]     render the 28-plane stacks (CPU; missing CT through the guard, 15 GB floor)
  python read268.py consume NAME[,NAME...]     read them on the Arc (arc.lock by exclusive create, FP32), save maps (float16 npz); the factory's own PNG is written too
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common268 as C  # noqa: E402

PROB = C.PROB
OUT = C.RENDER
sys.path.insert(0, str(PROB / "work/t245-picture-factory"))
sys.path.insert(0, str(PROB / "work/t237-new-candidates"))


def patched():
    import factory245 as F
    Q = {"volumes": json.loads((F.FAC / "queue.json").read_text())["volumes"], "surfaces": json.loads((OUT / "queue.json").read_text())["surfaces"]}
    F.FAC_SHARED = F.FAC
    F.FAC = OUT
    F.STACKS, F.MAPS, F.PNG, F.LOGS = OUT / "stacks", OUT / "maps", OUT / "png", OUT / "logs"
    F.STATUS, F.INDEX, F.PSTATUS = OUT / "job-status.json", OUT / "index.json", OUT / "producer-status.json"
    F.OWN_CT = F.FAC_SHARED / "ctcache"
    F.queue = lambda: Q
    F.LOOKAHEAD = 10 ** 6
    F.ME, F.AGENT = "experimenter-38", "experimenter-38 (T-268 X-212 read)"
    F.LOGS.mkdir(parents=True, exist_ok=True)
    return F


class A:
    def __init__(self, names, exp):
        self.names, self.limit, self.experiment, self.gate_chains = names, None, exp, False


if __name__ == "__main__":
    cmd, names = sys.argv[1], sys.argv[2]
    F = patched()
    if cmd == "produce":
        F.produce(A(names, "X-212"))
    elif cmd == "consume":
        F.consume(A(names, "X-212"))
