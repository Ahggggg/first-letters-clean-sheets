"""support276.py (T-276, POST-HOC, not registered in X-218): does each kept X-213 patch lie on the public m7 prediction it was grown on?  The UNCHANGED support268 logic (as X-212's post-hoc check and as named in X-213's
registration but never run by that chain), on the T-269 sites, written to data/t276-xsec/support_t269.json (nothing is written into data/t269-trace/).  CPU only.  python support276.py"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROB = HERE.parents[1]
sys.path.insert(0, str(PROB / "work/t269-trace"))
import common269 as C  # noqa: E402  (installs the common268 alias)
sys.path.insert(0, str(PROB / "work/t268-trace"))
import support268 as S  # noqa: E402

out = {}
for site in C.sites():
    d = C.scroll_dir(site) / "state.json"
    if not d.exists():
        continue
    st = json.loads(d.read_text())
    if not st["patches"]:
        continue
    sm = S.Sampler(site)
    for p in st["patches"]:
        r = S.patch_support(sm, p["mesh"])
        out[p["name"]] = {**r, "site": site}
        print(p["name"], round(p["cm2"], 2), "on", round(r["on_prediction_share"], 3), "near", round(r["near_prediction_share_pm4"], 3), flush=True)
(PROB / "data/t276-xsec/support_t269.json").write_text(json.dumps(out, indent=1))
