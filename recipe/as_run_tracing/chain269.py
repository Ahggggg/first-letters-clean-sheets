"""T-269 (X-213) chains.  Start each with tools/bg.py (no window).  Every stage is its own tools/run.py record under --experiment X-213.
  python chain269.py fetch SITE [SITE...]   per site: list, meta, (control: level3, block), plan, fetch, volume, seeds   (fetches run beside the tracer rounds; every guard use is under common269.guard_lock)
  python chain269.py drive SITE [SITE...]   per site: wait until it is ready, run drive269.py to its target trace count (sites269.py: target_traces), summarize
  python chain269.py extend                 sites other than box 1 that have >= 8 cm2 kept after 16 traces continue to 32 traces (registered EXTENSION rule)
  python chain269.py prep                   box 1: copy X-212's seed list into data/t269-trace/PHerc0826/ (same seeds, same order)
A STOP file (data/t269-trace/STOP) ends the drive between rounds and the chains between stages.  A site whose fetch stage failed gets data/t269-trace/<site>/FAILED and is skipped by the drive."""
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common269 as C  # noqa: E402

ROOT = C.PROB.parents[1]
VPY = str(ROOT / ".venv-xpu/Scripts/python.exe")
W = "problems/herculaneum-ink/work/t269-trace"
LOGD = C.DATA / "logs"
STOP = C.DATA / "STOP"


def say(*a):
    print(time.strftime("%Y-%m-%d %H:%M:%S"), *a, flush=True)


def run(note, *cmd, log=None):
    args = [sys.executable, "tools/run.py", "--problem", "herculaneum-ink", "--experiment", "X-213", "--timeout", "604800", "--note", note, "--"] + list(cmd)
    lf = open(log, "a") if log else subprocess.DEVNULL
    r = subprocess.run(args, cwd=str(ROOT), stdout=lf, stderr=subprocess.STDOUT)
    if log:
        lf.close()
    say("run", note, "exit", r.returncode)
    return r.returncode


def prep():
    d = C.scroll_dir("PHerc0826")
    d.mkdir(parents=True, exist_ok=True)
    src = C.DATA268 / "PHerc0826" / "seed-candidates.json"
    dst = d / "seed-candidates.json"
    if not dst.exists():
        shutil.copy2(src, dst)
    say("prep: seed list copied", dst)


def fetch_site(sc):
    d = C.scroll_dir(sc)
    d.mkdir(parents=True, exist_ok=True)
    LOGD.mkdir(parents=True, exist_ok=True)
    stages = ["list", "meta"] + (["level3", "block"] if C.site(sc).get("control") else []) + ["plan", "fetch", "volume"]
    for st in stages:
        if STOP.exists():
            say("STOP file: fetch chain ends")
            return False
        done = {"list": C.DATA / f"listing-{sc}.json", "meta": C.DATA / f"receipt-meta-{sc}.json", "level3": C.DATA / f"receipt-l3-{sc}.json", "block": d / "box.json",
                "plan": C.DATA / f"plan-{sc}.json", "volume": d / "volume.done"}.get(st)
        if done is not None and done.exists():
            continue
        say("---", sc, st)
        if run(f"fetch269 {st} {sc}", VPY, f"{W}/fetch269.py", st, sc, log=LOGD / f"fetch-{sc}-{st}.log") != 0:
            (d / "FAILED").write_text(f"{st} {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            say("FAILED", sc, st)
            return False
    if not (d / "seed-candidates.json").exists():
        say("--- seeds", sc)
        if run(f"seeds268 {sc}", VPY, f"{W}/seeds269.py", sc, "--procs", "3", log=LOGD / f"seeds-{sc}.log") != 0:
            (d / "FAILED").write_text(f"seeds {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            return False
    return True


def ready(sc):
    d = C.scroll_dir(sc)
    return (d / "seed-candidates.json").exists() and ((not C.site(sc).get("fetch")) or (d / "volume.done").exists())


def drive_site(sc, target):
    d = C.scroll_dir(sc)
    while not ready(sc):
        if STOP.exists() or (d / "FAILED").exists():
            say("site", sc, "not ready:", "STOP" if STOP.exists() else "FAILED")
            return False
        time.sleep(60)
    if STOP.exists():
        return False
    say("=== drive", sc, "target", target)
    run(f"drive269 {sc} to {target} traces", VPY, "-u", f"{W}/drive269.py", "--site", sc, "--max-traces", str(target), log=LOGD / f"drive-{sc}.log")
    run(f"summarize269 {sc}", VPY, f"{W}/summarize269.py", *[s for s in C.sites() if (C.scroll_dir(s) / "state.json").exists()], log=LOGD / "summarize.log")
    return True


def kept_cm2(sc):
    f = C.scroll_dir(sc) / "state.json"
    if not f.exists():
        return 0.0, 0
    st = json.loads(f.read_text())
    return sum(p["cm2"] for p in st["patches"]), sum(1 for t in st["traces"].values() if t.get("status") not in ("skipped_on_traced_sheet", "round_error"))


if __name__ == "__main__":
    cmd, sites = sys.argv[1], sys.argv[2:]
    if cmd == "prep":
        prep()
    elif cmd == "fetch":
        for sc in sites:
            if STOP.exists():
                break
            fetch_site(sc)
        say("fetch chain finished")
    elif cmd == "drive":
        prep()
        for sc in sites:
            drive_site(sc, int(C.site(sc)["target_traces"]))
        say("drive chain finished")
    elif cmd == "extend":
        for sc, s in sorted(C.sites().items(), key=lambda kv: kv[1]["order"]):
            if sc == "PHerc0826" or STOP.exists():
                continue
            cm2, n = kept_cm2(sc)
            if n >= 16 and cm2 >= 8.0 and n < 32:
                say("extension rule met for", sc, round(cm2, 2), "cm2 after", n, "traces")
                drive_site(sc, 32)
        say("extend finished")
