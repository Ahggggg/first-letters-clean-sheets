#!/usr/bin/env python3
"""trace_round.py - run villa's vc_grow_seg_from_seed from several seeds in ONE Docker container, from a Docker named volume (not a bind mount).  MIT licence.

  python trace_round.py volume MIRROR_DIR --name t-0813 --voxel-um 9.362 --shape Z Y X        # once per mirror: write meta.json, create the volume, copy the mirror into it
  python trace_round.py round --name t-0813 --voxel-um 9.362 --out OUT_DIR --seeds "x y z" "x y z" ...   [--generations 75] [--cpus 5] [--memory 14g]

Why a volume: through a Windows bind mount the tracer sat in state D at about 0.3 of a core (25 generations in 28 minutes, aborted); from a named volume (a copy on the Docker VM's own disk;
3-6 GB took 6.6 min to copy) 75 generations of the official recipe take 14-20 minutes per tracer with 4 tracers side by side on 5 CPUs.
Why meta.json: the tracer reads the volume root's meta.json {type vol, uuid, width, height, slices, voxelsize}.  Without a voxel size it computes area 0, applies min_area_cm and DELETES the finished
output; a public metadata.json without a 'scan' key makes it throw.  So the voxel size (9.362 for the PHerc0125/0191/0211/0257/0358/0813/0826/0846/1203/1545 scans) is written there.
Entry point: the image's default entry point wants AWS credentials; it is overridden (--entrypoint) so the container needs no network and no credentials.
Image: ghcr.io/scrollprize/villa/volume-cartographer:main (the runs used image sha256:bad516f66001abca759454cc43e4fd11e5b19aa55d36bdc2043817291c8083c4, built 2026-05-13).
Params = the official PHerc0826 auto_grown recipe: mode seed, generations 75, thread_limit 1, use_cuda false, voxelsize, min_area_cm 0.3, cache_size 3e9.  Two runs of one seed differ by about 2.4 voxels
(median), so traces are not bit-exact.  The side convention of the output grid (sign of the frame flag) depends on the seed: always measure it (code/tifxyz_audit.py).
Output: OUT_DIR/<tid>/seg/{x,y,z,generations}.tif + meta.json per seed (tid = s000, s001, ...) and OUT_DIR/<tid>.log.

STATUS: adapted from the project's fetch268.py (volume) and drive269.py (round), which ran this 36 rounds; THIS CLEANED COPY was not re-run (a tracer job was using the machine when it was assembled).
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

IMAGE = "ghcr.io/scrollprize/villa/volume-cartographer:main"


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def cmd_volume(a):
    src = Path(a.mirror).resolve()
    z, y, x = a.shape
    meta = {"uuid": f"{a.name}-m7", "name": f"{a.name} m7 surface prediction (local mirror)", "type": "vol", "format": "zarr", "width": x, "height": y, "slices": z, "voxelsize": a.voxel_um}
    (src / "pred.zarr" / "meta.json").write_text(json.dumps(meta))
    vol = a.name
    subprocess.run(["docker", "volume", "create", vol], check=True)
    script = ("set -e; mkdir -p /data/pred.zarr /data/ngrid; cp -r /mirror/pred.zarr/0 /data/pred.zarr/; cp /mirror/pred.zarr/meta.json /data/pred.zarr/; "
              "cp -r /mirror/ngrid/. /data/ngrid/; du -sh /data/pred.zarr /data/ngrid")
    r = run(["docker", "run", "--rm", "-v", f"{vol}:/data", "-v", f"{src}:/mirror:ro", "--entrypoint", "sh", IMAGE, "-c", script])
    print(r.stdout[-1500:], r.stderr[-800:])
    if r.returncode:
        sys.exit("volume copy failed")


def cmd_round(a):
    out = Path(a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    params = {"cache_root": "/tmp/cache", "cache_size": 3000000000, "generations": a.generations, "min_area_cm": 0.3, "mode": "seed", "normal_grid_path": "/data/ngrid", "thread_limit": 1,
              "use_cuda": False, "voxelsize": a.voxel_um}
    (out / "params.json").write_text(json.dumps(params))
    lines = []
    for i, s in enumerate(a.seeds):
        tid = f"s{i:03d}"
        (out / tid).mkdir(exist_ok=True)
        lines.append(f"{tid} {s}")
    (out / "seeds.txt").write_text("\n".join(lines) + "\n")
    (out / "run.sh").write_text("#!/bin/sh\nwhile read tid x y z; do\n"
                                "  ( vc_grow_seg_from_seed -v /data/pred.zarr -t /out/$tid/seg -p /params.json --seed $x $y $z --segment-name seg > /out/$tid.log 2>&1 ; echo exit $? >> /out/$tid.log ) &\n"
                                "done < /seeds.txt\nwait\n")
    name = f"{a.name}-round"
    run(["docker", "rm", "-f", name])
    cmd = ["docker", "run", "-d", "--name", name, "--cpus", str(a.cpus), "--memory", a.memory, "-v", f"{a.name}:/data:ro", "-v", f"{out / 'params.json'}:/params.json:ro",
          "-v", f"{out / 'seeds.txt'}:/seeds.txt:ro", "-v", f"{out / 'run.sh'}:/run.sh:ro", "-v", f"{out}:/out", "--entrypoint", "sh", IMAGE, "/run.sh"]
    r = run(cmd)
    if r.returncode:
        sys.exit("docker run failed: " + r.stderr[-500:])
    print("started", name, "with", len(a.seeds), "tracers; waiting")
    run(["docker", "wait", name])
    print(run(["docker", "inspect", name, "--format", "exit {{.State.ExitCode}} oom {{.State.OOMKilled}}"]).stdout.strip())
    run(["docker", "rm", name])


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("volume")
    v.add_argument("mirror")
    v.add_argument("--name", required=True)
    v.add_argument("--voxel-um", type=float, required=True)
    v.add_argument("--shape", nargs=3, type=int, required=True, metavar=("Z", "Y", "X"))
    v.set_defaults(f=cmd_volume)
    r = sub.add_parser("round")
    r.add_argument("--name", required=True)
    r.add_argument("--voxel-um", type=float, required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--seeds", nargs="+", required=True, help='each "x y z" in level-0 voxels')
    r.add_argument("--generations", type=int, default=75)
    r.add_argument("--cpus", type=float, default=5)
    r.add_argument("--memory", default="14g")
    r.set_defaults(f=cmd_round)
    a = ap.parse_args()
    a.f(a)
