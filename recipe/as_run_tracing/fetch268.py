"""T-268 (X-212) data for the NEXT scrolls: listing, level-3 pyramid, mirror block choice, planned fetch and execution - every byte through the UNCHANGED download guard via work/t267-trace/fetch267.py code
(listed sizes asserted against each GET Content-Length, SHA-256 receipts, packs of <= 4000 objects, 15 GB free-space floor).  fetch267's module globals (output folder, scroll -> volume id) are re-pointed to
data/t268-trace/; nothing else is changed.
  python fetch268.py list   SCROLL            S3 ListObjectsV2 of level-0 chunk keys and normal-grid files (metadata, not metered) + level 3 keys  -> data/t268-trace/listing-SCROLL.json
  python fetch268.py meta   SCROLL            small metadata objects (guarded)
  python fetch268.py level3 SCROLL            the level-3 prediction chunks (about 0.1 GB) + its .zarray (guarded)
  python fetch268.py block  SCROLL            choose the mirror block from level 3 (registered rule: the 3456-voxel cube of 18 x 18 x 18 chunks inside the volume with the largest sum of level-3 values) -> SCROLL/box.json
  python fetch268.py plan   SCROLL            plan of the level-0 chunks + normal-grid slices of the block -> plan-SCROLL.json (no download)
  python fetch268.py fetch  SCROLL            execute the plan (guarded)
  python fetch268.py volume SCROLL            copy the mirror into the Docker named volume t268-<scroll> (helper container, removed afterwards) and write pred.zarr/meta.json
"""
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common268 as C  # noqa: E402
import fetch267 as F  # noqa: E402

F.OUT = C.DATA
for sc in C.queue_volumes():
    F.SC[sc] = C.queue_volumes()[sc]["volume"]
BLOCK_CH = 18


def cmd_list(sc):
    p = F.C.load_r119().p
    b = F.base(sc)
    out = {"base": b}
    _, k = p.listing(b + ".zarr/0/", delimiter=False)
    out["chunks"] = [[x[0][len(b + ".zarr/0/"):], x[1]] for x in k]
    for d in ("xy", "xz", "yz"):
        _, k = p.listing(b + f".normal-grids/{d}/", delimiter=False)
        out["ng_" + d] = [[x[0].rsplit("/", 1)[1], x[1]] for x in k]
    _, k = p.listing(b + ".zarr/3/", delimiter=False)
    out["l3"] = [[x[0][len(b + ".zarr/3/"):], x[1]] for x in k]
    C.jwrite(C.DATA / f"listing-{sc}.json", out)
    print(sc, "level-0 chunks", len(out["chunks"]), round(sum(s for _, s in out["chunks"]) / 2 ** 30, 2), "GiB; normal grids",
          {d: (len(out["ng_" + d]), round(sum(s for _, s in out["ng_" + d]) / 2 ** 30, 2)) for d in ("xy", "xz", "yz")}, "; level 3", len(out["l3"]), round(sum(s for _, s in out["l3"]) / 2 ** 20, 1), "MiB")


def cmd_meta(sc):
    F.cmd_meta(sc)


def cmd_level3(sc):
    lst = json.loads((C.DATA / f"listing-{sc}.json").read_text())
    b = F.base(sc)
    objs = [(b + ".zarr/3/" + k, s) for k, s in lst["l3"]]
    rec, staged, dg, retries, size, m0, m1 = F.guarded_pack(f"l3-{sc}", objs, b)
    receipts = F.split(objs, staged, dg, sc)
    C.jwrite(C.DATA / f"receipt-l3-{sc}.json", {"bytes": size, "meter_before": m0, "meter_after": m1, "pack_sha256": rec["dest_sha256"], "retries": retries, "receipts": receipts})
    print("level3", sc, size, "bytes; metered", m1 - m0)


def cmd_block(sc):
    import numcodecs
    vi = C.vol_info(sc)
    zarray = json.loads((C.DATA / sc / "pred.zarr/3/.zarray").read_text())
    shape3, ch3 = zarray["shape"], zarray["chunks"]
    codec = numcodecs.Blosc()
    sub = 24                                                    # level-3 voxels per level-0 chunk (192 / 8)
    nz, ny, nx = (int(np.ceil(vi["shape"][i] / C.CH)) for i in range(3))          # volume extent in level-0 chunks (z, y, x)
    occ = np.zeros((nz, ny, nx), np.float64)
    for f in (C.DATA / sc / "pred.zarr/3").glob("*/*/*"):
        cz, cy, cx = (int(v) for v in f.parts[-3:])
        a = np.frombuffer(codec.decode(f.read_bytes()), np.uint8)
        a = a.reshape(ch3[0], ch3[1], ch3[2]).astype(np.float64)
        zz, yy, xx = (min(ch3[i], shape3[i] - (cz, cy, cx)[i] * ch3[i]) for i in range(3))
        a = a[:zz, :yy, :xx]
        pad = [(0, (-s) % sub) for s in a.shape]
        a = np.pad(a, pad)
        a = a.reshape(a.shape[0] // sub, sub, a.shape[1] // sub, sub, a.shape[2] // sub, sub).sum((1, 3, 5))
        z0, y0, x0 = cz * ch3[0] // sub, cy * ch3[1] // sub, cx * ch3[2] // sub
        sz = occ[z0:z0 + a.shape[0], y0:y0 + a.shape[1], x0:x0 + a.shape[2]].shape
        occ[z0:z0 + sz[0], y0:y0 + sz[1], x0:x0 + sz[2]] += a[:sz[0], :sz[1], :sz[2]]
    from scipy.ndimage import uniform_filter
    n = BLOCK_CH
    cs = uniform_filter(occ, n, mode="constant") * n ** 3        # sum over the window centred at each chunk
    lo = n // 2
    best, arg = -1.0, None
    for z in range(lo, nz - (n - lo) + 1):
        for y in range(lo, ny - (n - lo) + 1):
            row = cs[z, y, lo:nx - (n - lo) + 1]
            if row.size and row.max() > best:
                best, arg = float(row.max()), (z, y, lo + int(row.argmax()))
    z, y, x = arg
    zc0, yc0, xc0 = z - lo, y - lo, x - lo
    box = [xc0 * C.CH, (xc0 + n) * C.CH, yc0 * C.CH, (yc0 + n) * C.CH, zc0 * C.CH, (zc0 + n) * C.CH]
    box = [box[0], min(box[1], vi["shape"][2]), box[2], min(box[3], vi["shape"][1]), box[4], min(box[5], vi["shape"][0])]
    C.jwrite(C.DATA / sc / "box.json", {"scroll": sc, "box_xyz": box, "chunk_origin_zyx": [zc0, yc0, xc0], "window_sum": best, "rule": "3456-voxel cube (18 chunks) inside the volume with the largest sum of level-3 values",
                                       "level3_total": float(occ.sum())})
    print("block", sc, box, "window sum", best, "share of all level-3 mass", best / max(occ.sum(), 1))


def cmd_plan(sc):
    box = C.box_of(sc)
    F.cmd_plan(sc, sc, *box)


def cmd_fetch(sc):
    F.cmd_fetch(sc)


def cmd_volume(sc):
    vi = C.vol_info(sc)
    src = C.DATA / sc
    vol = C.docker_volume(sc)
    meta = {"uuid": f"t268-{sc}-{vi['volume']}-m7", "name": f"{sc} {vi['volume']} m7 surface prediction (T-268 local mirror)", "type": "vol", "format": "zarr", "width": vi["shape"][2], "height": vi["shape"][1],
            "slices": vi["shape"][0], "voxelsize": vi["voxel_um"]}
    (src / "pred.zarr" / "meta.json").write_bytes(json.dumps(meta).encode())
    import os
    os.environ["PATH"] = C.DOCKER_BIN.replace("/", "\\") + os.pathsep + os.environ.get("PATH", "")
    subprocess.run(["docker", "volume", "create", vol], check=True)
    t0 = time.time()
    script = "set -e; mkdir -p /data/pred.zarr /data/ngrid; cp -r /mirror/pred.zarr/0 /data/pred.zarr/; cp /mirror/pred.zarr/meta.json /data/pred.zarr/; cp -r /mirror/ngrid/. /data/ngrid/; du -sh /data/pred.zarr /data/ngrid; ls /data/ngrid"
    r = subprocess.run(["docker", "run", "--rm", "-v", f"{vol}:/data", "-v", f"{src}:/mirror:ro", "--entrypoint", "sh", C.IMAGE, "-c", script], capture_output=True, text=True)
    print(r.stdout[-1500:], r.stderr[-800:], "seconds", round(time.time() - t0))
    if r.returncode != 0:
        raise SystemExit("volume copy failed")
    (src / "volume.done").write_text(time.strftime("%Y-%m-%d %H:%M:%S") + "\n")


if __name__ == "__main__":
    {"list": cmd_list, "meta": cmd_meta, "level3": cmd_level3, "block": cmd_block, "plan": cmd_plan, "fetch": cmd_fetch, "volume": cmd_volume}[sys.argv[1]](sys.argv[2])
