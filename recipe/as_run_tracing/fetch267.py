"""T-267 data fetch: the public surface-prediction zarr chunks and normal-grid slices a tracer run needs (and nothing else), every byte through the UNCHANGED
download guard (work/download_guard_d005.py guarded_download, continuing meter configs/transfer-meter-d005.json, 15 GB free-space floor on the data drive), packs of
<= 4000 objects, listed sizes asserted against each GET Content-Length, SHA-256 receipt per object.  No other fetch path.
  python fetch267.py meta  <scroll>                                   small metadata objects (.zgroup .zattrs metadata.json 0/.zarray, normal-grids/metadata.json)
  python fetch267.py plan  <name> <scroll> xlo xhi ylo yhi zlo zhi    voxel box (level 0, xyz) -> plan json (chunks present on S3 + grid slices), no download
  python fetch267.py fetch <name>                                     execute the plan
Local layout: data/t267-trace/<scroll>/pred.zarr/{.zgroup,.zattrs,metadata.json,0/.zarray,0/z/y/x} and data/t267-trace/<scroll>/ngrid/{metadata.json,xy|xz|yz/NNNNNN.grid}.
Chunks absent from the S3 listing are the zarr fill value (0); they are simply not written."""
import collections
import hashlib
import itertools
import json
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROB = HERE.parents[1]
OUT = PROB / "data/t267-trace"
sys.path.insert(0, str(PROB / "work/t261-single-sheet"))
import common261 as C  # noqa: E402

S3 = "https://vesuvius-challenge-open-data.s3.amazonaws.com/"
SC = {"PHerc0826": "20250821151701", "PHerc0139": "20250728140407"}
CH = 192


def base(sc):
    return f"{sc}/representations/predictions/surfaces/{SC[sc]}-surface-20260413222639-surface-m7-L0-th0.2"


def dest_of(sc, key):
    b = base(sc)
    if key.startswith(b + ".zarr/"):
        return OUT / sc / "pred.zarr" / key[len(b + ".zarr/"):]
    return OUT / sc / "ngrid" / key[len(b + ".normal-grids/"):]


def guarded_pack(name, objs, tag):
    """objs: list of (s3 key, listed size)."""
    r119 = C.load_r119()
    acq = r119.acq
    cfg = acq.load(acq.CONFIG)
    staging = OUT / "staging"
    staging.mkdir(parents=True, exist_ok=True)
    staged = staging / f"{name}.bin"
    assert not staged.exists() and not Path(str(staged) + ".partial").exists(), ("leftover staging", staged)
    size = sum(s for _, s in objs)
    digests, retries = [], []

    def get(ks):
        key, sz = ks
        req = urllib.request.Request(S3 + key, headers={"Accept-Encoding": "identity", "User-Agent": "Lacuna-T267"})
        for attempt in range(3):
            try:
                with urllib.request.urlopen(req, timeout=120) as resp:
                    assert int(resp.headers["Content-Length"]) == sz, (key, resp.headers["Content-Length"], sz)
                    data = resp.read()
                assert len(data) == sz
                return data
            except Exception as e:  # noqa: BLE001
                retries.append({"key": key, "attempt": attempt, "error": type(e).__name__})
                if attempt == 2:
                    raise
                time.sleep(2 + 4 * attempt)

    def blocks():
        with ThreadPoolExecutor(8) as pool:
            pending, source = collections.deque(), iter(objs)
            for ks in itertools.islice(source, 16):
                pending.append(pool.submit(get, ks))
            while pending:
                data = pending.popleft().result()
                ks = next(source, None)
                if ks is not None:
                    pending.append(pool.submit(get, ks))
                digests.append(hashlib.sha256(data).hexdigest())
                for st in range(0, len(data), 1 << 20):
                    yield data[st:st + (1 << 20)]

    meter_before = acq.meter.total(acq.METER)
    rec = acq.guard.guarded_download(
        item={"id": f"T267-{name}", "root": "C", "kind": "download", "basis": "listed_object_bytes",
              "listing_source": "S3 ListObjectsV2 sizes per object; each GET Content-Length asserted", "bytes": size,
              "restorable": {"kind": "refetch", "source": f"{S3}{tag} (plan json)"}},
        dest=staged, allowed_root=staging, meter_path=acq.METER, cfg=cfg, chunks=blocks(),
        c_drive=acq.storage.measure_drive(acq.ROOT), d_drive=acq.storage.measure_drive(cfg["d_data_root"]),
        d_args=acq.d_args(cfg), project_bytes_c=acq.fast_project_bytes(), project_measure=acq.fast_project_bytes)
    assert rec["status"] == "completed", (name, rec["status"], rec.get("reason"))
    assert acq.meter.total(acq.METER) - meter_before == size and len(digests) == len(objs)
    return rec, staged, digests, retries, size, meter_before, acq.meter.total(acq.METER)


def split(objs, staged, digests, sc):
    raw = staged.read_bytes()
    assert len(raw) == sum(s for _, s in objs)
    off, rec = 0, {}
    for (key, sz), dg in zip(objs, digests):
        blob = raw[off:off + sz]
        off += sz
        assert hashlib.sha256(blob).hexdigest() == dg
        f = dest_of(sc, key)
        f.parent.mkdir(parents=True, exist_ok=True)
        tmp = Path(str(f) + ".part")
        tmp.write_bytes(blob)
        tmp.replace(f)
        rec[key[len(base(sc)):]] = [sz, dg]
    staged.unlink()
    return rec


def cmd_meta(sc):
    b = base(sc)
    cand = [b + ".zarr/.zgroup", b + ".zarr/.zattrs", b + ".zarr/metadata.json", b + ".zarr/0/.zarray", b + ".normal-grids/metadata.json"]
    p = C.load_r119().p
    sizes = {}
    for pre in (b + ".zarr/", b + ".zarr/0/", b + ".normal-grids/"):
        _, k = p.listing(pre)
        sizes.update(dict(k))
    objs = [(k, sizes[k]) for k in cand]
    rec, staged, dg, retries, size, m0, m1 = guarded_pack(f"meta-{sc}", objs, b)
    receipts = split(objs, staged, dg, sc)
    (OUT / f"receipt-meta-{sc}.json").write_text(json.dumps({"bytes": size, "meter_before": m0, "meter_after": m1, "pack_sha256": rec["dest_sha256"],
                                                              "retries": retries, "receipts": receipts}, indent=1))
    print("meta", sc, size, "bytes; metered", m1 - m0)


def cmd_plan(name, sc, xlo, xhi, ylo, yhi, zlo, zhi):
    xlo, xhi, ylo, yhi, zlo, zhi = map(int, (xlo, xhi, ylo, yhi, zlo, zhi))
    b = base(sc)
    lst = json.loads((OUT / f"listing-{sc}.json").read_text())
    ck = {tuple(map(int, k.split("/"))): s for k, s in lst["chunks"] if k != ".zarray"}
    want = [(z, y, x) for z in range(zlo // CH, (zhi - 1) // CH + 1) for y in range(ylo // CH, (yhi - 1) // CH + 1) for x in range(xlo // CH, (xhi - 1) // CH + 1)]
    have = [k for k in want if k in ck]
    objs = [(b + ".zarr/0/" + "/".join(map(str, k)), ck[k]) for k in have]
    rngs = {"xy": (zlo, zhi), "xz": (ylo, yhi), "yz": (xlo, xhi)}
    n_grid = {}
    for d, (lo, hi) in rngs.items():
        g = {int(k.split(".")[0]): s for k, s in lst["ng_" + d]}
        sel = [i for i in range(lo, hi) if i in g]
        objs += [(b + f".normal-grids/{d}/{i:06d}.grid", g[i]) for i in sel]
        n_grid[d] = [len(sel), sum(g[i] for i in sel)]
    plan = {"name": name, "scroll": sc, "box_xyz": [xlo, xhi, ylo, yhi, zlo, zhi], "chunks_wanted": len(want), "chunks_on_s3": len(have),
            "chunk_bytes": sum(ck[k] for k in have), "grids": n_grid, "total_bytes": sum(s for _, s in objs), "objects": objs}
    (OUT / f"plan-{name}.json").write_text(json.dumps(plan))
    print(name, {k: v for k, v in plan.items() if k != "objects"}, f"{plan['total_bytes'] / 2**30:.2f} GiB")


def cmd_fetch(name):
    plan = json.loads((OUT / f"plan-{name}.json").read_text())
    sc = plan["scroll"]
    objs = [tuple(o) for o in plan["objects"] if not dest_of(sc, o[0]).exists()]
    print("objects to fetch", len(objs), "of", len(plan["objects"]), flush=True)
    for i in range(0, len(objs), 4000):
        part = objs[i:i + 4000]
        rec, staged, dg, retries, size, m0, m1 = guarded_pack(f"{name}-{i // 4000:03d}", part, base(sc))
        receipts = split(part, staged, dg, sc)
        info = {"objects": len(part), "bytes": size, "meter_before": m0, "meter_after": m1, "pack_sha256": rec["dest_sha256"], "retries": retries}
        (OUT / f"receipt-{name}-{i // 4000:03d}.json").write_text(json.dumps({**info, "receipts": receipts}))
        print(json.dumps(info), flush=True)
    print("done", name)


if __name__ == "__main__":
    {"meta": cmd_meta, "plan": cmd_plan, "fetch": cmd_fetch}[sys.argv[1]](*sys.argv[2:])
