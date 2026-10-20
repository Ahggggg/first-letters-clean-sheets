#!/usr/bin/env python3
"""fetch_mirror.py - download the part of the public m7 surface prediction (level 0) and its normal grids that one tracer box needs.  MIT licence.

  python fetch_mirror.py PHerc0813 20250821151723 --box XLO XHI YLO YHI ZLO ZHI --out MIRROR_DIR [--dry-run]

Layout written (what vc_grow_seg_from_seed needs):  MIRROR_DIR/pred.zarr/{.zgroup,.zattrs,metadata.json,0/.zarray,0/<cz>/<cy>/<cx>}  and  MIRROR_DIR/ngrid/{metadata.json,xy|xz|yz/NNNNNN.grid}.
Chunks that the bucket does not list are the zarr fill value 0 and are simply not written.  Chunk = 192^3 voxels (blosc); a normal-grid file is one whole slice (xy: z index, xz: y index, yz: x index).
Typical size: 0.4-1.5 GB of prediction chunks plus 3-5 GB of normal-grid slices for a 3.5-5 k-voxel box.  Every GET is checked against its listed size (Content-Length); existing files of the right size are kept.

STATUS: adapted from the project's fetch267.py / fetch268.py, which downloaded 3.2 GB (w035 box), 6.2 GB (PHerc0826 box) and the five later mirror blocks this way, with the project's own download meter and
SHA-256 receipts added.  THIS CLEANED COPY (no meter, no receipts) was not re-run against the bucket when the repository was assembled.  Public data: Vesuvius Challenge open data (CC BY-NC 4.0).
"""
import argparse
import json
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

S3 = "https://vesuvius-challenge-open-data.s3.amazonaws.com/"
CH = 192


def base(scroll, volume):
    return f"{scroll}/representations/predictions/surfaces/{volume}-surface-20260413222639-surface-m7-L0-th0.2"


def listing(prefix):
    """[(key, size)] under prefix (S3 ListObjectsV2, all pages)."""
    out, token = [], None
    while True:
        q = {"list-type": "2", "prefix": prefix}
        if token:
            q["continuation-token"] = token
        with urllib.request.urlopen(S3 + "?" + urllib.parse.urlencode(q), timeout=120) as r:
            root = ET.fromstring(r.read())
        ns = {"s": root.tag.split("}")[0].strip("{")}
        for c in root.findall("s:Contents", ns):
            out.append((c.find("s:Key", ns).text, int(c.find("s:Size", ns).text)))
        if root.find("s:IsTruncated", ns).text == "true":
            token = root.find("s:NextContinuationToken", ns).text
        else:
            return out


def get(key, size, dest):
    if dest.exists() and dest.stat().st_size == size:
        return 0
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(4):
        try:
            req = urllib.request.Request(S3 + key, headers={"Accept-Encoding": "identity"})
            with urllib.request.urlopen(req, timeout=120) as r:
                if int(r.headers["Content-Length"]) != size:
                    raise OSError(f"{key}: Content-Length {r.headers['Content-Length']} != listed {size}")
                data = r.read()
            if len(data) != size:
                raise OSError(f"{key}: got {len(data)} of {size}")
            tmp = Path(str(dest) + ".part")
            tmp.write_bytes(data)
            tmp.replace(dest)
            return size
        except Exception as e:  # noqa: BLE001
            if attempt == 3:
                raise
            time.sleep(2 + 4 * attempt)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("scroll")
    ap.add_argument("volume")
    ap.add_argument("--box", nargs=6, type=int, required=True, metavar=("XLO", "XHI", "YLO", "YHI", "ZLO", "ZHI"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    b = base(a.scroll, a.volume)
    out = Path(a.out)
    xlo, xhi, ylo, yhi, zlo, zhi = a.box
    chunks = {tuple(int(v) for v in k[len(b + ".zarr/0/"):].split("/")): s for k, s in listing(b + ".zarr/0/") if k[len(b + ".zarr/0/"):][0].isdigit()}
    want = [(z, y, x) for z in range(zlo // CH, (zhi - 1) // CH + 1) for y in range(ylo // CH, (yhi - 1) // CH + 1) for x in range(xlo // CH, (xhi - 1) // CH + 1) if (z, y, x) in chunks]
    objs = [(f"{b}.zarr/0/{z}/{y}/{x}", chunks[(z, y, x)], out / "pred.zarr/0" / str(z) / str(y) / str(x)) for z, y, x in want]
    for k, s in listing(b + ".zarr/"):
        n = k[len(b + ".zarr/"):]
        if n in (".zgroup", ".zattrs", "metadata.json", "0/.zarray"):
            objs.append((k, s, out / "pred.zarr" / n))
    ranges = {"xy": (zlo, zhi), "xz": (ylo, yhi), "yz": (xlo, xhi)}
    for k, s in listing(b + ".normal-grids/"):
        n = k[len(b + ".normal-grids/"):]
        if n == "metadata.json":
            objs.append((k, s, out / "ngrid" / n))
        else:
            d, f = n.split("/")
            lo, hi = ranges[d]
            if lo <= int(f.split(".")[0]) < hi:
                objs.append((k, s, out / "ngrid" / n))
    total = sum(s for _, s, _ in objs)
    print(f"{len(objs)} objects, {total / 2 ** 30:.2f} GiB ({len(want)} of the box's chunks exist)")
    if a.dry_run:
        return
    done = 0
    with ThreadPoolExecutor(8) as pool:
        for n in pool.map(lambda o: get(*o), objs):
            done += n
    print("downloaded", round(done / 2 ** 30, 2), "GiB")


if __name__ == "__main__":
    main()
