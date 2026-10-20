#!/usr/bin/env python3
"""xsec_patch.py - CT cross-sections of one patch and the layer-tilt measure of the README ("CT cross-sections").  MIT licence.

  python xsec_patch.py data/PHerc0813/PHerc0813_t269_s000.tifxyz --ct CT_CHUNK_DIR [--ct DIR2 ...] [--npz sections.npz]

CT_CHUNK_DIR = level-0 chunks of the scan as raw uint8 files <dir>/<cz>/<cy>/<cx> (128^3 voxels each; the public volume's level-0 zarr is stored uncompressed in exactly this layout, so a local copy of the
needed chunks, or a mount of the bucket, works).  Chunks that are not present are treated as unknown and the columns that need them are marked invalid (never read as zero).
Per patch: 5 straightened sections along each direction (canvas rows / columns at 10, 30, 50, 70, 90 % of the patch), depth -32..+32 voxels along the surface normal (+ = away from the scroll axis),
tilt (degrees; 0 = layer structure parallel to the surface), coherence, noise threshold (99th percentile of 20 fully shuffled copies, numpy seed 276), implied drift = L tan(tilt) in voxels.
The tilt of a section is the dominant orientation of the structure tensor (Gaussian-derivative gradients, sigma 2 voxels, window |d| <= 24 voxels) pooled over the section; see xsec.py.
Note: the noise-null 'informative' flag was true for every one of the 900 patch sections; it does not screen anything.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import xsec as X  # noqa: E402

DEPTHS = np.arange(-32, 33)
FRACS = (0.1, 0.3, 0.5, 0.7, 0.9)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tifxyz")
    ap.add_argument("--ct", action="append", required=True)
    ap.add_argument("--npz")
    ap.add_argument("--json")
    a = ap.parse_args()
    geom = X.load_tifxyz(a.tifxyz)
    meta = geom["meta"]
    shape = meta["volume_shape_zyx"]
    vol = X.ChunkVolume(a.ct, shape)
    side, minority = X.frame_side(geom, shape)
    rng = np.random.default_rng(276)
    H, W = X.canvas_size(geom)
    secs, arrays = [], {}
    for kind in ("row", "col"):
        n_idx, n_s = (H, W) if kind == "row" else (W, H)
        for f in FRACS:
            idx = int(round(f * (n_idx - 1)))
            S, valid, L, _ = X.section(vol, geom, kind, idx, (0, n_s), DEPTHS, side)
            m = X.measure(S, valid, DEPTHS, rng)
            rec = {"kind": kind, "index": idx, "frac": f, "arclen_vox": L, **m}
            if np.isfinite(m["tilt"]):
                rec["drift_vox"] = float(L * np.tan(np.radians(m["tilt"])))
            secs.append(rec)
            arrays[f"{kind}_{idx}"] = S.astype(np.float16)
            arrays[f"{kind}_{idx}_valid"] = valid
            print(kind, idx, f"tilt {m['tilt']:+.2f} deg  coherence {m['coherence']:.2f}  drift {rec.get('drift_vox', float('nan')):+.0f} vox  uncached columns {int((~valid).sum())}")
    if a.npz:
        np.savez_compressed(a.npz, **arrays)
    if a.json:
        Path(a.json).write_text(json.dumps({"side": side, "frame_minority": minority, "sections": secs}, indent=1))


if __name__ == "__main__":
    main()
