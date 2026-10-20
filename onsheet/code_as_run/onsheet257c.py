"""T-257 / X-201: X-199 with the depth term replaced by C = depth-profile contrast of the tile (percent). Imports onsheet257 UNCHANGED and replaces only score_row's B output."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import onsheet257 as m  # noqa: E402

_orig = m.score_row


def score_row_c(st, valid, r0, nx, rng=None, null=None):
    r = _orig(st, valid, r0, nx, rng, "spatial" if null == "spatial" else None)
    T, BLK = m.T, m.BLK
    nb = T // BLK
    vb = valid[r0:r0 + T, :nx * T].reshape(nb, BLK, nx * nb, BLK)
    cnt = vb.sum((1, 3)).astype(np.float32)
    bm = np.zeros((st.shape[0], nb, nx * nb), np.float32)
    for k in range(st.shape[0]):
        pk = np.asarray(st[k, r0:r0 + T, :nx * T]).astype(np.float32) * valid[r0:r0 + T, :nx * T]
        bm[k] = pk.reshape(nb, BLK, nx * nb, BLK).sum((1, 3)) / np.maximum(cnt, 1)
    if null == "depth":
        bm = rng.permuted(bm, axis=0)
    good = (cnt >= 0.5 * BLK * BLK).astype(np.float32) * cnt
    C = np.full(nx, np.nan, np.float32)
    for j in range(nx):
        sl = slice(j * nb, (j + 1) * nb)
        w = good[:, sl]
        if (w > 0).sum() >= 0.5 * nb * nb:
            prof = (bm[:, :, sl] * w).sum((1, 2)) / w.sum()
            C[j] = (prof.max() - prof.min()) / prof.mean() * 100
    r["B"] = C
    return r


m.score_row = score_row_c
m.OUT = m.FAC / "onsheet_c"
m.TILES = m.OUT / "tiles"
if __name__ == "__main__":
    m.main()
