"""T-257 / X-199: on-sheet quality of the 120 picture-factory surfaces (descriptive, stored stacks only, CPU, no Arc, no model).
Per 256 x 256 canvas-px tile (>= 0.9 valid):
  A = fourth-harmonic orientation coherence of plane 14 (Gaussian sigma 1.5, w = |grad|^2): |sum w e^{4 i theta}| / sum w
  B = share of 8x8 blocks whose bright-layer depth (argmax over the 28 planes of the [1,2,1]-smoothed block mean) is within 2 planes of the tile median.
Tile on-sheet iff A >= tA and B >= tB; tA, tB = 10th percentile over valid tiles of known-writing stack w00 (nothing else).
Order: (1) w00 -> thresholds.json (frozen, sha256) (2) nulls on w00 (3) held-out known writing a144 a174 p500 (4) 120 factory surfaces (5) ranking + pictures of the top 20.
Seed 20261257 for the two null shuffles only.
python onsheet257.py [--only NAME ...] [--skip-pictures]"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

HERE = Path(__file__).resolve().parent
PROB = HERE.parents[1]
FAC = PROB / "data/picture-factory"
OUT = FAC / "onsheet"
TILES = OUT / "tiles"
REFDIR = PROB / "data/t240-heldout/stacks"
sys.path.insert(0, str(PROB / "work/t245-picture-factory"))
import pic245  # noqa: E402

SEED = 20261257
T = 256
BLK = 8
PLANE = 14
SIGMA = 1.5
PAD = 5
VALID_TILE = 0.9
GAP, TOP, BOT, PADT = 30, 44, 34, 6
NAVY = (20, 20, 60)
BG = (22, 22, 22)
FLAT = ["PHerc0800_z15264_w100", "PHerc0826_auto_grown_20260902140900360"]
SWIRL = ["PHerc0826_v3_s0", "PHerc0358_v3_s7", "PHerc0257_v2_s0", "PHerc1203_v3_s0"]
HELD = ["a144", "a174", "p500"]


def kernel(sigma):
    r = int(np.ceil(3 * sigma))
    x = np.arange(-r, r + 1, dtype=np.float32)
    k = np.exp(-0.5 * (x / sigma) ** 2)
    return (k / k.sum()).astype(np.float32)


KER = kernel(SIGMA)


def conv_sep(a, k):
    r = len(k) // 2
    ap = np.pad(a, ((r, r), (0, 0)), mode="edge")
    out = np.zeros_like(a)
    for i, kv in enumerate(k):
        out += kv * ap[i:i + a.shape[0]]
    ap = np.pad(out, ((0, 0), (r, r)), mode="edge")
    out2 = np.zeros_like(a)
    for i, kv in enumerate(k):
        out2 += kv * ap[:, i:i + a.shape[1]]
    return out2


def erode(v, r):
    out = v.copy()
    for dy in range(-r, r + 1):
        for dx in range(-r, r + 1):
            if dy == 0 and dx == 0:
                continue
            sh = np.zeros_like(v)
            ys = slice(max(dy, 0), v.shape[0] + min(dy, 0))
            yd = slice(max(-dy, 0), v.shape[0] + min(-dy, 0))
            xs = slice(max(dx, 0), v.shape[1] + min(dx, 0))
            xd = slice(max(-dx, 0), v.shape[1] + min(-dx, 0))
            sh[yd, xd] = v[ys, xs]
            out &= sh
    return out


def score_row(st, valid, r0, nx, rng=None, null=None):
    """Tiles of tile-row starting at canvas row r0: returns dict of arrays over nx tiles."""
    H, W = valid.shape
    ra, rb = max(r0 - PAD, 0), min(r0 + T + PAD, H)
    v = valid[ra:rb, :nx * T]
    p = np.asarray(st[PLANE, ra:rb, :nx * T]).astype(np.float32)
    if null == "spatial":
        # permute plane-14 pixels inside each tile (padding rows kept as is: scored px are the tile interior)
        off = r0 - ra
        for j in range(nx):
            blk = p[off:off + T, j * T:(j + 1) * T]
            fl = blk.reshape(-1).copy()
            rng.shuffle(fl)
            p[off:off + T, j * T:(j + 1) * T] = fl.reshape(T, T)
    vf = v.astype(np.float32)
    num = conv_sep(p * vf, KER)
    den = conv_sep(vf, KER)
    sm = np.where(den > 1e-3, num / np.maximum(den, 1e-3), 0.0)
    gy = np.zeros_like(sm)
    gx = np.zeros_like(sm)
    gy[1:-1] = 0.5 * (sm[2:] - sm[:-2])
    gx[:, 1:-1] = 0.5 * (sm[:, 2:] - sm[:, :-2])
    ok = erode(v, PAD - 1)
    z = (gx + 1j * gy)
    w = (gx * gx + gy * gy)
    ok &= w > 1e-8
    zn = np.where(ok, z, 0)
    z4 = zn ** 4 / np.maximum(w, 1e-8)
    off = r0 - ra
    h = min(T, H - r0)
    A = np.full(nx, np.nan, np.float32)
    vfrac = np.zeros(nx, np.float32)
    for j in range(nx):
        sl = (slice(off, off + h), slice(j * T, (j + 1) * T))
        sw = float(np.where(ok[sl], w[sl], 0).sum())
        if sw > 0:
            A[j] = abs(complex(z4[sl].sum())) / sw
        vfrac[j] = float(valid[r0:r0 + h, j * T:(j + 1) * T].mean())
    # B: block means of valid px for each plane
    nb = T // BLK
    vb = valid[r0:r0 + T, :nx * T].reshape(T // BLK, BLK, nx * nb, BLK)
    cnt = vb.sum((1, 3)).astype(np.float32)                            # (nb, nx*nb)
    bm = np.zeros((st.shape[0], nb, nx * nb), np.float32)
    for k in range(st.shape[0]):
        pk = np.asarray(st[k, r0:r0 + T, :nx * T]).astype(np.float32) * valid[r0:r0 + T, :nx * T]
        bm[k] = pk.reshape(nb, BLK, nx * nb, BLK).sum((1, 3)) / np.maximum(cnt, 1)
    good = cnt >= 0.5 * BLK * BLK
    if null == "depth":
        bm = rng.permuted(bm, axis=0)
    pad = np.concatenate([bm[:1], bm, bm[-1:]], 0)
    sm_d = 0.25 * pad[:-2] + 0.5 * pad[1:-1] + 0.25 * pad[2:]
    d = sm_d.argmax(0).astype(np.float32)                              # (nb, nx*nb)
    B = np.full(nx, np.nan, np.float32)
    med = np.full(nx, np.nan, np.float32)
    for j in range(nx):
        g = good[:, j * nb:(j + 1) * nb]
        dd = d[:, j * nb:(j + 1) * nb][g]
        if dd.size >= 0.5 * nb * nb:
            m = float(np.median(dd))
            B[j] = float((np.abs(dd - m) <= 2).mean())
            med[j] = m
    return {"A": A, "B": B, "med_d": med, "vfrac": vfrac}


def score_stack(st, valid, rng=None, null=None):
    H, W = valid.shape
    ny, nx = H // T, W // T
    out = {k: np.full((ny, nx), np.nan, np.float32) for k in ("A", "B", "med_d", "vfrac")}
    for i in range(ny):
        r = score_row(st, valid, i * T, nx, rng, null)
        for k in out:
            out[k][i] = r[k]
    return out


def load(name, kind):
    if kind == "factory":
        st = np.load(FAC / "stacks" / f"{name}.npy", mmap_mode="r")
        valid = np.load(FAC / "stacks" / f"{name}.valid.npy")
    else:
        st = np.load(REFDIR / f"{name}.npy", mmap_mode="r")
        valid = np.load(REFDIR / f"{name}.valid.npy")
    assert st.shape[0] == 28 and st.shape[1:] == valid.shape, (name, st.shape, valid.shape)
    return st, valid.astype(bool)


def vt(o):
    return np.isfinite(o["A"]) & np.isfinite(o["B"]) & (o["vfrac"] >= VALID_TILE)


def status(msg, **kw):
    OUT.mkdir(parents=True, exist_ok=True)
    d = {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "msg": msg}
    d.update(kw)
    (OUT / "status.json").write_text(json.dumps(d, indent=1))
    print(msg, flush=True)


def flat_summary(name, o, tA, tB, um):
    m = vt(o)
    on = m & (o["A"] >= tA) & (o["B"] >= tB)
    nv = int(m.sum())
    return {"surface": name, "um_per_px": um, "tiles_valid": nv, "tiles_on_sheet": int(on.sum()), "share": (float(on.sum() / nv) if nv else None),
            "area_cm2": float(on.sum() * (T * um / 1e4) ** 2), "valid_area_cm2": float(nv * (T * um / 1e4) ** 2),
            "A_median": float(np.nanmedian(o["A"][m])) if nv else None, "B_median": float(np.nanmedian(o["B"][m])) if nv else None,
            "share_A_only": float((m & (o["A"] >= tA))[m].mean()) if nv else None, "share_B_only": float((m & (o["B"] >= tB))[m].mean()) if nv else None}


def auc(pos, neg):
    s = [(1.0 if p > n else 0.5 if p == n else 0.0) for p in pos for n in neg]
    return float(np.mean(s)) if s else None


def tile_um(name):
    z = np.load(FAC / "maps" / f"{name}.npz")
    return float(z["um_per_px"])


def pictures(rank_rows, frows, on_by_name):
    """[scan plane 14 | writing-order map] at 40 um/px, non-on-sheet tiles and off-surface painted navy; panels cut from the stored X-198 PNG."""
    fr = {r["surface"]: r for r in frows}
    pdir = OUT / "top20"
    pdir.mkdir(parents=True, exist_ok=True)
    made = []
    for rank, row in enumerate(rank_rows[:20], 1):
        name = row["surface"]
        f = fr[name]
        w, h = f["panel_px"]
        um = f["um_per_canvas_px"]
        scale = um / pic245.TARGET_UM
        z = np.load(FAC / "maps" / f"{name}.npz")
        valid = z["valid"]
        H, W = valid.shape
        oh, ow = max(1, int(round(H * scale))), max(1, int(round(W * scale)))
        Fs = pic245._resample(np.where(valid, pic245.maps_to_min(z, "F"), np.nan), valid, scale)
        Rs = pic245._resample(np.where(valid, pic245.maps_to_min(z, "R"), np.nan), valid, scale)
        m = np.isfinite(Fs) | np.isfinite(Rs)
        ys, xs = np.where(m)
        y0, x0 = ys.min(), xs.min()
        del Fs, Rs, m
        assert ys.max() + 1 - y0 == h and xs.max() + 1 - x0 == w, (name, h, w)
        flag = on_by_name[name]
        ycan = (np.arange(y0, y0 + h) + 0.5) * (H / oh) - 0.5
        xcan = (np.arange(x0, x0 + w) + 0.5) * (W / ow) - 0.5
        ti = np.clip(np.floor(ycan / T).astype(int), 0, None)
        tj = np.clip(np.floor(xcan / T).astype(int), 0, None)
        inside = (ti < flag.shape[0])[:, None] & (tj < flag.shape[1])[None, :]
        fm = np.zeros((h, w), bool)
        fm[inside] = flag[np.minimum(ti, flag.shape[0] - 1)[:, None], np.minimum(tj, flag.shape[1] - 1)[None, :]][inside]
        src = Image.open(FAC / "framed" / f"{name}.png").convert("RGB")
        p1 = np.array(src.crop((0, TOP, w, TOP + h)))
        p3 = np.array(src.crop((2 * w + 2 * GAP, TOP, 3 * w + 2 * GAP, TOP + h)))
        p1[~fm] = NAVY
        p3[~fm] = NAVY
        Wimg = max(2 * w + GAP, pic245.BAR_PX + 520)
        img = Image.new("RGB", (Wimg, TOP + h + BOT), BG)
        img.paste(Image.fromarray(p1), (0, TOP))
        img.paste(Image.fromarray(p3), (w + GAP, TOP))
        d = ImageDraw.Draw(img)
        d.text((PADT, 3), f"#{rank} {name}  on-sheet {row['area_cm2']:.2f} cm2 = {100 * row['share']:.0f}% of valid  panel {w * 0.04:.1f} x {h * 0.04:.1f} mm  40 um per pixel",
               fill=(235, 235, 235), font=pic245.font(14))
        d.text((PADT, 22), "scan plane 14, on-sheet tiles only (navy = not on-sheet / off-surface)", fill=(235, 235, 235), font=pic245.font(12))
        d.text((w + GAP + PADT, 22), "writing order (X-198): F where outward, R where inward, grey undetermined; on-sheet tiles only", fill=(235, 235, 235), font=pic245.font(12))
        yb = TOP + h + 14
        d.rectangle([PADT, yb, PADT + pic245.BAR_PX - 1, yb + 3], fill=(255, 255, 255))
        d.text((PADT + pic245.BAR_PX + 8, yb - 6), "1 cm  (40 um per pixel)", fill=(235, 235, 235), font=pic245.font(14))
        p = pdir / f"{rank:02d}_{name}.png"
        img.save(p)
        made.append({"rank": rank, "surface": name, "png": str(p), "png_px": list(img.size)})
        status(f"picture {rank} {name}")
    return made


def overview(rows, flags, path, cell=3):
    names = [r["surface"] for r in rows]
    ncol = 10
    mx = max(f[0].shape[1] for f in flags.values())
    my = max(f[0].shape[0] for f in flags.values())
    cw, ch = mx * cell + 8, my * cell + 22
    nrow = (len(names) + ncol - 1) // ncol
    img = Image.new("RGB", (ncol * cw, nrow * ch), BG)
    d = ImageDraw.Draw(img)
    for k, n in enumerate(names):
        fl, vm = flags[n]
        ox, oy = (k % ncol) * cw + 4, (k // ncol) * ch + 18
        a = np.full(fl.shape + (3,), NAVY, np.uint8)
        a[vm & ~fl] = (150, 40, 40)
        a[fl] = (230, 230, 230)
        a = np.kron(a, np.ones((cell, cell, 1), np.uint8))
        img.paste(Image.fromarray(a), (ox, oy))
        d.text((ox, oy - 14), f"{k + 1} {n[:30]}", fill=(235, 235, 235), font=pic245.font(10))
    img.save(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--skip-pictures", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    TILES.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    # (1) thresholds from w00
    thp = OUT / "thresholds.json"
    ref = {}
    if thp.exists():
        th = json.loads(thp.read_text())
        tA, tB = th["tA"], th["tB"]
        status("thresholds loaded (frozen)", tA=tA, tB=tB)
    else:
        status("scoring w00")
        st, valid = load("w00", "ref")
        o = score_stack(st, valid)
        m = vt(o)
        nfin = int(np.isfinite(o["A"][o["vfrac"] >= VALID_TILE]).sum())
        assert int(m.sum()) >= 100, "kill (a): fewer than 100 valid tiles"
        tA, tB = float(np.percentile(o["A"][m], 10)), float(np.percentile(o["B"][m], 10))
        th = {"tA": tA, "tB": tB, "w00_valid_tiles": int(m.sum()), "w00_tiles_nonfinite": int((o["vfrac"] >= VALID_TILE).sum() - m.sum()),
              "A_pctiles": {str(q): float(np.percentile(o["A"][m], q)) for q in (1, 5, 10, 25, 50, 75, 90, 99)},
              "B_pctiles": {str(q): float(np.percentile(o["B"][m], q)) for q in (1, 5, 10, 25, 50, 75, 90, 99)},
              "tile_px": T, "block_px": BLK, "plane": PLANE, "sigma": SIGMA, "made": time.strftime("%Y-%m-%d %H:%M:%S")}
        thp.write_text(json.dumps(th, indent=1))
        th["sha256"] = hashlib.sha256(thp.read_bytes()).hexdigest()
        np.savez(TILES / "w00.npz", **o)
        ref["w00"] = flat_summary("w00", o, tA, tB, 9.366)
        (OUT / "thresholds.sha256").write_text(th["sha256"])
        status("thresholds fixed", tA=tA, tB=tB, sha=th["sha256"])
        # (2) nulls
        for null in ("spatial", "depth"):
            on = score_stack(st, valid, np.random.default_rng(SEED), null)
            ref["w00_null_" + null] = flat_summary("w00_null_" + null, on, tA, tB, 9.366)
            ref["w00_null_" + null]["frac_A_below"] = float((on["A"][vt(on)] < tA).mean())
            ref["w00_null_" + null]["frac_B_below"] = float((on["B"][vt(on)] < tB).mean())
            status(f"null {null}: on-sheet share {ref['w00_null_' + null]['share']:.4f}")
        del st, valid
        # (3) held-out known writing
        for n in HELD:
            st, valid = load(n, "ref")
            o = score_stack(st, valid)
            np.savez(TILES / f"{n}.npz", **o)
            ref[n] = flat_summary(n, o, tA, tB, 9.37)
            status(f"held-out {n}: share {ref[n]['share']:.3f}")
            del st, valid
        (OUT / "reference.json").write_text(json.dumps(ref, indent=1))
    # (4) factory
    fidx = json.loads((FAC / "framed/index.json").read_text())["surfaces"]
    rows = []
    t0 = time.time()
    for k, f in enumerate(fidx):
        n = f["surface"]
        if a.only and n not in a.only:
            continue
        tp = TILES / f"{n}.npz"
        um = f["um_per_canvas_px"]
        if tp.exists():
            o = dict(np.load(tp))
        else:
            st, valid = load(n, "factory")
            z = np.load(FAC / "maps" / f"{n}.npz")
            assert tuple(z["canvas"]) == valid.shape, n
            o = score_stack(st, valid)
            np.savez(tp, **o)
            del st, valid
        r = flat_summary(n, o, tA, tB, um)
        r["scroll"] = f["scroll"]
        r["order"] = k
        rows.append(r)
        status(f"{k + 1}/120 {n} share {r['share']} area {r['area_cm2']:.2f} cm2", done=len(rows), seconds=round(time.time() - t0))
    if a.only:
        (OUT / "rows_partial.json").write_text(json.dumps(rows, indent=1))
        return
    # (5) ranking
    prim = sorted(rows, key=lambda r: (-r["area_cm2"], -(r["share"] or 0)))
    sec = sorted(rows, key=lambda r: (-(r["share"] or 0), -r["area_cm2"]))
    for i, r in enumerate(prim):
        r["rank_area"] = i + 1
    for i, r in enumerate(sec):
        r["rank_share"] = i + 1
    ov = len({r["surface"] for r in prim[:20]} & {r["surface"] for r in sec[:20]})
    by = {r["surface"]: r for r in rows}
    v2 = {"flat": {n: by[n]["share"] for n in FLAT if n in by}, "swirl": {n: by[n]["share"] for n in SWIRL if n in by}}
    v2["auc"] = auc(list(v2["flat"].values()), list(v2["swirl"].values()))
    shares = np.array([r["share"] for r in rows if r["share"] is not None])
    summ = {"made": time.strftime("%Y-%m-%d %H:%M:%S"), "experiment": "X-199", "tA": tA, "tB": tB, "n": len(rows), "median_share": float(np.median(shares)),
            "share_quartiles": [float(np.percentile(shares, q)) for q in (10, 25, 50, 75, 90)], "n_share_ge_0.5": int((shares >= 0.5).sum()),
            "total_on_sheet_cm2": float(sum(r["area_cm2"] for r in rows)), "total_valid_cm2": float(sum(r["valid_area_cm2"] for r in rows)),
            "top20_overlap_area_vs_share": ov, "V2": v2, "rows": prim}
    (OUT / "ranking.json").write_text(json.dumps(summ, indent=1))
    md = ["| rank | surface | on-sheet cm2 | valid cm2 | share | A med | B med | rank by share |", "|---|---|---|---|---|---|---|---|"]
    for r in prim:
        md.append(f"| {r['rank_area']} | {r['surface']} | {r['area_cm2']:.2f} | {r['valid_area_cm2']:.2f} | {r['share']:.2f} | {r['A_median']:.2f} | {r['B_median']:.2f} | {r['rank_share']} |")
    (OUT / "ranking.md").write_text("\n".join(md))
    flags = {}
    for r in prim:
        o = np.load(TILES / f"{r['surface']}.npz")
        vm = vt(o)
        flags[r["surface"]] = ((vm & (o["A"] >= tA) & (o["B"] >= tB)), vm)
    overview(prim, flags, OUT / "onsheet_map_all.png")
    status("ranking written", median_share=summ["median_share"], auc=v2["auc"])
    if a.skip_pictures:
        return
    made = pictures(prim, fidx, {n: v[0] for n, v in flags.items()})
    (OUT / "top20.json").write_text(json.dumps(made, indent=1))
    status("done", done=len(rows), pictures=len(made))


if __name__ == "__main__":
    main()
