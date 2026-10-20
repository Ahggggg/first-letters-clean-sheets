"""T-259 / X-202: X-201 on-sheet score with a two-sided depth-contrast band and a texture floor, all thresholds from known writing w00 only.
Reuses X-201's stored per-tile A, C(label 'B'), vfrac (data/picture-factory/onsheet_c/tiles) for w00, a144, a174, p500 and the 120 surfaces; computes only the new
texture statistic Tx per tile (plane 14, one plane read) and, for the w00 nulls and a reproduction check, re-runs the X-201 scorer (onsheet257c, imported unchanged).
  Tx = median over the tile's interior valid px of |grad| of the sigma-1.5 smoothed plane 14, divided by the mean of plane 14 over the tile's valid px (tile >= half valid px interior).
  on-sheet iff A >= tA and tC <= C <= cHi and Tx >= tT; tA, tC from X-201 thresholds.json (w00 p10), cHi = w00 p99 of C, tT = w00 p10 of Tx.
Order: (1) w00 -> thresholds259.json (frozen, sha256) (2) reproduction check and nulls on w00 (3) held-out a144 a174 p500 (4) 120 factory surfaces (5) ranking, X-201 comparison, top-20 pictures.
Seed 20261259 for the Tx spatial-null shuffle; X-201 scorer nulls use its own seed 20261257.
"""
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROB = HERE.parents[1]
sys.path.insert(0, str(PROB / "work/t257-onsheet"))
import onsheet257c as c  # noqa: E402  (patches onsheet257.score_row; imported unchanged)
import onsheet257 as m  # noqa: E402

FAC = m.FAC
SRC = FAC / "onsheet_c"
OUT = FAC / "onsheet_d"
TILES = OUT / "tiles"
m.OUT = OUT
m.TILES = TILES
SEED = 20261259
T = m.T
PLANE = m.PLANE


def status(msg, **kw):
    OUT.mkdir(parents=True, exist_ok=True)
    d = {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "msg": msg}
    d.update(kw)
    (OUT / "status.json").write_text(json.dumps(d, indent=1))
    print(msg, flush=True)


def texture(st, valid, rng=None):
    """Tx per 256 px tile of the canvas grid (H//T x W//T)."""
    H, W = valid.shape
    ny, nx = H // T, W // T
    p = np.asarray(st[PLANE]).astype(np.float32)
    if rng is not None:
        for i in range(ny):
            for j in range(nx):
                blk = p[i * T:(i + 1) * T, j * T:(j + 1) * T]
                fl = blk.reshape(-1).copy()
                rng.shuffle(fl)
                p[i * T:(i + 1) * T, j * T:(j + 1) * T] = fl.reshape(T, T)
    vf = valid.astype(np.float32)
    num = m.conv_sep(p * vf, m.KER)
    den = m.conv_sep(vf, m.KER)
    sm = np.where(den > 1e-3, num / np.maximum(den, 1e-3), 0.0)
    gy = np.zeros_like(sm)
    gx = np.zeros_like(sm)
    gy[1:-1] = 0.5 * (sm[2:] - sm[:-2])
    gx[:, 1:-1] = 0.5 * (sm[:, 2:] - sm[:, :-2])
    g = np.sqrt(gx * gx + gy * gy)
    ok = m.erode(valid, m.PAD - 1)
    out = np.full((ny, nx), np.nan, np.float32)
    for i in range(ny):
        for j in range(nx):
            sl = (slice(i * T, (i + 1) * T), slice(j * T, (j + 1) * T))
            o = ok[sl]
            v = valid[sl]
            if o.sum() >= 0.5 * T * T and v.sum() > 0:
                mu = float(p[sl][v].mean())
                if mu > 0:
                    out[i, j] = float(np.median(g[sl][o])) / mu
    return out


def vt(o):
    return np.isfinite(o["A"]) & np.isfinite(o["B"]) & (o["vfrac"] >= m.VALID_TILE)


def onmask(o, th):
    ok = vt(o) & np.isfinite(o["Tx"])
    return ok & (o["A"] >= th["tA"]) & (o["B"] >= th["tC"]) & (o["B"] <= th["cHi"]) & (o["Tx"] >= th["tT"])


def summ(name, o, th, um):
    v = vt(o) & np.isfinite(o["Tx"])
    on = onmask(o, th)
    nv = int(v.sum())
    a = lambda x: float(x[v].mean()) if nv else None
    return {"surface": name, "um_per_px": um, "tiles_valid": nv, "tiles_on_sheet": int(on.sum()), "share": float(on.sum() / nv) if nv else None,
            "area_cm2": float(on.sum() * (T * um / 1e4) ** 2), "valid_area_cm2": float(nv * (T * um / 1e4) ** 2),
            "share_x201": float((v & (o["A"] >= th["tA"]) & (o["B"] >= th["tC"]))[v].mean()) if nv else None,
            "share_A": a(o["A"] >= th["tA"]), "share_Clow": a(o["B"] >= th["tC"]), "share_Chigh_dropped": a(o["B"] > th["cHi"]), "share_Tx": a(o["Tx"] >= th["tT"]),
            "Tx_median": float(np.nanmedian(o["Tx"][v])) if nv else None, "C_median": float(np.nanmedian(o["B"][v])) if nv else None}


def auc(pos, neg):
    s = [(1.0 if p > n else 0.5 if p == n else 0.0) for p in pos for n in neg]
    return float(np.mean(s)) if s else None


def main():
    skip_pics = "--skip-pictures" in sys.argv
    OUT.mkdir(parents=True, exist_ok=True)
    TILES.mkdir(parents=True, exist_ok=True)
    x201 = json.loads((SRC / "thresholds.json").read_text())
    tA, tC, cHi = x201["tA"], x201["tB"], x201["B_pctiles"]["99"]
    ref = {}
    thp = OUT / "thresholds259.json"
    st, valid = m.load("w00", "ref")
    if thp.exists():
        th = json.loads(thp.read_text())
        status("thresholds loaded (frozen)")
    else:
        status("w00: reproduction check and Tx")
        o = dict(np.load(SRC / "tiles/w00.npz"))
        rec = m.score_stack(st, valid)  # patched X-201 scorer
        diff = {k: float(np.nanmax(np.abs(rec[k] - o[k]))) for k in ("A", "B")}
        same_nan = bool((np.isnan(rec["A"]) == np.isnan(o["A"])).all() and (np.isnan(rec["B"]) == np.isnan(o["B"])).all())
        tx = texture(st, valid)
        o["Tx"] = tx
        mm = vt(o) & np.isfinite(tx)
        assert int(mm.sum()) >= 100, "kill: fewer than 100 valid w00 tiles"
        tT = float(np.percentile(tx[mm], 10))
        th = {"tA": tA, "tC": tC, "cHi": cHi, "tT": tT, "w00_valid_tiles": int(mm.sum()), "reproduction_max_abs_diff": diff, "reproduction_nan_pattern_same": same_nan,
              "Tx_pctiles": {str(q): float(np.percentile(tx[mm], q)) for q in (1, 5, 10, 25, 50, 75, 90, 99)},
              "source_thresholds_sha256": hashlib.sha256((SRC / "thresholds.json").read_bytes()).hexdigest(), "made": time.strftime("%Y-%m-%d %H:%M:%S")}
        thp.write_text(json.dumps(th, indent=1))
        (OUT / "thresholds259.sha256").write_text(hashlib.sha256(thp.read_bytes()).hexdigest())
        np.savez(TILES / "w00.npz", **o)
        status("thresholds fixed", **{k: th[k] for k in ("tA", "tC", "cHi", "tT")})
        ref["w00"] = summ("w00", o, th, 9.366)
        # nulls
        for null in ("spatial", "depth"):
            on = m.score_stack(st, valid, np.random.default_rng(m.SEED), null)
            if null == "spatial":
                on["Tx"] = texture(st, valid, np.random.default_rng(SEED))
            else:
                on["Tx"] = tx  # depth shuffle leaves plane 14 order of planes only; Tx unaffected by construction (plane 14 is the same image)
            r = summ("w00_null_" + null, on, th, 9.366)
            v = vt(on)
            r["frac_A_below"] = float((on["A"][v] < tA).mean())
            r["frac_C_below_tC"] = float((on["B"][v] < tC).mean())
            r["frac_C_above_cHi"] = float((on["B"][v] > cHi).mean())
            ref["w00_null_" + null] = r
            status(f"null {null}: combined on-sheet share {r['share']:.4f}")
        del st, valid
        for n in m.HELD:
            o = dict(np.load(SRC / f"tiles/{n}.npz"))
            s2, v2 = m.load(n, "ref")
            o["Tx"] = texture(s2, v2)
            np.savez(TILES / f"{n}.npz", **o)
            ref[n] = summ(n, o, th, 9.37)
            status(f"held-out {n}: share {ref[n]['share']:.3f} (X-201 rule {ref[n]['share_x201']:.3f})")
            del s2, v2
        (OUT / "reference.json").write_text(json.dumps(ref, indent=1))
    # factory
    fidx = json.loads((FAC / "framed/index.json").read_text())["surfaces"]
    rows = []
    t0 = time.time()
    for k, f in enumerate(fidx):
        n = f["surface"]
        tp = TILES / f"{n}.npz"
        if tp.exists():
            o = dict(np.load(tp))
        else:
            o = dict(np.load(SRC / f"tiles/{n}.npz"))
            s2, v2 = m.load(n, "factory")
            assert o["A"].shape == (v2.shape[0] // T, v2.shape[1] // T), n
            o["Tx"] = texture(s2, v2)
            np.savez(tp, **o)
            del s2, v2
        r = summ(n, o, th, f["um_per_canvas_px"])
        r["scroll"] = f["scroll"]
        r["order"] = k
        rows.append(r)
        status(f"{k + 1}/120 {n} share {r['share']} area {r['area_cm2']:.2f} cm2", done=len(rows), seconds=round(time.time() - t0))
    x201r = {r["surface"]: r for r in json.loads((SRC / "ranking.json").read_text())["rows"]}
    prim = sorted(rows, key=lambda r: (-r["area_cm2"], -(r["share"] or 0)))
    sec = sorted(rows, key=lambda r: (-(r["share"] or 0), -r["area_cm2"]))
    for i, r in enumerate(prim):
        r["rank_area"] = i + 1
        r["rank_area_x201"] = x201r[r["surface"]]["rank_area"]
    for i, r in enumerate(sec):
        r["rank_share"] = i + 1
    by = {r["surface"]: r for r in rows}
    v2 = {"flat": {n: by[n]["share"] for n in m.FLAT}, "swirl": {n: by[n]["share"] for n in m.SWIRL}}
    v2["auc_all6"] = auc(list(v2["flat"].values()), list(v2["swirl"].values()))
    # trio never looked at before X-201's diagnostic (the 4 diagnostic surfaces are 0826_..900360, 0257_v2_s0, 1203_v3_s0 and w00): 0800, 0826_v3_s0, 0358_v3_s7
    v2["auc_unseen_pair"] = auc([by["PHerc0800_z15264_w100"]["share"]], [by["PHerc0826_v3_s0"]["share"], by["PHerc0358_v3_s7"]["share"]])
    shares = np.array([r["share"] for r in rows if r["share"] is not None])
    top20 = {r["surface"] for r in prim[:20]}
    top20_x = {r["surface"] for r in sorted(x201r.values(), key=lambda r: r["rank_area"])[:20]}
    sm = {"made": time.strftime("%Y-%m-%d %H:%M:%S"), "experiment": "X-202", "thresholds": {k: th[k] for k in ("tA", "tC", "cHi", "tT")}, "n": len(rows),
          "median_share": float(np.median(shares)), "share_quartiles": [float(np.percentile(shares, q)) for q in (10, 25, 50, 75, 90)],
          "n_share_ge_0.5": int((shares >= 0.5).sum()), "total_on_sheet_cm2": float(sum(r["area_cm2"] for r in rows)),
          "total_valid_cm2": float(sum(r["valid_area_cm2"] for r in rows)), "top20_overlap_with_x201": len(top20 & top20_x), "top20_overlap_area_vs_share": len(top20 & {r["surface"] for r in sec[:20]}),
          "x201_top5_new_ranks": {n: by[n]["rank_area"] for n in ("PHerc0211_v3_s9", "PHerc0257_v3_s9", "PHerc0826_v3_s9", "PHerc0175A_v2_s0", "PHerc0125_v3_s7")},
          "V2": v2, "rows": prim}
    (OUT / "ranking.json").write_text(json.dumps(sm, indent=1))
    md = ["| rank | surface | on-sheet cm2 | valid cm2 | share | share X-201 rule | C med | Tx med | X-201 rank | rank by share |", "|---|---|---|---|---|---|---|---|---|---|"]
    for r in prim:
        md.append(f"| {r['rank_area']} | {r['surface']} | {r['area_cm2']:.2f} | {r['valid_area_cm2']:.2f} | {r['share']:.2f} | {r['share_x201']:.2f} | {r['C_median']:.1f} | {r['Tx_median']:.4f} | {r['rank_area_x201']} | {r['rank_share']} |")
    (OUT / "ranking.md").write_text("\n".join(md))
    flags = {}
    for r in prim:
        o = dict(np.load(TILES / f"{r['surface']}.npz"))
        flags[r["surface"]] = (onmask(o, th), vt(o))
    m.overview(prim, flags, OUT / "onsheet_map_all.png")
    status("ranking written", median_share=sm["median_share"], auc=v2["auc_all6"])
    if skip_pics:
        return
    made = m.pictures(prim, fidx, {n: v[0] for n, v in flags.items()})
    (OUT / "top20.json").write_text(json.dumps(made, indent=1))
    status("done", done=len(rows), pictures=len(made))


if __name__ == "__main__":
    main()
