"""T-268 (X-212) pictures and reading summaries from the stored Arc maps (CPU).  Descriptive only; model ink is not letters.
  picture(name, patch)   data/t268-trace/render/maps/<name>.npz -> outputs/t268-trace/<name>.png : 40 um per pixel EXACTLY (same area-filter resampling as pic245), left = min(Reader v2 F, d9v4C F) normal
                         order, right = min(... R) reversed order, 1 cm bar = 250 px; the FRAME-CORRECT panel (X-198 rule: outward patch -> F, inward patch -> R) carries an orange frame and the words FRAME-CORRECT.
  summary(name)          shares of valid pixels above 0.5 / 0.8, mean, P99 of min(Reader v2, d9v4C) per order, reader-to-reader Pearson r (canvas pixels, as X-211).
  contact(scroll, names, out_png)   frame-correct panels only, each reduced to fit a 420 px tile, with its own reduced scale bar (1 cm) and name, cm2, shares.
"""
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common268 as C  # noqa: E402

sys.path.insert(0, str(C.PROB / "work/t245-picture-factory"))
import pic245  # noqa: E402

ORANGE = (255, 150, 0)
TARGET_UM = pic245.TARGET_UM
BAR_PX = pic245.BAR_PX


def load_maps(name):
    z = np.load(C.RENDER / "maps" / f"{name}.npz")
    valid = z["valid"].astype(bool)
    m = {k: z[k].astype(np.float32) for k in ("rv2__F", "rv2__R", "d9v4C__F", "d9v4C__R")}
    F = np.minimum(m["rv2__F"], m["d9v4C__F"])
    R = np.minimum(m["rv2__R"], m["d9v4C__R"])
    return m, F, R, valid


def summary(name):
    m, F, R, valid = load_maps(name)
    out = {}
    for tag, mm, a, b in (("F", F, m["rv2__F"], m["d9v4C__F"]), ("R", R, m["rv2__R"], m["d9v4C__R"])):
        v = mm[valid]
        out[tag] = {"share_gt_0.5": float((v > 0.5).mean()), "share_gt_0.8": float((v > 0.8).mean()), "mean": float(v.mean()), "p99": float(np.percentile(v, 99)),
                    "reader_r": float(np.corrcoef(a[valid], b[valid])[0, 1])}
    out["valid_px"] = int(valid.sum())
    return out


def panels(name, um_per_px):
    m, F, R, valid = load_maps(name)
    f = um_per_px / TARGET_UM
    Fs = pic245._resample(np.where(valid, F, np.nan), valid, f)
    Rs = pic245._resample(np.where(valid, R, np.nan), valid, f)
    mk = np.isfinite(Fs) | np.isfinite(Rs)
    ys, xs = np.where(mk)
    y0, y1, x0, x1 = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1
    g = lambda a: (np.clip(np.nan_to_num(a[y0:y1, x0:x1], nan=0.0), 0, 1) * 255).astype(np.uint8)
    return g(Fs), g(Rs)


def picture(name, patch, out_png):
    """patch: dict with side (+1 outward, -1 inward), um_per_canvas_px, scroll, voxel_um, cm2"""
    pf, pr = panels(name, patch["um_per_canvas_px"])
    h, w = pf.shape
    fc = "F" if patch["side"] == 1 else "R"
    size_mm = (round(w * TARGET_UM / 1000, 1), round(h * TARGET_UM / 1000, 1))
    gap, top, bot, pad = 30, 44, 34, 6
    Wimg = max(2 * w + gap, BAR_PX + 120, 700)
    img = Image.new("RGB", (Wimg, top + h + bot), (22, 22, 22))
    d = ImageDraw.Draw(img)
    img.paste(Image.fromarray(pf).convert("RGB"), (0, top))
    img.paste(Image.fromarray(pr).convert("RGB"), (w + gap, top))
    d.text((pad, 4), f"{name}  ({size_mm[0]:g} x {size_mm[1]:g} mm, {patch['cm2']:.2f} cm2)  {patch['scroll']}  {patch['voxel_um']} um scan", fill=(235, 235, 235), font=pic245.font(14))
    d.text((pad, 24), "left: normal order (F)", fill=(200, 200, 200), font=pic245.font(14))
    d.text((w + gap + pad, 24), "right: reversed order (R)", fill=(200, 200, 200), font=pic245.font(14))
    x0 = 0 if fc == "F" else w + gap
    d.rectangle([x0 - 2, top - 2, x0 + w + 1, top + h + 1], outline=ORANGE, width=3)
    tag = "FRAME-CORRECT"
    d.rectangle([x0 + 3, top + 3, x0 + 3 + 8 * len(tag) + 8, top + 3 + 18], fill=(22, 22, 22))
    d.text((x0 + 7, top + 5), tag, fill=ORANGE, font=pic245.font(14))
    yb = top + h + 14
    d.rectangle([pad, yb, pad + BAR_PX - 1, yb + 3], fill=(255, 255, 255))
    d.text((pad + BAR_PX + 8, yb - 6), "1 cm  (40 um per pixel)", fill=(235, 235, 235), font=pic245.font(14))
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    img.save(out_png)
    return {"png": str(out_png), "size_mm": list(size_mm), "panel_px": [int(w), int(h)], "image_px": list(img.size), "frame_correct": fc}


def contact(items, out_png, title, tile=420, cols=4, which="fc"):
    """items: list of dict(name, patch, label). Frame-correct panel only (40 um/px reduced by a per-tile factor stated on the tile)."""
    tiles = []
    for it in items:
        pf, pr = panels(it["name"], it["patch"]["um_per_canvas_px"])
        pan = pf if (it["patch"]["side"] == 1) == (which == "fc") else pr
        h, w = pan.shape
        k = min(1.0, (tile - 8) / max(h, w))
        im = Image.fromarray(pan)
        if k < 1.0:
            im = im.resize((max(1, int(round(w * k))), max(1, int(round(h * k)))), Image.BOX)
        tiles.append((im, k, it))
    rows = (len(tiles) + cols - 1) // cols
    cell_h = tile + 62
    img = Image.new("RGB", (cols * (tile + 8) + 8, 36 + rows * cell_h), (22, 22, 22))
    d = ImageDraw.Draw(img)
    d.text((8, 8), title + ("   (frame-correct order only; " if which == "fc" else "   (the OTHER order only, not frame-correct; ") + "each tile reduced by the factor shown, with its own 1 cm bar)", fill=(235, 235, 235), font=pic245.font(14))
    for i, (im, k, it) in enumerate(tiles):
        cx, cy = 8 + (i % cols) * (tile + 8), 36 + (i // cols) * cell_h
        d.text((cx, cy), it["label"], fill=(235, 235, 235), font=pic245.font(12))
        d.text((cx, cy + 14), it["label2"], fill=ORANGE, font=pic245.font(12))
        img.paste(im.convert("RGB"), (cx, cy + 30))
        yb = cy + 30 + im.size[1] + 8
        bar = max(2, int(round(BAR_PX * k)))
        d.rectangle([cx, yb, cx + bar - 1, yb + 3], fill=(255, 255, 255))
        d.text((cx + bar + 6, yb - 5), f"1 cm   (x{k:.2f}, {TARGET_UM / k:.0f} um/px)", fill=(220, 220, 220), font=pic245.font(11))
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    img.save(out_png)
    return str(out_png)
