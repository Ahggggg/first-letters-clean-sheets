#!/usr/bin/env python3
"""make_qa.py - draw the QA picture of one patch.  MIT licence.  (numpy + Pillow only; the arrays come from xsec.py / the stored renders, see README.)

One PNG per patch, top to bottom:
  A  flattened scan render: the CT sampled on the traced surface (scan plane 14 of the 28-plane render = the surface itself, 9.4 um per pixel; shown at 1/k if wider than MAXW),
     with the two section lines drawn in colour.  Papyrus fibres show as crossed horizontal/vertical texture where the trace follows one sheet; swirled, streaked texture where it crosses layers.
  B, C  two straightened CT cross-sections through the patch, one along each line of A: horizontal = along the line, vertical = distance from the traced surface along its normal
     (-32..+32 voxels = +-0.30 mm, stretched x3 vertically; up = away from the scroll axis).  The traced surface is the middle row (marked by triangles at both ends).  If the trace stays on one
     papyrus layer the layer structure runs horizontally through the middle row; a layer that drifts diagonally across the section means the trace crosses layers.
  D  the frame-correct Reader v2 / d9v4C map: min(Reader v2, d9v4C) in the reading order that matches the scroll's winding direction at this patch (F or R), at 40 um per pixel with a 1 cm bar.
     Model output, not ink and not letters.
"""
import numpy as np
from PIL import Image, ImageDraw, ImageFont

BG = (22, 22, 22)
FG = (235, 235, 235)
DIM = (170, 170, 170)
ORANGE = (255, 150, 0)
CYAN = (60, 200, 255)
NAVY = (30, 30, 120)
MAXW = 1600
VSTRETCH = 3


def font(n):
    try:
        return ImageFont.load_default(size=n)
    except TypeError:
        return ImageFont.load_default()


def stretch(a, valid=None, lo=1, hi=99):
    v = a[valid] if valid is not None and valid.any() else a[np.isfinite(a)]
    p1, p2 = np.percentile(v, [lo, hi])
    return np.clip((a - p1) / max(p2 - p1, 1e-6), 0, 1)


def box_down(a, k):
    """integer block-mean downsample by k along both axes of a 2-D array (partial blocks dropped)."""
    if k == 1:
        return a
    h, w = a.shape[0] // k * k, a.shape[1] // k * k
    return a[:h, :w].reshape(h // k, k, w // k, k).mean((1, 3))


def to_rgb(g01, bad=None):
    g = (np.clip(g01, 0, 1) * 255).astype(np.uint8)
    rgb = np.stack([g, g, g], -1)
    if bad is not None:
        rgb[bad] = NAVY
    return Image.fromarray(rgb)


def section_image(S, valid, d_lo, d_hi, clim, k):
    """S (n_depth, n_s) -> RGB image, vertical stretch x3, horizontal 1/k (block mean), invalid columns navy."""
    g = np.clip((S - clim[0]) / max(clim[1] - clim[0], 1e-6), 0, 1)
    if k > 1:
        n = g.shape[1] // k * k
        gk = g[:, :n].reshape(g.shape[0], -1, k).mean(2)
    else:
        gk = g
    if k > 1:
        n = valid.shape[0] // k * k
        vk = valid[:n].reshape(-1, k).mean(1) > 0.5
    else:
        vk = valid
    img = np.kron(gk, np.ones((VSTRETCH, 1)))
    bad = np.zeros(img.shape, bool)
    bad[:, ~vk[:img.shape[1]]] = True
    return to_rgb(img, bad)


def render_qa(out_png, title, subtitle, scan, scan_valid, lines, sections, ink, depths_vox, um_per_vox, footer):
    """scan: plane-14 canvas (H, W) float; scan_valid bool (H, W); lines: [("row", canvas_row), ("col", canvas_col)] drawn on the scan;
    sections: two dicts {S, valid, label_text}; ink: dict {img (uint8 grey, or RGB), size_mm, order, bar_px}; depths_vox: the section depth axis (voxels)."""
    H, W = scan.shape
    k = max(1, int(np.ceil(W / MAXW)))
    sc = stretch(scan, scan_valid)
    sc = np.where(scan_valid, sc, 0.0)
    img_scan = to_rgb(box_down(sc, k), ~(box_down(scan_valid.astype(float), k) > 0.5))
    d = ImageDraw.Draw(img_scan)
    cols = (ORANGE, CYAN)
    for (kind, idx), c, tag in zip(lines, cols, ("1", "2")):
        if kind == "row":
            y = idx // k
            d.line([(0, y), (img_scan.width, y)], fill=c, width=2)
            d.text((4, max(0, y - 18)), tag, fill=c, font=font(16))
        else:
            x = idx // k
            d.line([(x, 0), (x, img_scan.height)], fill=c, width=2)
            d.text((x + 4, 2), tag, fill=c, font=font(16))
    # section contrast: shared over both sections of the patch
    allv = np.concatenate([s["S"][:, s["valid"]].ravel() for s in sections if s["valid"].any()])
    clim = np.percentile(allv, [1, 99])
    secs = []
    for s in sections:
        ks = max(1, int(np.ceil(s["S"].shape[1] / MAXW)))
        secs.append((section_image(s["S"][::-1], s["valid"], depths_vox[0], depths_vox[-1], clim, ks), ks, s))       # flipped: top row = +32 = outward
    ink_img = ink["img"]
    margin, lab = 62, 22
    panel_w = max(img_scan.width, max(im.width for im, _, _ in secs), ink_img.width + 20, 560)
    Wimg = margin + panel_w + 14
    # heights
    h_head = 52
    h_scan = lab + img_scan.height + 30
    h_secs = sum(lab + im.height + 24 for im, _, _ in secs)
    h_ink = lab + ink_img.height + 44
    h_foot = 46
    canvas = Image.new("RGB", (Wimg, h_head + h_scan + h_secs + h_ink + h_foot), BG)
    dr = ImageDraw.Draw(canvas)
    dr.text((8, 6), title, fill=FG, font=font(16))
    dr.text((8, 28), subtitle, fill=DIM, font=font(13))
    y = h_head
    dr.text((8, y), f"A  flattened scan render (CT on the traced surface, plane 14, {um_per_vox * 1:.2f} um per voxel; shown 1/{k}), section lines drawn", fill=FG, font=font(14))
    y += lab
    canvas.paste(img_scan, (margin, y))
    yb = y + img_scan.height + 12
    bar = int(round(5000 / (um_per_vox * 1.0) / k))
    dr.rectangle([margin, yb, margin + bar - 1, yb + 3], fill=(255, 255, 255))
    dr.text((margin + bar + 8, yb - 6), "5 mm", fill=FG, font=font(13))
    y = yb + 18
    for (im, ks, s), tag, c in zip(secs, ("B", "C"), cols):
        dr.text((8, y), f"{tag}  CT cross-section {'1 (orange)' if tag == 'B' else '2 (cyan)'}: {s['label_text']}   [horizontal 1/{ks}, vertical x{VSTRETCH}; up = away from scroll axis]", fill=c, font=font(14))
        y += lab
        canvas.paste(im, (margin, y))
        mid = y + im.height // 2
        for xx, sgn in ((margin - 10, 1), (margin + im.width + 2, -1)):
            dr.polygon([(xx, mid - 5), (xx, mid + 5), (xx + 8 * sgn, mid)], fill=(255, 60, 60))
        for dv in (-32, -16, 0, 16, 32):
            yy = y + int((depths_vox[-1] - dv + 0.5) * VSTRETCH)
            dr.text((margin - 40, yy - 7), f"{dv:+d}" if dv else "0", fill=DIM, font=font(11))
        dr.text((2, y + 2), "vox", fill=DIM, font=font(10))
        y += im.height + 24
    dr.text((8, y), f"D  frame-correct reading order {ink['order']}: min(Reader v2, d9v4C), grey = model output 0..1 (not ink, not letters); {ink['size_mm'][0]:g} x {ink['size_mm'][1]:g} mm at 40 um per pixel", fill=FG, font=font(14))
    y += lab
    canvas.paste(ink_img if ink_img.mode == "RGB" else ink_img.convert("RGB"), (margin, y))
    yb = y + ink_img.height + 14
    dr.rectangle([margin, yb, margin + ink["bar_px"] - 1, yb + 3], fill=(255, 255, 255))
    dr.text((margin + ink["bar_px"] + 8, yb - 6), "1 cm  (40 um per pixel)", fill=FG, font=font(13))
    y = yb + 22
    dr.text((8, y), footer, fill=DIM, font=font(12))
    canvas.save(out_png, optimize=True)
    return canvas.size


def ink_panel(F, R, valid, side_order, um_per_canvas_px, nonfinite=None, target_um=40.0):
    """frame-correct ink panel: F / R are min(Reader v2, d9v4C) maps (float, NaN allowed) on the canvas, valid bool; returns dict(img RGB, size_mm, order, bar_px)."""
    m = F if side_order == "F" else R
    f = um_per_canvas_px / target_um
    H, W = m.shape
    oh, ow = max(1, int(round(H * f))), max(1, int(round(W * f)))
    fin = np.isfinite(m) & valid
    v0 = np.where(fin, m, 0).astype(np.float32)
    a = np.asarray(Image.fromarray(v0, "F").resize((ow, oh), Image.BOX))
    b = np.asarray(Image.fromarray(fin.astype(np.float32), "F").resize((ow, oh), Image.BOX))
    vv = np.asarray(Image.fromarray(valid.astype(np.float32), "F").resize((ow, oh), Image.BOX))
    with np.errstate(invalid="ignore", divide="ignore"):
        val = np.where(b > 0.5, a / np.maximum(b, 1e-6), 0.0)
    ys, xs = np.where(vv > 0.5)
    y0, y1, x0, x1 = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1
    g = np.clip(val[y0:y1, x0:x1], 0, 1)
    bad = (vv[y0:y1, x0:x1] > 0.5) & ~(b[y0:y1, x0:x1] > 0.5)
    img = to_rgb(g, bad)
    return {"img": img, "size_mm": (round(img.width * target_um / 1000, 1), round(img.height * target_um / 1000, 1)), "order": side_order, "bar_px": int(round(10000 / target_um))}
