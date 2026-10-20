"""T-245 picture factory: the one-PNG-per-surface picture at EXACTLY 40 um per pixel (coordinator format, reference data/t237-new-candidates/lookfirst/scale_01.png).
Left panel = min(Reader v2 F, d9v4C F) (normal order), right panel = min(Reader v2 R, d9v4C R) (reversed order); grey = probability 0..1 (clip, x255), invalid px black; both panels cropped to the
same valid bounding box; header text = surface name, size in mm of the cropped panel (W x H), voxel size; scale bar 250 px = 1 cm at the bottom. Never shrinks to fit: the canvas is resampled by
(micrometres per canvas px) / 40 with an area (BOX) filter, where micrometres per canvas px = voxel_um x (median 3-D grid-node spacing in voxels x mesh scale) is MEASURED from the mesh.
python pic245.py --test NAME_OF_X185_MAP (format test on a stored X-185 map, no model)"""
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

TARGET_UM = 40.0
BAR_PX = 250          # 1 cm at 40 um / px


def px_voxels(points, scale):
    """Voxels per canvas px = median 3-D distance between adjacent valid grid nodes x scale (nodes 1/scale canvas px apart)."""
    d = []
    for ax in (0, 1):
        a = points[:-1] if ax == 0 else points[:, :-1]
        b = points[1:] if ax == 0 else points[:, 1:]
        e = np.linalg.norm(b - a, axis=-1)
        e = e[np.isfinite(e)]
        if len(e):
            d.append(np.median(e))
    return float(np.mean(d) * scale)


def _resample(val, w, f):
    """val float32 (H, W) with NaN invalid -> (value, weight) at the output size round(H*f) x round(W*f), area filter (weight = share of valid input px)."""
    H, W = val.shape
    oh, ow = max(1, int(round(H * f))), max(1, int(round(W * f)))
    m = np.isfinite(val)
    v0 = np.where(m, val, 0).astype(np.float32)
    a = np.asarray(Image.fromarray(v0, "F").resize((ow, oh), Image.BOX))
    b = np.asarray(Image.fromarray(m.astype(np.float32), "F").resize((ow, oh), Image.BOX))
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(b > 0.5, a / np.maximum(b, 1e-6), np.nan)
    return out.astype(np.float32)


def font(n):
    try:
        return ImageFont.load_default(size=n)
    except TypeError:
        return ImageFont.load_default()


def make_picture(name, F, R, valid, um_per_px, out_png, subtitle=""):
    """F, R: min over the two readers (float32, NaN where not read), valid bool mask (canvas coordinates). um_per_px: canvas px size in um. Returns dict with size_mm, panel px and the scale used."""
    f = um_per_px / TARGET_UM
    F = np.where(valid, F, np.nan)
    R = np.where(valid, R, np.nan)
    Fs, Rs = _resample(F, valid, f), _resample(R, valid, f)
    m = np.isfinite(Fs) | np.isfinite(Rs)
    if not m.any():
        raise ValueError("nothing valid")
    ys, xs = np.where(m)
    y0, y1, x0, x1 = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1
    g = lambda a: (np.clip(np.nan_to_num(a[y0:y1, x0:x1], nan=0.0), 0, 1) * 255).astype(np.uint8)
    pf, pr = g(Fs), g(Rs)
    h, w = pf.shape
    size_mm = (round(w * TARGET_UM / 1000, 1), round(h * TARGET_UM / 1000, 1))
    gap, top, bot, pad = 30, 24, 34, 6
    Wimg = max(2 * w + gap, BAR_PX + 120, 560)
    img = Image.new("L", (Wimg, top + h + bot), 22)
    d = ImageDraw.Draw(img)
    img.paste(Image.fromarray(pf), (0, top))
    img.paste(Image.fromarray(pr), (w + gap, top))
    d.text((pad, 4), f"{name}  ({size_mm[0]:g} x {size_mm[1]:g} mm)  left normal order  right reversed  {subtitle}", fill=235, font=font(14))
    yb = top + h + 14
    d.rectangle([pad, yb, pad + BAR_PX - 1, yb + 3], fill=255)
    d.text((pad + BAR_PX + 8, yb - 6), "1 cm  (40 um per pixel)", fill=235, font=font(14))
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    img.save(out_png)
    return {"png": str(out_png), "size_mm": list(size_mm), "panel_px": [int(w), int(h)], "image_px": list(img.size), "scale_f": f, "um_per_px_native": um_per_px}


def maps_to_min(z, arm):
    a = z["rv2__" + arm].astype(np.float32)
    b = z["d9v4C__" + arm].astype(np.float32)
    return np.minimum(a, b)


if __name__ == "__main__":
    PROB = Path(__file__).resolve().parents[2]
    nm = sys.argv[2]
    grp = sys.argv[3] if len(sys.argv) > 3 else "B27"
    z = np.load(PROB / "data/t237-new-candidates/maps" / grp / f"{nm}.npz")
    um = float(sys.argv[4]) if len(sys.argv) > 4 else 9.362
    print(make_picture(nm, maps_to_min(z, "F"), maps_to_min(z, "R"), z["valid"], um, PROB / "data/picture-factory/test" / f"{nm}.png"))
