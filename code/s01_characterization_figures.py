"""
Characterization figures, one image per figure (Cleaner Chemical Engineering version)
=====================================================================================
  fig02_xrd.png       -- XRD patterns of the six treatment controls vs JCPDS 00-004-0477
  fig03_ftir.png      -- FTIR spectra of the same six samples
  fig04_sem_15kx.png  -- SEM micrograph of the design-centre sample, 15,000x
  fig05_sem_30kx.png  -- SEM micrograph of the design-centre sample, 30,000x

Source data (author-supplied, D:/Thanapon/2569/TiO2/Creat_manuscript/New folder):
  - 496398.png : XRD patterns, samples coded N7-N12
  - 496397.png : FTIR spectra, samples coded N7-N12
  - SEM TiO2/.../SEM/TiO2-No.13-*.bmp : design-centre sample (1.250 g glucose, 0.050 g CuSO4, 150 min)

N7-N12 are the treatment controls K1-K6 (confirmed by the authors, 29 October 2026); each code in the supplied
images is covered and replaced by the treatment name. The images are otherwise unchanged. The band gap and the
diffuse-reflectance spectra are no longer reported, so this script no longer draws them.
"""

import os

import numpy as np
from PIL import Image, ImageDraw, ImageFont

FIGDIR = os.environ.get("FIGDIR", "../figures")
os.makedirs(FIGDIR, exist_ok=True)
SRC = "../data/"

# sample code in the supplied images -> treatment control (K1-K6 in the scripts and result files)
TREATMENT = {7: "Heat only (dark)", 8: "Light only", 9: "Untreated TiO₂", 10: "Glucose only", 11: "CuSO₄ only",
             12: "Glucose + CuSO₄"}
# (code, box covering the old label, top-left of the new label), in pixels of the 3201 x 2400 source images
XRD_LABELS = [(12, (736, 352, 912, 448), (504, 356)), (11, (736, 600, 912, 700), (504, 604)),
              (10, (736, 856, 912, 955), (504, 860)), (9, (736, 1104, 872, 1203), (504, 1108)),
              (8, (736, 1356, 872, 1453), (504, 1360)), (7, (736, 1616, 872, 1717), (504, 1620))]
FTIR_LABELS = [(7, (512, 392, 640, 474), (512, 394)), (8, (512, 620, 640, 700), (512, 698)),
               (9, (512, 928, 640, 1012), (512, 944)), (10, (512, 1296, 684, 1380), (512, 1298)),
               (11, (512, 1472, 680, 1556), (512, 1474)), (12, (512, 1768, 684, 1852), (512, 1770))]


def font(size):
    for name in ("times.ttf", "DejaVuSerif.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    raise OSError("no serif font found")


def relabel(path, labels, max_width):
    im = Image.open(path).convert("RGB")
    draw = ImageDraw.Draw(im)
    for code, box, xy in labels:
        px = np.asarray(im.crop(box)).reshape(-1, 3)
        color = tuple(int(c) for c in px[np.argmin(px.sum(axis=1))])  # the label's own ink
        draw.rectangle(box, fill="white")
        size = 72
        while draw.textlength(TREATMENT[code], font=font(size)) > max_width and size > 44:
            size -= 2
        draw.text(xy, TREATMENT[code], fill=color, font=font(size))
    return im


def trim(im, pad=40):
    """Crop the white margin of a supplied plot image."""
    a = np.asarray(im.convert("L"))
    ys, xs = np.nonzero(a < 245)
    return im.crop((max(xs.min() - pad, 0), max(ys.min() - pad, 0),
                    min(xs.max() + pad, im.width), min(ys.max() + pad, im.height)))


def save(im, name, dpi):
    im.save(f"{FIGDIR}/{name}.png", dpi=(dpi, dpi))
    im.save(f"{FIGDIR}/{name}.tif", compression="tiff_lzw", dpi=(dpi, dpi))
    print(f"  saved {FIGDIR}/{name}.png + .tif  {im.size}")


save(trim(relabel(SRC + "496398_0.png", XRD_LABELS, max_width=650)), "fig02_xrd", 600)
save(trim(relabel(SRC + "496397_0.png", FTIR_LABELS, max_width=900)), "fig03_ftir", 600)
SEM = SRC + "SEM TiO2/69.06.09 (IT800-0128)/SEM/TiO2-No.13-{}kx-1.bmp"
save(Image.open(SEM.format(15)).convert("RGB"), "fig04_sem_15kx", 300)
save(Image.open(SEM.format(30)).convert("RGB"), "fig05_sem_30kx", 300)
