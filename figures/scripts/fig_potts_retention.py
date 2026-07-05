"""Figure: blind gain + selective retention writes J into the profiles (Potts, idea 4).

Gain in the coupled model is field-BLIND: a gene arrives by horizontal transfer regardless of
the couplings. The couplings act afterwards, on RETENTION. When a transferred family's partners
are already present, its local field is large, its loss rate is small, and the new gene is kept;
when the partners are absent, it has no protection and is purged fast. Two scenes make the
contrast: the same HGT event into F1, once with partners present (kept) and once with partners
absent (purged). Differential retention of transferred genes is what writes the coupling J into
the observed presence/absence profiles.

House style: colour (didactic), one centered bold title, ASCII text, legend clear of the data.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_potts_retention.py
"""

from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))   # zombi_style

import cairosvg
import drawsvg as draw

from zombi_style import FONT, INK, MUTED, MODULE_COLORS, COOCCUR, AVOID, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent
NAME = "potts_retention"

W, H = 1120, 640
MODA = MODULE_COLORS[0]                           # module 1 (same palette as the panel)
POS, NEG = COOCCUR, AVOID                         # kept / purged

CELL = 62
GAP = 26


def cell(d, x, y, present, color, s=CELL):
    if present:
        d.append(draw.Rectangle(x, y, s, s, fill=color, stroke=INK, stroke_width=2.0))
    else:
        d.append(draw.Rectangle(x, y, s, s, fill="white", stroke="#bfc6cf", stroke_width=2.0))


def ghost(d, x, y, color, s=CELL):
    """A transferred-but-purged family: faint, dashed outline."""
    d.append(draw.Rectangle(x, y, s, s, fill=color, fill_opacity=0.28, stroke="#bfc6cf",
                            stroke_width=2.0, stroke_dasharray="6,4"))


def check(d, x, y, color):
    d.append(draw.Lines(x - 10, y + 2, x - 2, y + 12, x + 13, y - 10, close=False, fill="none",
                        stroke=color, stroke_width=5, stroke_linecap="round",
                        stroke_linejoin="round"))


def cross(d, x, y, color):
    for a, b in (((-9, -9), (9, 9)), ((-9, 9), (9, -9))):
        d.append(draw.Line(x + a[0], y + a[1], x + b[0], y + b[1], stroke=color, stroke_width=5,
                           stroke_linecap="round"))


def scene(d, x0, y0, partners_present, kept, label):
    """Three-cell genome F0 [F1 arrives by HGT] F2 with an outcome marker on the right."""
    # family 1
    cell(d, x0, y0, partners_present, MODA)
    d.append(draw.Text("1", FS_TICK, x0 + CELL / 2, y0 + CELL + 26, font_family=FONT,
                       text_anchor="middle", fill=INK))
    # family 2 (the transferred family)
    xm = x0 + CELL + GAP
    if kept:
        cell(d, xm, y0, True, MODA)
    else:
        ghost(d, xm, y0, MODA)
    d.append(draw.Text("2", FS_TICK, xm + CELL / 2, y0 + CELL + 26, font_family=FONT,
                       text_anchor="middle", fill=INK))
    # family 3
    x2 = x0 + 2 * (CELL + GAP)
    cell(d, x2, y0, partners_present, MODA)
    d.append(draw.Text("3", FS_TICK, x2 + CELL / 2, y0 + CELL + 26, font_family=FONT,
                       text_anchor="middle", fill=INK))

    # HGT arrow into F1 from above
    ax = xm + CELL / 2
    d.append(draw.Line(ax, y0 - 52, ax, y0 - 10, stroke=INK, stroke_width=3))
    d.append(draw.Lines(ax, y0 - 6, ax - 9, y0 - 20, ax + 9, y0 - 20, close=True, fill=INK))
    d.append(draw.Text("HGT", FS_TICK, ax + 22, y0 - 34, font_family=FONT, text_anchor="start",
                       fill=MUTED))

    # outcome marker
    ox = x2 + CELL + 74
    oy = y0 + CELL / 2
    if kept:
        check(d, ox, oy, POS)
        d.append(draw.Text("kept", FS_LABEL, ox + 24, oy, font_family=FONT, text_anchor="start",
                           dominant_baseline="central", font_weight="bold", fill=POS))
    else:
        cross(d, ox, oy, NEG)
        d.append(draw.Text("purged", FS_LABEL, ox + 24, oy, font_family=FONT, text_anchor="start",
                           dominant_baseline="central", font_weight="bold", fill=NEG))

    # per-scene caption
    d.append(draw.Text(label, FS_TICK, x0, y0 + CELL + 60, font_family=FONT, text_anchor="start",
                       fill=MUTED))


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Blind gain, selective retention", FS_TITLE, W / 2, 52,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))
    d.append(draw.Text("gain is field-blind HGT; the coupling acts on whether the new gene is kept",
                       FS_ANNOT, W / 2, 92, font_family=FONT, text_anchor="middle", fill=MUTED))

    x0 = 300
    scene(d, x0, 220, True, True,
          "module complete (partners 1, 3 present): large field, gene 2 is retained")
    scene(d, x0, 420, False, False,
          "module incomplete (partners 1, 3 absent): no field, gene 2 is lost fast")

    d.append(draw.Text("differential retention of transferred genes writes the coupling J into the profiles",
                       FS_TICK, W / 2, 588, font_family=FONT, text_anchor="middle", fill=MUTED))

    out = OUT_DIR / NAME
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{NAME}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{NAME}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{NAME}.svg / .png")


if __name__ == "__main__":
    render()
