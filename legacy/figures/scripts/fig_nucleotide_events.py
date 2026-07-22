"""Figure: nucleotide-level events on a 100-base genome (linear view of the circle).

Two events applied in sequence to the same genome, each indicated by an arrow (no
gradient ribbon): an INVERSION (a segment is reversed and its strands flip) and a
DUPLICATION (a segment is copied in tandem, so the genome grows). Bars are painted
by ancestral position, so an inversion shows as a reversed gradient and the two
duplicate copies share the same colours.

Rendered in colour and in black & white.
Run:  python figures/scripts/fig_nucleotide_events.py
"""

from __future__ import annotations

from pathlib import Path

import cairosvg
import drawsvg as draw

from fig_nucleotide import hexc, pos_color, ramp, strand_arrow  # noqa: F401
from zombi_style import FONT, INK

OUT_DIR = Path(__file__).resolve().parent.parent

CW, CH = 1040, 700                             # taller aspect (~20% larger on the page)
XL, XR = 150, 980
LMAX = 115                                     # longest bar (after the duplication)
BAR_H = 40
YA, YB, YC = 120, 360, 590                     # bar tops
INV_A, INV_B = 30, 55                          # inverted segment [a,b)
DUP_A, DUP_B = 10, 25                          # duplicated segment [c,d)


def x_of(p):
    return XL + p * (XR - XL) / LMAX


def invert(cells, a, b):
    return cells[:a] + [(sp, -st) for sp, st in reversed(cells[a:b])] + cells[b:]


def duplicate(cells, a, b):                    # tandem: copy [a,b) inserted right after b
    return cells[:b] + cells[a:b] + cells[b:]


def paint_bar(d, cells, y, mono):
    dx = (XR - XL) / LMAX
    for q, (sp, _st) in enumerate(cells):
        d.append(draw.Rectangle(x_of(q), y, dx + 0.6, BAR_H, fill=pos_color((sp + 0.5) / 100, mono)))
    d.append(draw.Rectangle(XL, y, len(cells) * dx, BAR_H, fill="none", stroke=INK, stroke_width=1.4))


def outline(d, q0, q1, y, w=2.4):
    d.append(draw.Rectangle(x_of(q0), y - 2.5, (q1 - q0) * (XR - XL) / LMAX, BAR_H + 5,
                            fill="none", stroke=INK, stroke_width=w))


def down_arrow(d, cx, y0, y1, label):
    d.append(draw.Line(cx, y0, cx, y1 - 9, stroke=INK, stroke_width=2.4))
    d.append(draw.Lines(cx, y1, cx - 6, y1 - 10, cx + 6, y1 - 10, close=True, fill=INK))
    d.append(draw.Text(label, 22, cx + 20, (y0 + y1) / 2, font_family=FONT, text_anchor="start",
                       dominant_baseline="middle", font_weight="bold", fill=INK))


def flip_arrow(d, x0, x1, y, up=16):
    """A curved 'reverse' arrow over a segment: right end wraps over to the left end."""
    p = draw.Path(fill="none", stroke=INK, stroke_width=2.2)
    p.M(x1, y).C(x1, y - up, x0, y - up, x0, y)
    d.append(p)
    d.append(draw.Lines(x0, y, x0 + 9, y - 8, x0 + 10, y + 3, close=True, fill=INK))


def render(mono):
    start = [(p, 1) for p in range(100)]
    inv = invert(start, INV_A, INV_B)
    dupd = duplicate(inv, DUP_A, DUP_B)

    d = draw.Drawing(CW, CH, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, CW, CH, fill="white"))

    # state labels on the left
    for y, txt in ((YA, "start"), (YB, "after inversion"), (YC, "after inversion\n+ duplication")):
        for k, line in enumerate(txt.split("\n")):
            d.append(draw.Text(line, 17, XL - 14, y + BAR_H / 2 - 10 + k * 20, font_family=FONT,
                               text_anchor="end", dominant_baseline="middle", fill=INK))

    paint_bar(d, start, YA, mono)
    outline(d, INV_A, INV_B, YA)
    flip_arrow(d, x_of(INV_A) + 4, x_of(INV_B) - 4, YA - 6)

    paint_bar(d, inv, YB, mono)
    outline(d, INV_A, INV_B, YB)               # the inverted segment
    outline(d, DUP_A, DUP_B, YB)               # the segment about to duplicate

    paint_bar(d, dupd, YC, mono)
    outline(d, DUP_A, DUP_B, YC)               # copy 1
    outline(d, DUP_B, DUP_B + (DUP_B - DUP_A), YC)   # copy 2 (tandem)

    down_arrow(d, x_of((INV_A + INV_B) / 2), YA + BAR_H + 6, YB - 6, "Inversion")
    down_arrow(d, x_of((DUP_A + DUP_B) / 2), YB + BAR_H + 6, YC - 6, "Duplication")

    d.append(draw.Text("Nucleotide-level events  (100-base genome, drawn linearly)", 22, XL - 40,
                       42, font_family=FONT, text_anchor="start", font_weight="bold", fill=INK))
    d.append(draw.Text("bars painted by ancestral position;  inversion = reversed gradient + "
                       "flipped strand;", 17, XL - 40, CH - 48, font_family=FONT,
                       text_anchor="start", fill="#555"))
    d.append(draw.Text("duplication = repeated colours, longer genome", 17, XL - 40, CH - 22,
                       font_family=FONT, text_anchor="start", fill="#555"))

    stem = OUT_DIR / ("nucleotide_events_bw" if mono else "nucleotide_events")
    stem.mkdir(parents=True, exist_ok=True)
    stem = stem / stem.name
    stem.with_suffix(".svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(stem.with_suffix(".png")), scale=300 / 72.0)


def main():
    for mono in (False, True):
        render(mono)
    print("wrote nucleotide_events (+_bw)")


if __name__ == "__main__":
    main()
