"""Figure: an ordered (circular) genome — the ZOMBI-1-style gene-order model.

ZOMBI2's OrderedGenome represents a chromosome as a *circular* sequence of genes,
each in a gene family and carrying a strand orientation. We draw a clean starting
genome of 12 genes, each in its own family, as a ring of block-arrows (colour =
family, arrow direction = strand).

Also provides the shared ring drawer (``draw_ring``) and the muted palette used by
the rearrangement-event figures.

Run:  python figures/scripts/fig_genome_circular.py
"""

from __future__ import annotations

import math
from collections import namedtuple
from pathlib import Path

import cairosvg
import drawsvg as draw

from zombi_style import FONT, INK

OUT_STEM = Path(__file__).resolve().parent.parent / "genome_circular" / "genome_circular"

W = H = 760
CX, CY, R = W / 2, H / 2 - 6, 212
BODY_HW, HEAD_HW, HEAD_LEN = 11.0, 18.0, 14.0        # arrow radial half-widths, head length

# muted / "dead" qualitative palette (14 dusty low-saturation colours)
MUTED = ["#a3807e", "#7e8aa3", "#83a380", "#a3967e", "#7e9ca3", "#957e95",
         "#a3877e", "#8a8a66", "#8f83a3", "#679380", "#a37e8c", "#7ea390",
         "#8a7f6a", "#6a7f8a"]

Gene = namedtuple("Gene", "family orientation")
ORIENTS = [1, -1, 1, 1, -1, 1, 1, -1, 1, -1, 1, 1]   # a fixed mix of strands


def initial_genes(n=12):
    """A clean starting genome: n genes, each its own family (1..n), mixed strands."""
    return [Gene(str(i + 1), ORIENTS[i % len(ORIENTS)]) for i in range(n)]


def arrow_points(cx_g, cy_g, a, orient, half_l):
    """Block-arrow polygon (list of x,y) tangent to the ring at angle a, pointing along strand."""
    tx, ty = -math.sin(a), math.cos(a)          # clockwise tangent
    rx, ry = math.cos(a), math.sin(a)           # outward radial
    ux, uy = tx * orient, ty * orient           # arrow points with the strand
    bh = half_l - HEAD_LEN
    local = [(-half_l, -BODY_HW), (bh, -BODY_HW), (bh, -HEAD_HW), (half_l, 0.0),
             (bh, HEAD_HW), (bh, BODY_HW), (-half_l, BODY_HW)]
    pts = []
    for lu, lv in local:
        pts += [cx_g + lu * ux + lv * rx, cy_g + lu * uy + lv * ry]
    return pts


def _annular_sector(cx, cy, ri, ro, a0, a1, fill):
    p0o = (cx + ro * math.cos(a0), cy + ro * math.sin(a0))
    p1o = (cx + ro * math.cos(a1), cy + ro * math.sin(a1))
    p1i = (cx + ri * math.cos(a1), cy + ri * math.sin(a1))
    p0i = (cx + ri * math.cos(a0), cy + ri * math.sin(a0))
    large = 1 if abs(a1 - a0) > math.pi else 0
    path = draw.Path(fill=fill, stroke="none")
    path.M(*p0o).A(ro, ro, 0, large, 1, *p1o).L(*p1i).A(ri, ri, 0, large, 0, *p0i).Z()
    return path


def draw_ring(d, genes, fill_of, cx=CX, cy=CY, r=R, center_text="Ordered genome",
              label_off=34, highlight=None, hl_fill="#ededed", fs_center=21, fs_label=14):
    """Draw the gene-arrows around a ring, with family labels. ``highlight`` is a set of
    (contiguous) gene indices to shade — used to mark the segment an event acts on."""
    n = len(genes)
    if highlight:
        idx = sorted(highlight)
        a0 = math.radians(-90 + (idx[0] - 0.5) * 360.0 / n)
        a1 = math.radians(-90 + (idx[-1] + 0.5) * 360.0 / n)
        d.append(_annular_sector(cx, cy, r - 26, r + 26, a0, a1, hl_fill))
    d.append(draw.Circle(cx, cy, r, fill="none", stroke="#e2e2e2", stroke_width=1.2))
    half_l = (2 * math.pi * r / n) * 0.40
    for i, gene in enumerate(genes):
        a = math.radians(-90 + i * 360.0 / n)
        cxg, cyg = cx + r * math.cos(a), cy + r * math.sin(a)
        d.append(draw.Lines(*arrow_points(cxg, cyg, a, gene.orientation, half_l),
                            close=True, fill=fill_of(gene.family), stroke=INK, stroke_width=1.3))
        lx, ly = cx + (r + label_off) * math.cos(a), cy + (r + label_off) * math.sin(a)
        d.append(draw.Text(str(gene.family), fs_label, lx, ly, font_family=FONT, text_anchor="middle",
                           dominant_baseline="middle", fill=INK))
    if center_text:
        d.append(draw.Text(center_text, fs_center, cx, cy, font_family=FONT, text_anchor="middle",
                           dominant_baseline="middle", font_weight="bold", fill=INK))


def color_map(genes):
    """Stable family -> muted colour over a list (or lists) of genes."""
    order = list(dict.fromkeys(g.family for g in genes))
    return {f: MUTED[i % len(MUTED)] for i, f in enumerate(order)}


def main():
    genes = initial_genes(12)
    fam_color = color_map(genes)

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    draw_ring(d, genes, lambda f: fam_color[f])

    OUT_STEM.parent.mkdir(parents=True, exist_ok=True)
    OUT_STEM.with_suffix(".svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(OUT_STEM.with_suffix(".png")),
                     scale=300 / 72.0)
    print(f"wrote {OUT_STEM}.svg / .png  ({len(genes)} genes, each its own family)")


if __name__ == "__main__":
    main()
