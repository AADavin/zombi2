"""Figure (experimental / language-model selection): the block-based genome pipeline (P2.6).

A real annotated genome is evolved down a species tree with selection on its genes:

  1. the root genome is partitioned by its GFF into GENE intervals and the INTERGENE gaps between them;
  2. the nucleotide model runs the structural simulation (inversion / duplication / loss / transfer),
     and the result is traced back into BLOCKS -- maximal never-cut intervals, each carrying its own
     gene tree. Design S keeps a gene to exactly one block, so a whole coding sequence evolves as one
     unit down one tree: a GENE block evolves under ESM2 codon selection, an INTERGENE block drifts
     neutrally;
  3. the evolved blocks are reassembled, in genome order, into the DNA at every node.

House style: one centered bold title, ASCII text, categorical colours for the two genes.
Run:  python figures/scripts/fig_selection_genome_blocks.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw

from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_TICK

OUT = Path(__file__).resolve().parent.parent.parent / "docs" / "img"
W, H = 1200, 640
GENE1, GENE2 = "#4477AA", "#7B5EA7"     # two genes (Paul Tol categorical)
INTER = "#c9c9c9"                        # intergene grey
BARH = 30


def _stage_title(d, x, y, text):
    d.append(draw.Text(text, FS_LABEL, x, y, font_family=FONT, text_anchor="start",
                       font_weight="bold", fill=INK))


def _genome_bar(d, x, y, w, segments):
    """segments: list of (fraction, color). Draws a segmented horizontal genome bar of width w."""
    cx = x
    for frac, color in segments:
        sw = w * frac
        d.append(draw.Rectangle(cx, y, sw, BARH, fill=color, stroke="white", stroke_width=1.5))
        cx += sw
    d.append(draw.Rectangle(x, y, w, BARH, fill="none", stroke=INK, stroke_width=1.2))


def _arrow(d, x0, y0, x1, y1):
    d.append(draw.Line(x0, y0, x1, y1, stroke=INK, stroke_width=2.6))
    import numpy as np
    ang = np.arctan2(y1 - y0, x1 - x0)
    for da in (0.5, -0.5):
        d.append(draw.Line(x1, y1, x1 - 14 * np.cos(ang + da), y1 - 14 * np.sin(ang + da),
                           stroke=INK, stroke_width=2.6))


def _minitree(d, x, y, w, h, color, tips=3):
    """A tiny 3-tip cladogram in [x,x+w] x [y,y+h], tips coloured squares on the right."""
    rootx = x
    midx = x + 0.42 * w
    tipx = x + w - 12
    ys = [y + h * k / (tips - 1) for k in range(tips)]
    # root -> split
    d.append(draw.Line(rootx, (ys[0] + ys[-1]) / 2, midx, (ys[0] + ys[-1]) / 2, stroke=INK, stroke_width=2.2))
    d.append(draw.Line(midx, ys[0], midx, ys[-1], stroke=INK, stroke_width=2.2))
    # a second internal split for the bottom two tips
    d.append(draw.Line(midx, ys[0], tipx, ys[0], stroke=INK, stroke_width=2.2))
    mid2 = x + 0.7 * w
    d.append(draw.Line(midx, (ys[1] + ys[2]) / 2, mid2, (ys[1] + ys[2]) / 2, stroke=INK, stroke_width=2.2))
    d.append(draw.Line(mid2, ys[1], mid2, ys[2], stroke=INK, stroke_width=2.2))
    d.append(draw.Line(mid2, ys[1], tipx, ys[1], stroke=INK, stroke_width=2.2))
    d.append(draw.Line(mid2, ys[2], tipx, ys[2], stroke=INK, stroke_width=2.2))
    for yy in ys:
        d.append(draw.Rectangle(tipx, yy - 6, 12, 12, fill=color, stroke="white", stroke_width=1))


def main():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Language-model selection on a real genome, block by block", FS_TITLE,
                       W / 2, 44, font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))

    # ---- Stage 1: root genome + GFF ----
    x1 = 50
    _stage_title(d, x1, 110, "1. Root genome + GFF")
    segs = [(0.12, INTER), (0.24, GENE1), (0.14, INTER), (0.26, GENE2), (0.24, INTER)]
    _genome_bar(d, x1, 135, 320, segs)
    d.append(draw.Text("gene A", FS_TICK, x1 + 320 * 0.24, 135 + BARH + 20, font_family=FONT,
                       text_anchor="middle", fill=GENE1, font_weight="bold"))
    d.append(draw.Text("gene B", FS_TICK, x1 + 320 * 0.63, 135 + BARH + 20, font_family=FONT,
                       text_anchor="middle", fill=GENE2, font_weight="bold"))
    d.append(draw.Text("grey = intergene", FS_TICK, x1 + 320 + 16, 135 + BARH / 2 + 6,
                       font_family=FONT, text_anchor="start", fill=MUTED))

    # ---- arrow down to Stage 2 ----
    _arrow(d, x1 + 160, 210, x1 + 160, 252)
    d.append(draw.Text("structural evolution (inversion, duplication, loss, transfer),", FS_TICK,
                       x1 + 182, 230, font_family=FONT, text_anchor="start", fill=MUTED))
    d.append(draw.Text("then trace every extant position back to its ancestral block", FS_TICK,
                       x1 + 182, 252, font_family=FONT, text_anchor="start", fill=MUTED))

    # ---- Stage 2: each block on its own gene tree (labels INSIDE each box) ----
    _stage_title(d, x1, 300, "2. Each block on its own gene tree")
    BOXW = 250
    # gene block (selection)
    bx, by = x1, 320
    d.append(draw.Rectangle(bx, by, BOXW, 150, rx=10, fill="#f4f7fb", stroke=GENE1, stroke_width=1.8))
    d.append(draw.Text("GENE block", FS_TICK, bx + 14, by + 26, font_family=FONT, text_anchor="start",
                       fill=GENE1, font_weight="bold"))
    _minitree(d, bx + 20, by + 40, 150, 78, GENE1)
    d.append(draw.Text("ESM2 codon selection", FS_TICK, bx + BOXW / 2 + 20, by + 140, font_family=FONT,
                       text_anchor="middle", fill=INK, font_weight="bold"))
    # intergene block (neutral)
    ix, iy = x1, 500
    d.append(draw.Rectangle(ix, iy, BOXW, 120, rx=10, fill="#f6f6f6", stroke=MUTED, stroke_width=1.8))
    d.append(draw.Text("INTERGENE block", FS_TICK, ix + 14, iy + 24, font_family=FONT,
                       text_anchor="start", fill=MUTED, font_weight="bold"))
    _minitree(d, ix + 20, iy + 36, 150, 56, INTER)
    d.append(draw.Text("neutral drift", FS_TICK, ix + BOXW / 2 + 20, iy + 110, font_family=FONT,
                       text_anchor="middle", fill=INK, font_weight="bold"))

    # ---- arrow across to Stage 3 ----
    ax0, ax1 = x1 + BOXW + 20, x1 + BOXW + 120
    _arrow(d, ax0, 452, ax1, 452)
    for k, t in enumerate(("evolve each", "block, then", "reassemble")):
        d.append(draw.Text(t, FS_TICK, (ax0 + ax1) / 2, 388 + k * 20, font_family=FONT,
                           text_anchor="middle", fill=MUTED))

    # ---- Stage 3: genomes at every node ----
    x3 = x1 + BOXW + 150
    _stage_title(d, x3, 300, "3. Genomes at every node")
    tips = [
        [(0.10, INTER), (0.24, GENE1), (0.16, INTER), (0.26, GENE2), (0.24, INTER)],
        [(0.14, INTER), (0.26, GENE2), (0.12, INTER), (0.22, GENE1), (0.26, INTER)],   # rearranged
        [(0.12, INTER), (0.24, GENE1), (0.40, INTER), (0.24, GENE2)],                  # gene B moved
    ]
    labels = ["node (ancestor)", "tip 1", "tip 2"]
    for k, (segs3, lab) in enumerate(zip(tips, labels)):
        yy = 345 + k * 92
        _genome_bar(d, x3, yy, 360, segs3)
        d.append(draw.Text(lab, FS_TICK, x3, yy - 8, font_family=FONT, text_anchor="start", fill=MUTED))
    d.append(draw.Text("the root reproduces the input genome exactly", FS_TICK, x3, 345 + 3 * 92 + 4,
                       font_family=FONT, text_anchor="start", fill=MUTED, font_style="italic"))

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "selection_genome_blocks.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(),
                     write_to=str(OUT / "selection_genome_blocks.png"), scale=200 / 72.0)
    print(f"wrote {OUT}/selection_genome_blocks.svg / .png")


if __name__ == "__main__":
    main()
