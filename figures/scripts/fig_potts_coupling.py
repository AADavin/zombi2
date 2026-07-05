"""Figure: couplings J tie families into modules (Potts, idea 2).

The pairwise couplings J are the heart of the Potts model. A positive coupling (J > 0) makes
two families tend to co-occur; a negative coupling (J < 0) makes them avoid each other. Draw
the families as nodes and the couplings as edges and the panel organises itself into MODULES:
tight cliques of mutually positive couplings (here module A and module B), optionally linked by
a weak negative coupling that keeps the two modules apart.

House style: colour (didactic), one centered bold title, ASCII text, legend clear of the data.

NOTE ON THE FILE STEM: the stem "potts_modules" is already taken by the (distinct) emergent
co-occurrence-signal figure, so this coupling-graph figure uses the stem "potts_coupling".

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_potts_coupling.py
"""

from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))   # zombi_style

import cairosvg
import drawsvg as draw

from zombi_style import (FONT, INK, MUTED, MODULE_COLORS, COOCCUR, AVOID,
                         FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK)

OUT_DIR = Path(__file__).resolve().parent.parent
NAME = "potts_coupling"

W, H = 1080, 620
MODA, MODB = MODULE_COLORS[0], MODULE_COLORS[1]   # module 1 / module 2 (same palette as panel)
POS, NEG = COOCCUR, AVOID                          # J>0 co-occur, J<0 avoid
COL = {0: MODA, 1: MODA, 2: MODA, 3: MODB, 4: MODB, 5: MODB}
LBL = {0: "1", 1: "2", 2: "3", 3: "1", 4: "2", 5: "3"}   # within-module family number

NR = 34                                            # node radius


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Couplings J tie families into modules", FS_TITLE, W / 2, 52,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))

    # legend, centered under the title, clear of the graph
    ly = 104
    cx = W / 2
    d.append(draw.Line(cx - 300, ly, cx - 262, ly, stroke=POS, stroke_width=4.5))
    d.append(draw.Text("J > 0  (partners co-occur)", FS_TICK, cx - 252, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))
    d.append(draw.Line(cx + 60, ly, cx + 98, ly, stroke=NEG, stroke_width=4.0,
                       stroke_dasharray="9,6"))
    d.append(draw.Text("J < 0  (partners avoid)", FS_TICK, cx + 108, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))

    # node positions: two triangular cliques, well separated
    pos = {0: (300, 260), 1: (230, 380), 2: (390, 378),
           3: (850, 260), 4: (780, 380), 5: (940, 378)}

    # positive edges within each module
    for a, b in [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5)]:
        d.append(draw.Line(*pos[a], *pos[b], stroke=POS, stroke_width=4.5))
    # a negative coupling linking the two modules (keeps them apart)
    d.append(draw.Line(*pos[2], *pos[4], stroke=NEG, stroke_width=4.0, stroke_dasharray="9,6"))

    # nodes on top of the edges
    for i, (x, y) in pos.items():
        d.append(draw.Circle(x, y, NR, fill=COL[i], stroke=INK, stroke_width=2.2))
        d.append(draw.Text(LBL[i], FS_LABEL, x, y, font_family=FONT, text_anchor="middle",
                           dominant_baseline="central", fill="white", font_weight="bold"))

    # spell out the graph's semantics: a node is a family, an edge is a coupling J
    d.append(draw.Line(300, 205, 300, 224, stroke=MUTED, stroke_width=1.2))
    d.append(draw.Text("a gene family", FS_TICK, 300, 196, font_family=FONT,
                       text_anchor="middle", fill=MUTED))
    d.append(draw.Line(150, 306, 262, 320, stroke=MUTED, stroke_width=1.2))
    d.append(draw.Text("a coupling J", FS_TICK, 146, 302, font_family=FONT,
                       text_anchor="end", fill=MUTED))

    # module labels below each clique
    d.append(draw.Text("module 1", FS_LABEL, 310, 460, font_family=FONT, text_anchor="middle",
                       fill=MODA, font_weight="bold"))
    d.append(draw.Text("module 2", FS_LABEL, 860, 460, font_family=FONT, text_anchor="middle",
                       fill=MODB, font_weight="bold"))

    # caption, centered along the bottom
    d.append(draw.Text("families are nodes, couplings are edges; mutually positive J form a module",
                       FS_TICK, W / 2, 540, font_family=FONT, text_anchor="middle", fill=MUTED))
    d.append(draw.Text("a weak J < 0 between modules keeps them apart",
                       FS_TICK, W / 2, 568, font_family=FONT, text_anchor="middle", fill=MUTED))

    out = OUT_DIR / NAME
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{NAME}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{NAME}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{NAME}.svg / .png")


if __name__ == "__main__":
    render()
