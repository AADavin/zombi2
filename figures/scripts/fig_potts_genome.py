"""Figure: a genome is a presence/absence vector over a fixed family panel (Potts, idea 1).

In ZOMBI2's coupled ("Potts") gene-family model, a genome is not a bag of independent genes:
it is a fixed-length present/absent vector sigma over a panel of N gene families. Family i is
either present (sigma_i = 1) or absent (sigma_i = 0). Families are coloured by the pathway
module they belong to; the couplings that tie a module together are the subject of the
companion "coupling graph" figure.

House style: colour (didactic), one centered bold title, ASCII text, legend clear of the data.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_potts_genome.py
"""

from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))   # zombi_style

import cairosvg
import drawsvg as draw

from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent
NAME = "potts_genome"

W, H = 1080, 520
MODA, MODB = "#4477AA", "#E08A3C"                 # two pathway modules
COL = {0: MODA, 1: MODA, 2: MODA, 3: MODB, 4: MODB, 5: MODB}
SIGMA = [1, 1, 0, 1, 0, 0]

CELL = 74                                          # panel-cell size (roomy)
GAP = 22


def cell(d, x, y, present, color, s=CELL):
    if present:
        d.append(draw.Rectangle(x, y, s, s, fill=color, stroke=INK, stroke_width=2.0))
    else:
        d.append(draw.Rectangle(x, y, s, s, fill="white", stroke="#bfc6cf", stroke_width=2.0))


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("A genome is a present / absent vector", FS_TITLE, W / 2, 52,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))

    # legend: present / absent chips, centered under the title, clear of the data
    ly = 96
    cx = W / 2
    cell(d, cx - 168, ly - 13, 1, "#9aa7b4", s=26)
    d.append(draw.Text("present (1)", FS_TICK, cx - 132, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))
    cell(d, cx + 32, ly - 13, 0, "#9aa7b4", s=26)
    d.append(draw.Text("absent (0)", FS_TICK, cx + 68, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))

    n = len(SIGMA)
    row_w = n * CELL + (n - 1) * GAP
    x0 = (W - row_w) / 2
    y = 190

    for i, on in enumerate(SIGMA):
        x = x0 + i * (CELL + GAP)
        cell(d, x, y, on, COL[i])
        # family label
        d.append(draw.Text(f"F{i}", FS_LABEL, x + CELL / 2, y + CELL + 34, font_family=FONT,
                           text_anchor="middle", fill=INK))
        # the bit itself
        d.append(draw.Text(str(on), FS_ANNOT, x + CELL / 2, y + CELL + 66, font_family=FONT,
                           text_anchor="middle", font_weight="bold", fill=MUTED))

    # sigma bracket / label to the left of the row
    d.append(draw.Text("sigma", FS_LABEL, x0 - 34, y + CELL / 2, font_family=FONT,
                       text_anchor="end", dominant_baseline="central", fill=INK,
                       font_style="italic"))

    # module colour key + caption, centered below
    my = y + CELL + 120
    d.append(draw.Rectangle(cx - 250, my - 11, 22, 22, fill=MODA, stroke=INK, stroke_width=1.4))
    d.append(draw.Text("module A", FS_TICK, cx - 220, my, font_family=FONT, text_anchor="start",
                       dominant_baseline="central", fill=INK))
    d.append(draw.Rectangle(cx - 60, my - 11, 22, 22, fill=MODB, stroke=INK, stroke_width=1.4))
    d.append(draw.Text("module B", FS_TICK, cx - 30, my, font_family=FONT, text_anchor="start",
                       dominant_baseline="central", fill=INK))
    d.append(draw.Text("colour = pathway module", FS_TICK, cx + 130, my, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=MUTED))
    d.append(draw.Text("a fixed panel of N gene families; each is present (1) or absent (0)",
                       FS_TICK, W / 2, my + 40, font_family=FONT, text_anchor="middle", fill=MUTED))

    out = OUT_DIR / NAME
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{NAME}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{NAME}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{NAME}.svg / .png")


if __name__ == "__main__":
    render()
