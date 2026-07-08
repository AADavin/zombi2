"""Figure: the timing null -- the same change, moved off the speciation events.

The companion to the null archetypes (docs/guide/coevolution_nulls.md), for the two at-speciation
edges (species:traits, species:genes). Two copies of the same tree; each teal tick is one unit of
change (a trait jump, or a gene gain/loss).

  * COUPLED: change happens AT each speciation -- the ticks cluster at the nodes, so sister tips
    differ sharply (the punctuational signature).
  * TIMING null: the SAME number of ticks, spread ALONG the branches -- sisters now differ only as
    much as their shared branch length allows. The amount of change is matched analytically.

House style: one centred title, ASCII text, colour + preserved B&W.

Run:  python figures/scripts/fig_coevolve_null_timing.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw

from zombi_style import FONT, INK, MUTED, STATE_ON, FS_TITLE, FS_LABEL, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1000, 470
GREY = "#9a9a9a"
TICK = STATE_ON       # a unit of change; -> INK in B&W


def _tree(d, ox):
    """A balanced 4-tip tree. Returns nothing; branches are horizontal, forks vertical."""
    def L(x1, y1, x2, y2):
        d.append(draw.Line(x1, y1, x2, y2, stroke=INK, stroke_width=2.8, stroke_linecap="round"))
    L(ox + 20, 220, ox + 70, 220)             # root stem
    L(ox + 70, 140, ox + 70, 300)             # root fork
    L(ox + 70, 140, ox + 140, 140)            # upper clade stem
    L(ox + 140, 100, ox + 140, 180)           # upper fork
    L(ox + 140, 100, ox + 250, 100)           # tips
    L(ox + 140, 180, ox + 250, 180)
    L(ox + 70, 300, ox + 140, 300)            # lower clade stem
    L(ox + 140, 260, ox + 140, 340)           # lower fork
    L(ox + 140, 260, ox + 250, 260)
    L(ox + 140, 340, ox + 250, 340)


def _tick(d, x, y):
    """One unit of change: a short teal bar across the (horizontal) branch."""
    d.append(draw.Line(x, y - 8, x, y + 8, stroke=TICK, stroke_width=5, stroke_linecap="round"))


def render(bw=False):
    global TICK
    TICK = INK if bw else STATE_ON

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Timing null: the same change, moved off the speciations",
                       FS_TITLE, W / 2, 44, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))

    oxA, oxB = 40, 540
    cxA, cxB = oxA + 145, oxB + 145

    # --- COUPLED: ticks cluster just after each node (change AT speciation) ---
    d.append(draw.Text("COUPLED (at speciation)", FS_LABEL, cxA, 96, font_family=FONT,
                       text_anchor="middle", fill=INK, font_weight="bold"))
    _tree(d, oxA)
    for x, y in [(oxA + 80, 140), (oxA + 80, 300),
                 (oxA + 150, 100), (oxA + 150, 180), (oxA + 150, 260), (oxA + 150, 340)]:
        _tick(d, x, y)
    d.append(draw.Text("change AT the nodes", FS_TICK, cxA, 396, font_family=FONT,
                       text_anchor="middle", fill=INK, font_style="italic"))
    d.append(draw.Text("sisters differ sharply (punctuational)", FS_TICK, cxA, 424,
                       font_family=FONT, text_anchor="middle", fill=MUTED))

    # --- TIMING null: same six ticks, spread along the branches ---
    d.append(draw.Text("TIMING null (along branches)", FS_LABEL, cxB, 96, font_family=FONT,
                       text_anchor="middle", fill=INK, font_weight="bold"))
    _tree(d, oxB)
    for x, y in [(oxB + 110, 140), (oxB + 110, 300),
                 (oxB + 205, 100), (oxB + 230, 180), (oxB + 195, 260), (oxB + 215, 340)]:
        _tick(d, x, y)
    d.append(draw.Text("same amount, spread ALONG branches", FS_TICK, cxB, 396, font_family=FONT,
                       text_anchor="middle", fill=INK, font_style="italic"))
    d.append(draw.Text("sisters differ only by shared branch length", FS_TICK, cxB, 424,
                       font_family=FONT, text_anchor="middle", fill=MUTED))

    name = "coevolve_null_timing"
    suffix = "_bw" if bw else ""
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}{suffix}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}{suffix}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{name}{suffix}.svg / .png")


if __name__ == "__main__":
    render(bw=False)
    render(bw=True)
