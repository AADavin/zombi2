"""Figure: the timing null -- the same change, moved off the speciation events.

The companion to the null archetypes (docs/guide/coevolution_nulls.md), for the two at-speciation
edges (species:traits, species:genes). Two copies of the same 6-tip tree; each teal tick is one
unit of change (a trait jump, or a gene gain/loss).

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

W, H = 1000, 480
GREY = "#9a9a9a"
BR_W = 2.6
Y_OFF = 40            # push the tree down, clear of the panel titles
TICK = STATE_ON       # a unit of change; -> INK in B&W


def _br(d, ox, x1, y1, x2, y2):
    d.append(draw.Line(ox + x1, y1 + Y_OFF, ox + x2, y2 + Y_OFF, stroke=INK, stroke_width=BR_W,
                       stroke_linecap="round"))


def _tick(d, ox, x, y):
    """One unit of change: a short heavy teal bar across the (horizontal) branch."""
    d.append(draw.Line(ox + x, y + Y_OFF - 9, ox + x, y + Y_OFF + 9, stroke=TICK, stroke_width=5.4,
                       stroke_linecap="round"))


# A balanced 6-tip tree: ((t1,t2),t3) over (t4,(t5,t6)). Coordinates are local to a panel origin.
def _tree(d, ox):
    _br(d, ox, 28, 215, 75, 215)          # root stem
    _br(d, ox, 75, 150, 75, 280)          # root fork
    _br(d, ox, 75, 150, 140, 150)         # -> clade 1
    _br(d, ox, 75, 280, 140, 280)         # -> clade 2
    _br(d, ox, 140, 120, 140, 180)        # clade-1 fork
    _br(d, ox, 140, 120, 210, 120)        # -> (t1,t2)
    _br(d, ox, 140, 180, 300, 180)        # -> t3
    _br(d, ox, 210, 100, 210, 140)        # (t1,t2) fork
    _br(d, ox, 210, 100, 300, 100)        # t1
    _br(d, ox, 210, 140, 300, 140)        # t2
    _br(d, ox, 140, 250, 140, 310)        # clade-2 fork
    _br(d, ox, 140, 250, 300, 250)        # -> t4
    _br(d, ox, 140, 310, 210, 310)        # -> (t5,t6)
    _br(d, ox, 210, 290, 210, 330)        # (t5,t6) fork
    _br(d, ox, 210, 290, 300, 290)        # t5
    _br(d, ox, 210, 330, 300, 330)        # t6


# Six ticks per panel (same count = same amount of change). AT-node positions sit right at a fork
# (a daughter branch's origin, +3 px past the node x, so the tick touches the split); ALONG-branch
# positions sit at branch midpoints, spread out. The internal node x-coords are 75 (root), 140
# (both deep clade forks) and 210 (both shallow forks).
_AT_NODES = [(78, 150), (78, 280), (143, 180), (213, 100), (143, 250), (213, 290)]
_ALONG = [(108, 150), (240, 180), (258, 100), (108, 280), (238, 250), (262, 330)]


def render(bw=False):
    global TICK
    TICK = INK if bw else STATE_ON

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Timing null: the same change, moved off the speciations",
                       FS_TITLE, W / 2, 46, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))

    oxA, oxB = 40, 540
    cxA, cxB = oxA + 165, oxB + 165

    d.append(draw.Text("COUPLED (at speciation)", FS_LABEL, cxA, 100, font_family=FONT,
                       text_anchor="middle", fill=INK, font_weight="bold"))
    _tree(d, oxA)
    for x, y in _AT_NODES:
        _tick(d, oxA, x, y)
    d.append(draw.Text("change AT the nodes", FS_TICK, cxA, 400, font_family=FONT,
                       text_anchor="middle", fill=INK, font_style="italic"))
    d.append(draw.Text("sisters differ sharply (punctuational)", FS_TICK, cxA, 428,
                       font_family=FONT, text_anchor="middle", fill=MUTED))

    d.append(draw.Text("TIMING null (along branches)", FS_LABEL, cxB, 100, font_family=FONT,
                       text_anchor="middle", fill=INK, font_weight="bold"))
    _tree(d, oxB)
    for x, y in _ALONG:
        _tick(d, oxB, x, y)
    d.append(draw.Text("same amount, spread ALONG branches", FS_TICK, cxB, 400, font_family=FONT,
                       text_anchor="middle", fill=INK, font_style="italic"))
    d.append(draw.Text("sisters differ only by shared branch length", FS_TICK, cxB, 428,
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
