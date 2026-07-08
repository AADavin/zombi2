"""Figure: the three archetypes of a coevolution null -- coupled vs neutral vs CID.

The didactic centrepiece for null models (docs/guide/coevolution_nulls.md). Three schematic trees
for the traits:species edge; a tip chip is the OBSERVED binary trait, a bushy/thick clade is
fast-diversifying.

  * COUPLED (the claim): the trait fills the fast clade, so a raw BiSSE fit reads it as causal.
  * NEUTRAL null (the strawman): a balanced tree -- no fast clade at all, nothing to explain.
  * CID null (the honest test): the SAME fast clade as the coupled panel, but the trait is spread
    across fast and slow, so a trustworthy detector should stay quiet. A and C are the same tree.

House style: one centred title, ASCII text, colour + preserved B&W. ON = trait present (filled).

Run:  python figures/scripts/fig_coevolve_nulls.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw

from zombi_style import (FONT, INK, MUTED, STATE_ON, STATE_OFF,
                         FS_TITLE, FS_LABEL, FS_TICK)

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1180, 640
GREY = "#9a9a9a"
W_FAST, W_SLOW, W_MED = 4.4, 2.2, 2.8

# Chip fill = observed trait (filled = state 1). ON swaps to INK in B&W (chips stay a separate
# channel from the branch ink either way, so the observed trait reads distinctly).
ON_COL = STATE_ON


def _line(d, x1, y1, x2, y2, w):
    d.append(draw.Line(x1, y1, x2, y2, stroke=INK, stroke_width=w, stroke_linecap="round"))


def _chip(d, cx, cy, on, s=11):
    d.append(draw.Rectangle(cx - s, cy - s, 2 * s, 2 * s,
                            fill=(ON_COL if on else "white"), stroke=INK, stroke_width=1.6))


def _bushy(d, ox, chips):
    """A fast (bushy, 6-tip, thick) clade over a slow (3-tip, thin) clade. `chips` = 9 bools."""
    _line(d, ox + 15, 250, ox + 50, 250, W_MED)              # root stem
    _line(d, ox + 50, 150, ox + 50, 315, W_MED)              # root fork
    _line(d, ox + 50, 150, ox + 118, 150, W_FAST)            # fast clade stem
    _line(d, ox + 118, 100, ox + 118, 200, W_FAST)
    fast_ys = [100, 120, 140, 160, 180, 200]
    for y in fast_ys:
        _line(d, ox + 118, y, ox + 205, y, W_FAST)
    _line(d, ox + 50, 315, ox + 105, 315, W_SLOW)            # slow clade stem
    _line(d, ox + 105, 290, ox + 105, 340, W_SLOW)
    slow_ys = [290, 315, 340]
    for y in slow_ys:
        _line(d, ox + 105, y, ox + 205, y, W_SLOW)
    for y, c in zip(fast_ys + slow_ys, chips):
        _chip(d, ox + 216, y, c)


def _balanced(d, ox, chips):
    """Two equal 4-tip clades, no bushiness -- no diversification structure. `chips` = 8 bools."""
    _line(d, ox + 15, 222, ox + 50, 222, W_MED)
    _line(d, ox + 50, 165, ox + 50, 278, W_MED)
    _line(d, ox + 50, 165, ox + 108, 165, W_MED)
    _line(d, ox + 108, 118, ox + 108, 210, W_MED)
    top = [118, 148, 178, 210]
    for y in top:
        _line(d, ox + 108, y, ox + 205, y, W_MED)
    _line(d, ox + 50, 278, ox + 108, 278, W_MED)
    _line(d, ox + 108, 235, ox + 108, 325, W_MED)
    bot = [235, 265, 295, 325]
    for y in bot:
        _line(d, ox + 108, y, ox + 205, y, W_MED)
    for y, c in zip(top + bot, chips):
        _chip(d, ox + 216, y, c)


def _caption(d, cx, key, muted):
    d.append(draw.Text(key, FS_TICK, cx, 402, font_family=FONT, text_anchor="middle",
                       fill=INK, font_style="italic"))
    d.append(draw.Text(muted, FS_TICK, cx, 430, font_family=FONT, text_anchor="middle", fill=MUTED))


def _header(d, ox, cx, letter, title):
    d.append(draw.Text(letter, FS_LABEL, ox + 8, 104, font_family=FONT, text_anchor="start",
                       fill=INK, font_weight="bold"))
    d.append(draw.Text(title, FS_LABEL, cx, 104, font_family=FONT, text_anchor="middle",
                       fill=INK, font_weight="bold"))


def render(bw=False):
    global ON_COL
    ON_COL = INK if bw else STATE_ON

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("A null keeps the tree's variation -- it cuts only the trait's grip on it",
                       FS_TITLE, W / 2, 48, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))

    # three panels: coupled / neutral / cid
    oxA, oxB, oxC = 60, 460, 860
    cxA, cxB, cxC = oxA + 120, oxB + 120, oxC + 120

    _header(d, oxA, cxA, "A", "COUPLED")
    _bushy(d, oxA, [True, True, True, True, True, True, False, False, False])   # trait fills fast clade
    _caption(d, cxA, "trait fills the fast clade", "the claim: looks causal")

    _header(d, oxB, cxB, "B", "NEUTRAL null")
    _balanced(d, oxB, [True, False, True, False, False, True, False, True])     # scattered, no fast clade
    _caption(d, cxB, "no fast clade at all", "the strawman: nothing to explain")

    _header(d, oxC, cxC, "C", "CID null")
    _bushy(d, oxC, [True, False, True, False, True, False, False, True, False])  # same tree, scattered
    _caption(d, cxC, "same fast clade, trait scattered", "the honest, worthy opponent")

    # legend
    ly = 486
    d.append(draw.Rectangle(W / 2 - 470, ly - 10, 18, 18, fill=ON_COL, stroke=INK, stroke_width=1.6))
    d.append(draw.Text("trait present", FS_TICK, W / 2 - 444, ly + 5, font_family=FONT,
                       text_anchor="start", fill=INK))
    d.append(draw.Rectangle(W / 2 - 250, ly - 10, 18, 18, fill="white", stroke=INK, stroke_width=1.6))
    d.append(draw.Text("trait absent", FS_TICK, W / 2 - 224, ly + 5, font_family=FONT,
                       text_anchor="start", fill=INK))
    _line(d, W / 2 - 40, ly - 1, W / 2 - 8, ly - 1, W_FAST)
    d.append(draw.Text("fast-diversifying (bushy) clade", FS_TICK, W / 2 + 4, ly + 5,
                       font_family=FONT, text_anchor="start", fill=INK))
    d.append(draw.Text("Panels A and C are the same tree: identical heterogeneity -- the only "
                       "difference is whether the trait tracks it.",
                       FS_TICK, W / 2, 548, font_family=FONT, text_anchor="middle", fill=MUTED))

    name = "coevolve_null_archetypes"
    suffix = "_bw" if bw else ""
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}{suffix}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}{suffix}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{name}{suffix}.svg / .png")


if __name__ == "__main__":
    render(bw=False)   # colour -> coevolve_null_archetypes.svg (embedded)
    render(bw=True)    # preserved B&W
