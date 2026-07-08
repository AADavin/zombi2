"""Figure: the three archetypes of a coevolution null -- coupled vs neutral vs CID.

The didactic centrepiece for null models (docs/guide/coevolution_nulls.md). Three schematic but
*proper bifurcating* trees for the traits:species edge; a tip chip is the OBSERVED binary trait,
and a fast-diversifying clade is drawn densely branched with heavy branches.

  * COUPLED (the claim): the trait fills the fast clade, so a raw BiSSE fit reads it as causal.
  * NEUTRAL null (the strawman): a balanced tree -- no fast clade at all, nothing to explain.
  * CID null (the honest test): the SAME tree as the coupled panel, but the trait is scattered
    across fast and slow clades, so a trustworthy detector should stay quiet.

House style: one centred title, ASCII text, colour + preserved B&W. ON = trait present (filled).

Run:  python figures/scripts/fig_coevolve_nulls.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw

from zombi_style import (FONT, INK, MUTED, STATE_ON, FS_TITLE, FS_LABEL, FS_TICK)

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1200, 560
GREY = "#9a9a9a"
W_FAST, W_SLOW, W_CONN = 4.2, 2.3, 2.3
Y_TOP, Y_GAP = 104, 31                     # 8 tips: y = 104 .. 104 + 7*31 = 321
ON_COL = STATE_ON                          # chip fill (observed trait); -> INK in B&W


# --------------------------------------------------------------------------- tree model + drawer
def _assign(node, ctr):
    """Post-order: give each leaf an increasing y-index (tip order) and each internal node the
    mean of its children; internal x = its stored depth, leaves sit at the tips (depth 1)."""
    if "kids" not in node:
        node["_x"], node["_y"] = 1.0, float(ctr[0])
        ctr[0] += 1
    else:
        node["_x"] = node["x"]
        for k in node["kids"]:
            _assign(k, ctr)
        node["_y"] = sum(k["_y"] for k in node["kids"]) / len(node["kids"])


def _leaves(node, out):
    if "kids" not in node:
        out.append(node)
    else:
        for k in node["kids"]:
            _leaves(k, out)


def _line(d, x1, y1, x2, y2, w):
    d.append(draw.Line(x1, y1, x2, y2, stroke=INK, stroke_width=w, stroke_linecap="round"))


def _chip(d, cx, cy, on, s=10):
    d.append(draw.Rectangle(cx - s, cy - s, 2 * s, 2 * s,
                            fill=(ON_COL if on else "white"), stroke=INK, stroke_width=1.6))


def _draw(d, node, parent_px, X, chip_x, fast=False):
    """Rectangular cladogram: a horizontal branch into each node + a vertical connector across an
    internal node's children. ``fast`` (inherited) thickens a fast-diversifying clade."""
    fast = fast or node.get("fast", False)
    x, y = X(node["_x"]), Y_TOP + node["_y"] * Y_GAP
    _line(d, parent_px, y, x, y, W_FAST if fast else W_SLOW)          # incoming branch
    if "kids" in node:
        ys = [Y_TOP + k["_y"] * Y_GAP for k in node["kids"]]
        _line(d, x, min(ys), x, max(ys), W_CONN)                     # vertical connector
        for k in node["kids"]:
            _draw(d, k, x, X, chip_x, fast)
    else:
        _chip(d, chip_x, y, node["chip"])


def _panel(d, ox, letter, title, tree, chips, cap1, cap2):
    _assign(tree, [0])
    leaves = []
    _leaves(tree, leaves)
    for lf, c in zip(leaves, chips):
        lf["chip"] = c
    cx = ox + 150
    X = lambda depth: ox + 44 + depth * 196                          # noqa: E731  depth -> px
    chip_x = ox + 44 + 196 + 16
    # header: panel letter far left, short title centred (kept clear of the top chip)
    d.append(draw.Text(letter, FS_LABEL, ox + 6, 76, font_family=FONT, text_anchor="start",
                       fill=INK, font_weight="bold"))
    d.append(draw.Text(title, FS_LABEL, cx, 76, font_family=FONT, text_anchor="middle",
                       fill=INK, font_weight="bold"))
    _draw(d, tree, ox + 18, X, chip_x)                               # root stem starts at ox+18
    d.append(draw.Text(cap1, FS_TICK, cx, 388, font_family=FONT, text_anchor="middle",
                       fill=INK, font_style="italic"))
    d.append(draw.Text(cap2, FS_TICK, cx, 416, font_family=FONT, text_anchor="middle", fill=MUTED))


# --------------------------------------------------------------------------- the three trees
def _ref_tree():
    """Bushy 5-tip fast clade (heavy, forks packed near the tips) over a sparse 3-tip slow clade."""
    return {"x": 0.14, "kids": [
        {"x": 0.36, "fast": True, "kids": [
            {"x": 0.66, "kids": [{}, {}]},
            {"x": 0.55, "kids": [{}, {"x": 0.80, "kids": [{}, {}]}]},
        ]},
        {"x": 0.30, "kids": [
            {},
            {"x": 0.52, "kids": [{}, {}]},
        ]},
    ]}


def _balanced_tree():
    def quad():
        return {"x": 0.46, "kids": [{"x": 0.73, "kids": [{}, {}]}, {"x": 0.73, "kids": [{}, {}]}]}
    return {"x": 0.14, "kids": [quad(), quad()]}


T = True
F = False


def render(bw=False):
    global ON_COL
    ON_COL = INK if bw else STATE_ON

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("A null keeps the tree's variation -- it cuts only the trait's grip on it",
                       FS_TITLE, W / 2, 44, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))

    _panel(d, 30, "A", "COUPLED", _ref_tree(), [T, T, T, T, T, F, F, F],
           "trait fills the fast clade", "the claim: looks causal")
    _panel(d, 430, "B", "NEUTRAL null", _balanced_tree(), [T, F, F, T, F, T, F, T],
           "no fast clade at all", "the strawman: nothing to explain")
    _panel(d, 830, "C", "CID null", _ref_tree(), [T, F, T, F, T, F, T, F],
           "same fast clade, trait scattered", "the honest, worthy opponent")

    # legend
    ly = 470
    d.append(draw.Rectangle(W / 2 - 470, ly - 10, 18, 18, fill=ON_COL, stroke=INK, stroke_width=1.6))
    d.append(draw.Text("trait present", FS_TICK, W / 2 - 444, ly + 5, font_family=FONT,
                       text_anchor="start", fill=INK))
    d.append(draw.Rectangle(W / 2 - 250, ly - 10, 18, 18, fill="white", stroke=INK, stroke_width=1.6))
    d.append(draw.Text("trait absent", FS_TICK, W / 2 - 224, ly + 5, font_family=FONT,
                       text_anchor="start", fill=INK))
    _line(d, W / 2 - 40, ly - 1, W / 2 - 8, ly - 1, W_FAST)
    d.append(draw.Text("fast-diversifying (heavy, bushy) clade", FS_TICK, W / 2 + 4, ly + 5,
                       font_family=FONT, text_anchor="start", fill=INK))
    d.append(draw.Text("Panels A and C are the same tree: identical heterogeneity -- the only "
                       "difference is whether the trait tracks it.",
                       FS_TICK, W / 2, 524, font_family=FONT, text_anchor="middle", fill=MUTED))

    name = "coevolve_null_archetypes"
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
