"""Figure: a DEC biogeographic history simulated on a tree.

Companion to the DEC schematic: here the model actually runs on a species tree. Each lineage carries
a geographic range (a set of areas); along a branch the range gains areas by dispersal and loses
them by local extinction, and at each speciation the ancestral range is passed to the daughters.
The root range and every tip range are drawn as rows of coloured area cells; dispersal (gain) and
extinction (loss) events are marked where they happen on the branches.

Areas use the same colours as the DEC schematic (categorical exception to the B&W house style).

Run:  python figures/scripts/fig_dec_tree.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw

from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.traits import DEC, simulate_biogeography

from fig_trait_pagel import _layout
from model_common import zombi_to_ete3
from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1460, 620
AREA_COLOR = {"A": "#4477AA", "B": "#EE6677", "C": "#228833"}
DISP, EXT = "#2f8f4e", "#cc4b3c"                 # dispersal (gain) green, extinction (loss) red
CELL = 24
FS_AREA = 14                                     # small area letter inside a range cell / legend swatch
N_TIPS, AGE, TREE_SEED, BIO_SEED = 9, 1.0, 3, 1


def _lum(h):
    return 0.299 * int(h[1:3], 16) + 0.587 * int(h[3:5], 16) + 0.114 * int(h[5:7], 16)


def range_cells(d, cx, cy, areas):
    """A geographic range as a centered row of coloured area cells."""
    areas = sorted(areas)
    x0 = cx - len(areas) * CELL / 2
    for i, a in enumerate(areas):
        x, fill = x0 + i * CELL, AREA_COLOR[a]
        d.append(draw.Rectangle(x, cy - CELL / 2, CELL, CELL, fill=fill, stroke=INK, stroke_width=1.2))
        # smaller area letter, centred in its cell (explicit y offset = reliable vertical centring)
        d.append(draw.Text(a, FS_AREA, x + CELL / 2, cy + 0.34 * FS_AREA, font_family=FONT,
                           text_anchor="middle", font_weight="bold",
                           fill="white" if _lum(fill) < 150 else INK))


def render():
    ztree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N_TIPS, age=AGE,
                                  direction="backward", seed=TREE_SEED)
    dec = DEC(areas=["A", "B", "C"], dispersal=0.7, extinction=0.45, max_range_size=3)
    res = simulate_biogeography(ztree, dec, root_state={"A"}, seed=BIO_SEED)
    rng_of = {n.name: res.full_label(i) for n, i in res.node_values.items()}
    tip_rng = {n.name: v for n, v in res.labeled_values().items()}
    events = {}
    for node, t, frm, to in res.changes():
        events.setdefault(node.name, []).append((t, "disp" if len(to) > len(frm) else "ext"))

    tree = zombi_to_ete3(ztree)
    tfo, present, ys, nleaf = _layout(tree)
    # wider tree: span most of the page width, leaving room for root cell (left) + tip cells (right)
    ox, oy, pw, ph = 60, 150, 1420, 400
    x_at = lambda t: ox + 90 + (t / present) * (pw - 320)      # noqa: E731
    y_at = lambda k: oy + 30 + (k / max(1, nleaf - 1)) * (ph - 60)

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("A DEC history on a tree", FS_TITLE, W / 2, 46, font_family=FONT,
                       text_anchor="middle", font_weight="bold", fill=INK))

    # tree branches (plain black) + vertical connectors
    for n in tree.traverse():
        if n.is_root():
            continue
        y = y_at(ys[n.name])
        d.append(draw.Line(x_at(tfo[n.up.name]), y, x_at(tfo[n.name]), y, stroke=INK, stroke_width=2.4))
    for n in tree.traverse("postorder"):
        if not n.is_leaf():
            x = x_at(tfo[n.name])
            yy = [y_at(ys[c.name]) for c in n.children]
            d.append(draw.Line(x, min(yy), x, max(yy), stroke=INK, stroke_width=2.4))

    # dispersal / extinction event marks on the branches
    for name, evs in events.items():
        y = y_at(ys[name])
        for t, kind in evs:
            x = x_at(t)
            if kind == "disp":                              # gain: green up-triangle
                d.append(draw.Lines(x, y - 8, x - 7, y + 5, x + 7, y + 5, close=True,
                                    fill=DISP, stroke=INK, stroke_width=0.8))
            else:                                           # loss: red cross
                for dx, dy in (((-6, -6), (6, 6)), ((-6, 6), (6, -6))):
                    d.append(draw.Line(x + dx[0], y + dx[1], x + dy[0], y + dy[1],
                                       stroke=EXT, stroke_width=2.6))

    # root range (left) and tip ranges (right)
    root = tree.get_tree_root()
    range_cells(d, x_at(0) - 46, y_at(ys[root.name]), rng_of[root.name])
    for lf in tree.get_leaves():
        y = y_at(ys[lf.name])
        range_cells(d, x_at(present) + 60, y, tip_rng[lf.name])

    # time axis
    base = oy + ph + 6
    d.append(draw.Line(x_at(0), base, x_at(present), base, stroke=INK, stroke_width=1.6))
    for k in range(5):
        tv = present * k / 4
        xx = x_at(tv)
        d.append(draw.Line(xx, base, xx, base + 6, stroke=INK, stroke_width=1.6))
        d.append(draw.Text(f"{tv:.2f}", FS_TICK, xx, base + 22, font_family=FONT,
                           text_anchor="middle", fill=INK))
    d.append(draw.Text("time (root to present)", FS_TICK, (x_at(0) + x_at(present)) / 2, base + 44,
                       font_family=FONT, text_anchor="middle", fill=MUTED))

    # legend (top-left): areas + event marks
    lx, ly = 60, 92
    for i, a in enumerate("ABC"):
        x = lx + i * 66
        d.append(draw.Rectangle(x, ly - 9, 18, 18, fill=AREA_COLOR[a], stroke=INK, stroke_width=1.1))
        # legend swatch letter: small and centred in its 18x18 square (square middle y = ly)
        d.append(draw.Text(a, FS_AREA, x + 9, ly + 0.34 * FS_AREA, font_family=FONT,
                           text_anchor="middle", font_weight="bold",
                           fill="white" if _lum(AREA_COLOR[a]) < 150 else INK))
    d.append(draw.Lines(lx + 230, ly - 8, lx + 223, ly + 5, lx + 237, ly + 5, close=True,
                        fill=DISP, stroke=INK, stroke_width=0.8))
    d.append(draw.Text("dispersal (gain)", FS_TICK, lx + 246, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))
    for dx, dy in (((-6, -6), (6, 6)), ((-6, 6), (6, -6))):
        d.append(draw.Line(lx + 430 + dx[0], ly + dx[1], lx + 430 + dy[0], ly + dy[1],
                           stroke=EXT, stroke_width=2.6))
    d.append(draw.Text("extinction (loss)", FS_TICK, lx + 446, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))

    name = "dec_tree"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png")


if __name__ == "__main__":
    render()
