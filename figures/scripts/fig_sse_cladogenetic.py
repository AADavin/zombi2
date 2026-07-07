"""Figure: where trait change happens -- anagenetic (BiSSE) vs cladogenetic (ClaSSE).

The two models share the *same* Gillespie machinery and, here, the *same* tree; the only
difference is WHAT a speciation does to the trait.

  * Panel A (anagenetic, as in BiSSE): the trait changes ALONG the branches at rate q; at a
    speciation the daughters are exact copies of the parent. Open circles mark the along-branch
    changes. The amount of change on a root-to-tip path scales with TIME (branch length), so
    sister tips are usually alike.
  * Panel B (cladogenetic, the ClaSSE addition): the trait changes only AT the splits -- each
    daughter's state is drawn as part of the speciation event (filled diamond at the node);
    branches are otherwise constant. Change scales with the NUMBER of speciations, so sister
    tips can differ sharply.

Both panels are drawn on one shared tree so that "where the change lives" (mid-branch vs at the
node) is the entire message.

House style: B&W, one centered title, ASCII text.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_sse_cladogenetic.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw
import numpy as np

from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.traits import Mk, Cladogenesis, simulate_traits

from fig_trait_pagel import _layout
from model_common import zombi_to_ete3
from zombi_style import (FONT, INK, MUTED, STATE_ON, STATE_OFF,
                         FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK)

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1220, 690
GREY = "#9a9a9a"                       # state 0: light / thin branch
N_TIPS, AGE, TREE_SEED = 11, 1.0, 3

# Colour version (default -> sse_cladogenetic.svg) + preserved B&W (*_bw.svg). ON = state 1 (heavy
# branch / filled chip), OFF = state 0. Event markers (open circle = anagenetic, diamond =
# cladogenetic) stay INK — they are shape-coded per the house gene-family grammar.
ON_COL, OFF_COL = STATE_ON, STATE_OFF


def chip(d, cx, cy, on, s=13):
    """Tip chip: filled = state 1, open = state 0.  Colour-aware."""
    d.append(draw.Rectangle(cx - s, cy - s, 2 * s, 2 * s,
                            fill=ON_COL if on else "white", stroke=ON_COL, stroke_width=1.6))

Q_ANA = np.array([[-1.4, 1.4], [1.4, -1.4]])   # anagenetic rate (state flips along branches)
CLADO_SHIFT = 0.45                             # per-daughter hop probability at each split


def _nv(res):
    """{node name -> observed state index 0/1} over every node."""
    return {n.name: int(i) for n, i in res.node_values.items()}


def _sib_pairs(tree):
    """Sibling *leaf* pairs (both children are tips) -- the tests of 'do sisters agree?'."""
    ete = zombi_to_ete3(tree)
    return [(n.children[0].name, n.children[1].name) for n in ete.traverse()
            if len(n.children) == 2 and n.children[0].is_leaf() and n.children[1].is_leaf()]


def _anagenetic(tree, seed):
    """Pick a seed with a few along-branch changes, a balanced tip mix, and sisters that agree
    (the gradual-change look: recently split tips share their state)."""
    sib = _sib_pairs(tree)
    best = None
    for s in range(1, 200):
        res = simulate_traits(tree, Mk(Q_ANA), seed=s)
        nv = _nv(res)
        nc = len(list(res.changes()))
        leaves = [nv[a] for a, _ in sib] + [nv[b] for _, b in sib]
        frac1 = sum(1 for n, i in res.node_values.items() if n.is_leaf() and i == 1) / N_TIPS
        agree = sum(1 for a, b in sib if nv[a] == nv[b]) / max(1, len(sib))
        if not (3 <= nc <= 6 and 0.3 <= frac1 <= 0.7):
            continue
        score = agree - abs(frac1 - 0.5)
        if best is None or score > best[0]:
            best = (score, s, res)
    return best[2] if best else simulate_traits(tree, Mk(Q_ANA), seed=seed)


def _cladogenetic(tree, seed):
    """Pick a seed with a handful of node-jumps AND at least one sister pair that differs."""
    sib = _sib_pairs(tree)
    best = None
    for s in range(1, 120):
        res = simulate_traits(tree, Mk(np.zeros((2, 2))),
                              cladogenesis=Cladogenesis(shift=CLADO_SHIFT), seed=s)
        nv = _nv(res)
        jumps = sum(1 for n in res.node_values
                    if n.parent is not None and nv[n.name] != nv[n.parent.name])
        diff_sib = any(nv[a] != nv[b] for a, b in sib)
        if 3 <= jumps <= 7 and diff_sib and (best is None or jumps > best[0]):
            best = (jumps, s, res)
    return best[2] if best else simulate_traits(
        tree, Mk(np.zeros((2, 2))), cladogenesis=Cladogenesis(shift=CLADO_SHIFT), seed=seed)


def _panel(d, ox, oy, pw, ph, tree, res, mode, header):
    nv = _nv(res)
    changes = {}
    for node, t, frm, to in res.changes():
        changes.setdefault(node.name, []).append((t, to))

    ete = zombi_to_ete3(tree)
    tfo, present, ys, nleaf = _layout(ete)
    x_at = lambda t: ox + 40 + (t / present) * (pw - 150)      # noqa: E731
    y_at = lambda k: oy + 44 + (k / max(1, nleaf - 1)) * (ph - 104)

    d.append(draw.Text(header, FS_LABEL, ox, oy - 6, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))

    def seg(x1, x2, y, on):
        d.append(draw.Line(x1, y, x2, y, stroke=ON_COL if on else OFF_COL,
                           stroke_width=5.2 if on else 2.4, stroke_linecap="butt"))

    for n in ete.traverse():
        if n.is_root():
            continue
        y = y_at(ys[n.name])
        if mode == "ana":                                     # paint the along-branch character map
            t0, cur = tfo[n.up.name], nv[n.up.name]
            for tt, to in sorted(changes.get(n.name, [])):
                seg(x_at(t0), x_at(tt), y, cur == 1)
                d.append(draw.Circle(x_at(tt), y, 5.0, fill="white", stroke=INK, stroke_width=2.0))
                t0, cur = tt, to
            seg(x_at(t0), x_at(tfo[n.name]), y, cur == 1)
        else:                                                 # constant colour; jump lives at the node
            seg(x_at(tfo[n.up.name]), x_at(tfo[n.name]), y, nv[n.name] == 1)

    # vertical connectors, plus (cladogenetic) a jump diamond where a daughter was born different
    for n in ete.traverse("postorder"):
        if n.is_leaf():
            continue
        x = x_at(tfo[n.name])
        yy = [y_at(ys[c.name]) for c in n.children]
        d.append(draw.Line(x, min(yy), x, max(yy), stroke=INK, stroke_width=2.4))
        if mode == "clado":
            for c in n.children:
                if nv[c.name] != nv[n.name]:
                    d._draw_shape_at(x, y_at(ys[c.name]), "square",
                                     INK if nv[c.name] == 1 else "white",
                                     r=7.5, stroke=INK, stroke_width=2.0, rotation=45.0) \
                        if hasattr(d, "_draw_shape_at") else _diamond(d, x, y_at(ys[c.name]),
                                                                      nv[c.name] == 1)

    # tip chips
    colc = ox + pw - 40
    for n in ete.get_leaves():
        chip(d, colc, y_at(ys[n.name]), nv[n.name] == 1)

    # time axis
    base = oy + ph - 8
    d.append(draw.Line(x_at(0), base, x_at(present), base, stroke=INK, stroke_width=1.6))
    for k in range(5):
        tv = present * k / 4
        xx = x_at(tv)
        d.append(draw.Line(xx, base, xx, base + 6, stroke=INK, stroke_width=1.6))
        d.append(draw.Text(f"{tv:.2f}", FS_TICK, xx, base + 22, font_family=FONT,
                           text_anchor="middle", fill=INK))
    d.append(draw.Text("time (root to present)", FS_TICK, (x_at(0) + x_at(present)) / 2, base + 44,
                       font_family=FONT, text_anchor="middle", fill=MUTED))


def _diamond(d, x, y, filled):
    d.append(draw.Lines(x, y - 7.5, x + 7.5, y, x, y + 7.5, x - 7.5, y,
                        close=True, fill=INK if filled else "white",
                        stroke=INK, stroke_width=2.0))


def render(bw=False):
    global ON_COL, OFF_COL
    ON_COL, OFF_COL = (INK, GREY) if bw else (STATE_ON, STATE_OFF)

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Where the change happens: anagenetic vs cladogenetic",
                       FS_TITLE, W / 2, 44, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))

    # shared legend, centered under the title. L#1: text baselines set explicitly so each label is
    # vertically centred on its marker, and the two items are pushed well apart so they never touch.
    ly = 80
    ty = ly + 0.34 * FS_TICK
    items_x = W / 2 - 400
    d.append(draw.Circle(items_x, ly, 5.0, fill="white", stroke=INK, stroke_width=2.0))
    d.append(draw.Text("change along a branch (anagenetic)", FS_TICK, items_x + 18, ty,
                       font_family=FONT, text_anchor="start", fill=INK))
    dx = W / 2 + 110
    _diamond(d, dx, ly, True)
    d.append(draw.Text("change at a split (cladogenetic)", FS_TICK, dx + 18, ty,
                       font_family=FONT, text_anchor="start", fill=INK))

    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N_TIPS, age=AGE,
                                 direction="backward", seed=TREE_SEED)
    res_a = _anagenetic(tree, seed=4)
    res_c = _cladogenetic(tree, seed=4)

    _panel(d, 70, 158, 520, 420, tree, res_a, "ana",
           "A   anagenetic (BiSSE): change accrues with time")
    _panel(d, 660, 158, 520, 420, tree, res_c, "clado",
           "B   cladogenetic (ClaSSE): change injected at the split")

    name = "sse_cladogenetic"
    suffix = "_bw" if bw else ""
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}{suffix}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}{suffix}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{name}{suffix}.svg / .png")


if __name__ == "__main__":
    render(bw=False)   # colour -> sse_cladogenetic.svg (embedded)
    render(bw=True)    # preserved B&W -> sse_cladogenetic_bw.svg
