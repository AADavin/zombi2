"""Figure: state-dependent diversification (BiSSE) -- two panels.

A binary trait sets each lineage's speciation and extinction rate, and the tree is grown jointly
with the trait. Here state 1 speciates three times faster than state 0 (equal extinction, equal
transitions), so state-1 lineages leave more descendants and come to dominate the standing tips --
the diversification signal is written into the shape of the tree itself.

  * Panel A (the model): the two states as nodes; curved arrows are the anagenetic transitions
    (width = rate); the fork under each state is its speciation rate (width = lambda). State 1's
    fork is visibly heavier -- it branches faster.
  * Panel B (a realization): the complete simulated tree. Branches are heavy where the lineage is
    in state 1 and light in state 0; lineages that go extinct end in a small cross; the extant
    tips carry chips. The heavy (state-1) lineages proliferate and fill the present-day tips.

House style: B&W, one centered title, ASCII text.

Run:  python figures/scripts/fig_sse.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw

from zombi2.coevolve import BiSSE, simulate_sse

from fig_trait_pagel import curved_arrow, rate_width, _layout, NR
from model_common import zombi_to_ete3
from zombi_style import (FONT, INK, MUTED, STATE_ON, STATE_OFF,
                         FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK)

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1200, 700
GREY = "#9a9a9a"


def state_node(d, cx, cy, label):
    """A Markov-state node whose label is geometrically centred in the disc (M#1): the text
    baseline is set explicitly rather than relying on dominant_baseline, which renders slightly
    high in rsvg/cairosvg."""
    d.append(draw.Circle(cx, cy, NR, fill="white", stroke=INK, stroke_width=2.2))
    d.append(draw.Text(label, FS_LABEL, cx, cy + 0.34 * FS_LABEL, font_family=FONT,
                       text_anchor="middle", fill=INK, font_weight="bold"))

# The figure ships in two flavours: a colour version (default; the tasteful pale STATE_ON/OFF
# two-tone) written to <name>.svg — the copy the manual embeds — and the original B&W kept as
# <name>_bw.svg.  ON_COL/OFF_COL are swapped between the two by render(bw=...).
ON_COL, OFF_COL = STATE_ON, STATE_OFF


def chip(d, cx, cy, on, s=13):
    """Tip chip: filled = state 1 (present), open = state 0.  Colour-aware."""
    d.append(draw.Rectangle(cx - s, cy - s, 2 * s, 2 * s,
                            fill=ON_COL if on else "white", stroke=ON_COL, stroke_width=1.6))

# the model: state 1 speciates 3x faster; equal extinction; symmetric slow transitions
L0, L1, MU, Q = 1.0, 3.0, 0.3, 0.12
N_TIPS = 12


def spec_fork(d, x, y, lam):
    """A little downward fork whose stroke width encodes the speciation rate lambda."""
    w = 2.0 + 2.1 * lam
    stem, arm = 16, 22
    d.append(draw.Line(x, y, x, y + stem, stroke=INK, stroke_width=w, stroke_linecap="round"))
    d.append(draw.Line(x, y + stem, x - arm, y + stem + arm, stroke=INK, stroke_width=w,
                       stroke_linecap="round"))
    d.append(draw.Line(x, y + stem, x + arm, y + stem + arm, stroke=INK, stroke_width=w,
                       stroke_linecap="round"))


# --------------------------------------------------------------------------- panel A: the model
def panel_model(d, cx0, cy0):
    gx = 250
    # P#1: panel letter at the top-left corner; panel title centred over the panel.
    d.append(draw.Text("A", FS_LABEL, cx0 - 60, cy0 - 150, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))
    d.append(draw.Text("the model", FS_LABEL, cx0 + gx / 2, cy0 - 150, font_family=FONT,
                       text_anchor="middle", fill=INK, font_weight="bold"))
    P = {"0": (cx0, cy0), "1": (cx0 + gx, cy0)}
    curved_arrow(d, P["0"], P["1"], +1, 20, rate_width(Q), f"{Q:g}")   # q01 bows up
    curved_arrow(d, P["1"], P["0"], -1, 20, rate_width(Q), f"{Q:g}")   # q10 bows down
    for lab, (x, y) in P.items():
        state_node(d, x, y, lab)
    spec_fork(d, cx0, cy0 + 46, L0)
    spec_fork(d, cx0 + gx, cy0 + 46, L1)
    d.append(draw.Text(f"lambda = {L0:g},  mu = {MU:g}", FS_TICK, cx0, cy0 + 118,
                       font_family=FONT, text_anchor="middle", fill=INK))
    d.append(draw.Text(f"lambda = {L1:g},  mu = {MU:g}", FS_TICK, cx0 + gx, cy0 + 118,
                       font_family=FONT, text_anchor="middle", fill=INK))
    d.append(draw.Text("nodes = state;  arrows = transitions;",
                       FS_TICK, cx0 + gx / 2, cy0 + 158, font_family=FONT, text_anchor="middle",
                       fill=MUTED))
    d.append(draw.Text("fork width = speciation rate",
                       FS_TICK, cx0 + gx / 2, cy0 + 186, font_family=FONT, text_anchor="middle",
                       fill=MUTED))
    d.append(draw.Text("state 1 branches 3x faster", FS_ANNOT, cx0 + gx / 2, cy0 - 96,
                       font_family=FONT, text_anchor="middle", fill=INK, font_style="italic"))


# --------------------------------------------------------------------------- panel B: a realization
def _pick(seed_range):
    """Grow a BiSSE tree that *tells the story*: state 1 dominates the extant tips, but state 0
    is still visibly present and a few lineages go extinct -- and the complete tree stays legible."""
    model = BiSSE(lambda0=L0, lambda1=L1, mu0=MU, mu1=MU, q01=Q, q10=Q)
    best = None
    for s in seed_range:
        res = simulate_sse(model, n_tips=N_TIPS, seed=s)
        nv = {n.name: int(i) for n, i in res.node_values.items()}
        ete = zombi_to_ete3(res.tree)
        leaves = ete.get_leaves()
        extant = [n for n in leaves if n.is_extant]
        n_ext = len(leaves) - len(extant)
        n0 = sum(1 for n in extant if nv[n.name] == 0)
        f1 = 1 - n0 / max(1, len(extant))
        if len(leaves) > 22 or not (2 <= n0 <= 4) or n_ext < 1:   # both states present; some extinction
            continue
        score = -abs(f1 - 0.72) + 0.04 * n_ext
        if best is None or score > best[0]:
            best = (score, s, res, ete, nv)
    if best is None:                                             # fall back: best f1 in range
        for s in seed_range:
            res = simulate_sse(model, n_tips=N_TIPS, seed=s)
            nv = {n.name: int(i) for n, i in res.node_values.items()}
            ete = zombi_to_ete3(res.tree)
            if len(ete.get_leaves()) <= 24:
                return (0.0, s, res, ete, nv)
    return best


def panel_realization(d, ox, oy, pw, ph):
    f1, seed, res, ete, nv = _pick(range(1, 200))
    changes = {}
    for node, t, frm, to in res.changes():
        changes.setdefault(node.name, []).append((t, to))

    tfo, present, ys, nleaf = _layout(ete)
    x_at = lambda t: ox + 44 + (t / present) * (pw - 150)      # noqa: E731
    y_at = lambda k: oy + 44 + (k / max(1, nleaf - 1)) * (ph - 104)

    # P#1: panel letter at the top-left corner; panel title centred over the panel. The italic
    # "state 1 fills ..." note drops to its own row below the title so it never collides with it.
    d.append(draw.Text("B", FS_LABEL, ox, oy - 6, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))
    d.append(draw.Text("a simulated realization", FS_LABEL, ox + pw / 2, oy - 6, font_family=FONT,
                       text_anchor="middle", fill=INK, font_weight="bold"))
    n_ext = sum(1 for n in ete.get_leaves() if n.is_extant)
    n1 = sum(1 for n in ete.get_leaves() if n.is_extant and nv[n.name] == 1)
    d.append(draw.Text(f"state 1 fills {n1} of {n_ext} extant tips", FS_ANNOT, ox + pw - 34,
                       oy + 24, font_family=FONT, text_anchor="end", fill=INK, font_style="italic"))

    def seg(x1, x2, y, on):
        d.append(draw.Line(x1, y, x2, y, stroke=ON_COL if on else OFF_COL,
                           stroke_width=5.2 if on else 2.4, stroke_linecap="butt"))

    for n in ete.traverse():
        if n.is_root():
            continue
        y = y_at(ys[n.name])
        t0, cur = tfo[n.up.name], nv[n.up.name]
        for tt, to in sorted(changes.get(n.name, [])):
            seg(x_at(t0), x_at(tt), y, cur == 1)
            t0, cur = tt, to
        seg(x_at(t0), x_at(tfo[n.name]), y, cur == 1)
        if n.is_leaf() and not n.is_extant:                    # extinct tip: cross cap, no chip
            xx, ah = x_at(tfo[n.name]), 6.0
            d.append(draw.Line(xx - ah, y - ah, xx + ah, y + ah, stroke=INK, stroke_width=2.2))
            d.append(draw.Line(xx - ah, y + ah, xx + ah, y - ah, stroke=INK, stroke_width=2.2))

    for n in ete.traverse("postorder"):
        if not n.is_leaf():
            x = x_at(tfo[n.name])
            yy = [y_at(ys[c.name]) for c in n.children]
            d.append(draw.Line(x, min(yy), x, max(yy), stroke=INK, stroke_width=2.4))

    colc = ox + pw - 34
    for n in ete.get_leaves():
        if n.is_extant:
            chip(d, colc, y_at(ys[n.name]), nv[n.name] == 1, s=9)

    base = oy + ph - 8
    d.append(draw.Line(x_at(0), base, x_at(present), base, stroke=INK, stroke_width=1.6))
    for k in range(5):
        tv = present * k / 4
        xx = x_at(tv)
        d.append(draw.Line(xx, base, xx, base + 6, stroke=INK, stroke_width=1.6))
        d.append(draw.Text(f"{tv:.1f}", FS_TICK, xx, base + 22, font_family=FONT,
                           text_anchor="middle", fill=INK))
    d.append(draw.Text("time (root to present)", FS_TICK, (x_at(0) + x_at(present)) / 2, base + 44,
                       font_family=FONT, text_anchor="middle", fill=MUTED))
    return seed


def render(bw=False):
    global ON_COL, OFF_COL
    ON_COL, OFF_COL = (INK, GREY) if bw else (STATE_ON, STATE_OFF)

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("State-dependent diversification (BiSSE)", FS_TITLE, W / 2, 46,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))

    # legend -- text baselines set explicitly (L#1) so each label sits vertically centred on its marker
    ly = 82
    ty = ly + 0.34 * FS_TICK
    d.append(draw.Line(W / 2 - 300, ly, W / 2 - 268, ly, stroke=ON_COL, stroke_width=5.2))
    d.append(draw.Text("in state 1", FS_TICK, W / 2 - 256, ty, font_family=FONT,
                       text_anchor="start", fill=INK))
    d.append(draw.Line(W / 2 - 150, ly, W / 2 - 118, ly, stroke=OFF_COL, stroke_width=2.4))
    d.append(draw.Text("in state 0", FS_TICK, W / 2 - 106, ty, font_family=FONT,
                       text_anchor="start", fill=INK))
    ah = 5.0
    d.append(draw.Line(W / 2 + 10 - ah, ly - ah, W / 2 + 10 + ah, ly + ah, stroke=INK, stroke_width=2.2))
    d.append(draw.Line(W / 2 + 10 - ah, ly + ah, W / 2 + 10 + ah, ly - ah, stroke=INK, stroke_width=2.2))
    d.append(draw.Text("extinct", FS_TICK, W / 2 + 24, ty, font_family=FONT,
                       text_anchor="start", fill=INK))

    panel_model(d, 210, 320)
    seed = panel_realization(d, 560, 150, 600, 460)

    name = "sse"
    suffix = "_bw" if bw else ""
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}{suffix}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}{suffix}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{name}{suffix}.svg / .png  (tree seed {seed})")


if __name__ == "__main__":
    render(bw=False)   # colour version -> sse.svg (embedded in the manual)
    render(bw=True)    # preserved B&W  -> sse_bw.svg
