"""Figure: correlated binary characters (Pagel 1994) — two panels.

Two binary traits X and Y evolve jointly over the 4-state space {00, 01, 10, 11}, with
only one trait flipping at a time. Each trait's gain/loss rate may depend on the *other*
trait's current state, and that dependence is exactly what "correlated evolution" means.

  * Panel A (the model): the four joint states as a square. Horizontal edges flip Y, vertical
    edges flip X. Arrow WIDTH is proportional to the rate, so the coupling is visible at a
    glance: with the example rates, Y is gained fast only when X = 1 and lost fast when X = 0
    (thick arrows on the bottom / top rows) — Y is driven toward X.
  * Panel B (a realization): the same model simulated on a tree. Branches are drawn heavy
    where X = 1 (present) and light where X = 0; the two tip columns give each leaf's (X, Y).
    Y tracks X: the clades where X switches on are the clades where Y is present.

This figure is produced in two palettes from one run: a soft, low-saturation COLOUR version
(the default, trait_pagel.svg/.png) and the original B&W version (trait_pagel_bw.svg/.png). In
colour, "present" (state 1) is drawn in a muted accent hue; in B&W it is near-black. The BW look
is preserved so the monochrome figure can always be regenerated.

Run:  python figures/scripts/fig_trait_pagel.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))   # zombi_style / model_common

import cairosvg
import drawsvg as draw

from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.traits import CorrelatedBinary, simulate_traits

from model_common import zombi_to_ete3
from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1200, 660
GREY = "#9a9a9a"                      # X = 0 (absent): light branch / open chip stroke

# "present" (state 1 / X = 1) is the colour-carrying element. Two palettes; BW preserves the
# original near-black look, COLOUR uses a soft dusty-blue accent. PRESENT is set per render mode.
PRESENT_BW    = INK
PRESENT_COLOR = "#5a86a0"             # dusty blue, matches the Mk / hidden-Mk colour set
PRESENT = PRESENT_COLOR              # active "present" colour (set per mode)

# --- the example model: X neutral & symmetric; Y strongly tracks X --------------------
MODEL = dict(x_gain_y0=0.5, x_gain_y1=0.5, x_loss_y0=0.5, x_loss_y1=0.5,
             y_gain_x0=0.05, y_gain_x1=2.0, y_loss_x0=2.0, y_loss_x1=0.05)
N_TIPS, AGE, TREE_SEED = 11, 1.0, 3


# --------------------------------------------------------------------------- helpers
def rate_width(r: float) -> float:
    """Map a rate to a stroke width so fast/slow transitions read at a glance."""
    return 1.6 + 3.6 * math.sqrt(r / 2.0)


def edge_point(cx, cy, tx, ty, r):
    a = math.atan2(ty - cy, tx - cx)
    return cx + r * math.cos(a), cy + r * math.sin(a)


def curved_arrow(d, c1, c2, side, bend, width, label):
    """Curved arrow c1 -> c2 (centres), width-coded by rate; bows to `side` (+1 / -1).

    The perpendicular is canonicalized to a fixed orientation (independent of which way the
    arrow points), so `side` reliably selects the same visual side of the edge for both
    directions of a pair. Each arc's rate label sits just beyond its own bulge, so the two
    labels of a bidirectional edge never collide.
    """
    x1, y1 = c1
    x2, y2 = c2
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    dx, dy = x2 - x1, y2 - y1
    L = math.hypot(dx, dy)
    px, py = -dy / L, dx / L                        # unit normal
    if py > 0 or (abs(py) < 1e-9 and px > 0):       # canonicalize: normal points up / left
        px, py = -px, -py
    ox, oy = px * side, py * side                    # chosen side, direction-independent
    cx, cy = mx + ox * bend, my + oy * bend          # quadratic control point
    sx, sy = edge_point(x1, y1, cx, cy, NR)
    ex, ey = edge_point(x2, y2, cx, cy, NR)
    p = draw.Path(fill="none", stroke=INK, stroke_width=width, stroke_linecap="round")
    p.M(sx, sy).Q(cx, cy, ex, ey)
    d.append(p)
    ang, ah = math.atan2(ey - cy, ex - cx), 11.0
    d.append(draw.Lines(ex, ey,
                        ex - ah * math.cos(ang - 0.42), ey - ah * math.sin(ang - 0.42),
                        ex - ah * math.cos(ang + 0.42), ey - ah * math.sin(ang + 0.42),
                        close=True, fill=INK))
    if label:
        d.append(draw.Text(label, FS_TICK, cx + ox * 16, cy + oy * 16, font_family=FONT,
                           text_anchor="middle", dominant_baseline="central", fill=MUTED))


NR = 34                                            # state-node radius


def state_node(d, cx, cy, label):
    d.append(draw.Circle(cx, cy, NR, fill="white", stroke=INK, stroke_width=2.2))
    d.append(draw.Text(label, FS_LABEL, cx, cy, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill=INK, font_weight="bold"))


def chip(d, cx, cy, on, s=13, fill_on=None):
    """A tip chip: filled = state 1 (present), open = state 0 (absent).

    `fill_on` sets the fill for the "on" state; it defaults to the active PRESENT colour so
    the colour / B&W palette flows through. Passing it explicitly lets other figures pick their
    own accent without touching this module's global.
    """
    fill = (fill_on if fill_on is not None else PRESENT) if on else "white"
    d.append(draw.Rectangle(cx - s, cy - s, 2 * s, 2 * s,
                            fill=fill, stroke=INK, stroke_width=1.6))


# --------------------------------------------------------------------------- panel A: the model
def panel_model(d, cx0, cy0):
    """4-state square. (cx0, cy0) is the top-left state's centre. Paired arrows bow apart."""
    gx = gy = 176
    # P#1: panel letter top-left of the panel; title horizontally centred over the panel
    d.append(draw.Text("A", FS_LABEL, cx0 - 40, cy0 - 146, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))
    d.append(draw.Text("the model", FS_LABEL, cx0 + gx / 2, cy0 - 146, font_family=FONT,
                       text_anchor="middle", fill=INK, font_weight="bold"))
    P = {"00": (cx0, cy0), "01": (cx0 + gx, cy0),
         "10": (cx0, cy0 + gy), "11": (cx0 + gx, cy0 + gy)}
    m = MODEL
    rw, B = rate_width, 15
    # horizontal edges flip Y (rate depends on X = the row); gain bows OUTSIDE, loss inside
    curved_arrow(d, P["00"], P["01"], +1, B, rw(m["y_gain_x0"]), f'{m["y_gain_x0"]:g}')
    curved_arrow(d, P["01"], P["00"], -1, B, rw(m["y_loss_x0"]), f'{m["y_loss_x0"]:g}')
    curved_arrow(d, P["10"], P["11"], -1, B, rw(m["y_gain_x1"]), f'{m["y_gain_x1"]:g}')
    curved_arrow(d, P["11"], P["10"], +1, B, rw(m["y_loss_x1"]), f'{m["y_loss_x1"]:g}')
    # vertical edges flip X (rate depends on Y = the column)
    curved_arrow(d, P["00"], P["10"], +1, B, rw(m["x_gain_y0"]), f'{m["x_gain_y0"]:g}')
    curved_arrow(d, P["10"], P["00"], -1, B, rw(m["x_loss_y0"]), f'{m["x_loss_y0"]:g}')
    curved_arrow(d, P["01"], P["11"], -1, B, rw(m["x_gain_y1"]), f'{m["x_gain_y1"]:g}')
    curved_arrow(d, P["11"], P["01"], +1, B, rw(m["x_loss_y1"]), f'{m["x_loss_y1"]:g}')
    for lab, (x, y) in P.items():
        state_node(d, x, y, lab)
    # which bit each edge flips (kept well clear of the rate labels)
    d.append(draw.Text("flip Y", FS_TICK, cx0 + gx / 2, cy0 - 92, font_family=FONT,
                       text_anchor="middle", fill=MUTED, font_style="italic"))
    d.append(draw.Text("flip X", FS_TICK, cx0 - 98, cy0 + gy / 2, font_family=FONT,
                       text_anchor="middle", dominant_baseline="central",
                       fill=MUTED, font_style="italic"))
    d.append(draw.Text("states = (X, Y);  arrow width = rate", FS_TICK, cx0 + gx / 2,
                       cy0 + gy + 96, font_family=FONT, text_anchor="middle", fill=MUTED))


# --------------------------------------------------------------------------- panel B: a realization
def _layout(tree):
    tfo, order = {}, []
    for n in tree.traverse("preorder"):
        tfo[n.name] = 0.0 if n.is_root() else tfo[n.up.name] + n.dist
    present = max(tfo.values())
    ys, i = {}, 0
    for n in tree.traverse("postorder"):
        if n.is_leaf():
            ys[n.name] = i
            i += 1
        else:
            ys[n.name] = sum(ys[c.name] for c in n.children) / len(n.children)
    return tfo, present, ys, i


def panel_realization(d, ox, oy, pw, ph):
    ztree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N_TIPS, age=AGE,
                                  direction="backward", seed=TREE_SEED)
    model = CorrelatedBinary(**MODEL)
    # pick a trait seed where the coupling is visibly expressed (many tips with X == Y)
    best = None
    for s in range(1, 60):
        res = simulate_traits(ztree, model, seed=s, root_state=(0, 0)) \
            if "root_state" in simulate_traits.__code__.co_varnames else \
            simulate_traits(ztree, model, seed=s)
        tips = res.labeled_values().values()
        agree = sum(1 for (x, y) in tips) and sum(1 for (x, y) in tips if x == y)
        xs = {x for (x, y) in tips}
        if len(xs) == 2 and (best is None or agree > best[0]):
            best = (agree, s, res)
    _, seed, res = best
    labeled = {n.name: xy for n, xy in res.labeled_values().items()}
    xstate = {n.name: model.states[i][0] for n, i in res.node_values.items()}
    xchg = {}
    for node, t, frm, to in res.changes():
        if frm[0] != to[0]:
            xchg.setdefault(node.name, []).append((t, to[0]))

    tree = zombi_to_ete3(ztree)
    tfo, present, ys, nleaf = _layout(tree)
    x_at = lambda t: ox + 60 + (t / present) * (pw - 220)      # noqa: E731
    y_at = lambda k: oy + 40 + (k / max(1, nleaf - 1)) * (ph - 90)

    # P#1: panel letter top-left; title centred over the panel's content (tree + chips)
    d.append(draw.Text("B", FS_LABEL, ox, oy - 24, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))
    d.append(draw.Text("a simulated realization", FS_LABEL, ox + pw / 2, oy - 24,
                       font_family=FONT, text_anchor="middle", fill=INK, font_weight="bold"))

    def seg(x1, x2, y, on):
        d.append(draw.Line(x1, y, x2, y, stroke=PRESENT if on else GREY,
                           stroke_width=5.2 if on else 2.4, stroke_linecap="butt"))

    # branches, painted by the X stochastic map
    for n in tree.traverse():
        if n.is_root():
            continue
        y = y_at(ys[n.name])
        t0, cur = tfo[n.up.name], xstate[n.up.name]
        for tt, to in sorted(xchg.get(n.name, [])):
            seg(x_at(t0), x_at(tt), y, cur)
            t0, cur = tt, to
        seg(x_at(t0), x_at(tfo[n.name]), y, cur)
    # vertical connectors
    for n in tree.traverse("postorder"):
        if not n.is_leaf():
            x = x_at(tfo[n.name])
            yy = [y_at(ys[c.name]) for c in n.children]
            d.append(draw.Line(x, min(yy), x, max(yy), stroke=INK, stroke_width=2.4))

    # tip chips: two columns X, Y
    colx = ox + pw - 120
    coly = ox + pw - 74
    d.append(draw.Text("X", FS_TICK, colx, oy + 14, font_family=FONT, text_anchor="middle",
                       fill=INK, font_weight="bold"))
    d.append(draw.Text("Y", FS_TICK, coly, oy + 14, font_family=FONT, text_anchor="middle",
                       fill=INK, font_weight="bold"))
    for n in tree.get_leaves():
        y = y_at(ys[n.name])
        x, yv = labeled[n.name]
        chip(d, colx, y, x)
        chip(d, coly, y, yv)

    # time axis
    base = oy + ph - 14
    d.append(draw.Line(x_at(0), base, x_at(present), base, stroke=INK, stroke_width=1.6))
    for k in range(5):
        tv = present * k / 4
        xx = x_at(tv)
        d.append(draw.Line(xx, base, xx, base + 6, stroke=INK, stroke_width=1.6))
        d.append(draw.Text(f"{tv:.2f}", FS_TICK, xx, base + 22, font_family=FONT,
                           text_anchor="middle", fill=INK))
    d.append(draw.Text("time (root to present)", FS_TICK, (x_at(0) + x_at(present)) / 2,
                       base + 44, font_family=FONT, text_anchor="middle", fill=MUTED))
    return seed


# --------------------------------------------------------------------------- render
def render_one(name):
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Correlated binary characters (Pagel)", FS_TITLE, W / 2, 46,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))
    # legend (clear of the data), centered under the title
    ly = 80
    # L#1: swatch middle y = (ly-12)+22/2 = ly-1; centre the label text on it
    ty = (ly - 1) + 0.34 * FS_TICK
    d.append(draw.Rectangle(W / 2 - 165, ly - 12, 22, 22, fill=PRESENT, stroke=INK, stroke_width=1.4))
    d.append(draw.Text("state 1 (present)", FS_TICK, W / 2 - 135, ty, font_family=FONT,
                       text_anchor="start", fill=INK))
    # second swatch pushed right so "(present)" never touches the white square
    d.append(draw.Rectangle(W / 2 + 78, ly - 12, 22, 22, fill="white", stroke=INK, stroke_width=1.4))
    d.append(draw.Text("state 0 (absent)", FS_TICK, W / 2 + 108, ty, font_family=FONT,
                       text_anchor="start", fill=INK))

    panel_model(d, 212, 272)                       # (cx0, cy0) = top-left state centre
    seed = panel_realization(d, 520, 150, 620, 470)

    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png  (trait seed {seed})")


def render():
    """Render both palettes: colour -> trait_pagel, B&W -> trait_pagel_bw."""
    global PRESENT
    PRESENT = PRESENT_COLOR
    render_one("trait_pagel")
    PRESENT = PRESENT_BW
    render_one("trait_pagel_bw")


if __name__ == "__main__":
    render()
