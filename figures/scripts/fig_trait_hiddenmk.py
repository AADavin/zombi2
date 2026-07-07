"""Figure: hidden rate classes (corHMM / HiddenStateMk) — two panels.

An observed binary character 0/1 evolves along the tree, but its transition rate depends on an
unobserved *hidden class* that itself switches between a slow and a fast regime. That is rate
heterogeneity a plain Mk cannot represent.

  * Panel A (the model): the four joint states (observed 0/1) x (hidden slow/fast). Horizontal
    edges flip the observed character — thin in the slow row, thick in the fast row (arrow width
    = rate) — and vertical edges switch the hidden class. So the SAME character flips fast or slow
    depending on the hidden regime.
  * Panel B (a realization): branches drawn heavy where the hidden class is fast and light where it
    is slow; small circles mark observed-state changes; tip chips give the observed 0/1. The
    observed changes cluster on the fast (heavy) branches — the signature of hidden rate classes.

This figure is produced in two palettes from one run: a soft, low-saturation COLOUR version
(the default, trait_hiddenmk.svg/.png) and the original B&W version (trait_hiddenmk_bw.svg/.png).
In colour the fast hidden class (and the observed-1 chips) is a muted accent hue; in B&W it is
near-black. The BW look is preserved so the monochrome figure can always be regenerated.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_trait_hiddenmk.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw

from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.traits import HiddenStateMk, simulate_traits

from fig_trait_pagel import curved_arrow, rate_width, state_node, chip, _layout
from model_common import zombi_to_ete3
from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1200, 660
GREY = "#9a9a9a"
# The fast hidden class (heavy branches) and the observed-1 chips carry the colour. Two palettes;
# BW preserves the original near-black look, COLOUR uses the shared dusty-blue accent. FAST_INK is
# set per render mode.
FAST_BW    = INK
FAST_COLOR = "#5a86a0"               # dusty blue, matches the Mk / Pagel colour set
FAST_INK = FAST_COLOR               # active "fast / present" colour (set per mode)
OBS_STATES, HID_STATES = ["0", "1"], ["slow", "fast"]
SLOW, FAST, HIDDEN_RATE = 0.4, 3.5, 0.7
N_TIPS, AGE, TREE_SEED, TRAIT_SEED = 11, 1.0, 3, 9


# --------------------------------------------------------------------------- panel A: the model
def panel_model(d, cx0, cy0):
    gx = gy = 176
    # P#1: panel letter top-left; title horizontally centred over the panel
    d.append(draw.Text("A", FS_LABEL, cx0 - 40, cy0 - 146, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))
    d.append(draw.Text("the model", FS_LABEL, cx0 + gx / 2, cy0 - 146, font_family=FONT,
                       text_anchor="middle", fill=INK, font_weight="bold"))
    P = {"0s": (cx0, cy0), "1s": (cx0 + gx, cy0),
         "0f": (cx0, cy0 + gy), "1f": (cx0 + gx, cy0 + gy)}
    B = 15
    # observed flips (horizontal): slow row thin, fast row thick
    curved_arrow(d, P["0s"], P["1s"], +1, B, rate_width(SLOW), f"{SLOW:g}")
    curved_arrow(d, P["1s"], P["0s"], -1, B, rate_width(SLOW), f"{SLOW:g}")
    curved_arrow(d, P["0f"], P["1f"], -1, B, rate_width(FAST), f"{FAST:g}")
    curved_arrow(d, P["1f"], P["0f"], +1, B, rate_width(FAST), f"{FAST:g}")
    # hidden-class switches (vertical)
    curved_arrow(d, P["0s"], P["0f"], +1, B, rate_width(HIDDEN_RATE), f"{HIDDEN_RATE:g}")
    curved_arrow(d, P["0f"], P["0s"], -1, B, rate_width(HIDDEN_RATE), "")
    curved_arrow(d, P["1s"], P["1f"], -1, B, rate_width(HIDDEN_RATE), f"{HIDDEN_RATE:g}")
    curved_arrow(d, P["1f"], P["1s"], +1, B, rate_width(HIDDEN_RATE), "")
    for key, (x, y) in P.items():
        state_node(d, x, y, key[0])                        # the observed digit
    # row / column hints
    d.append(draw.Text("slow", FS_TICK, cx0 - 92, cy0, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill=MUTED, font_style="italic"))
    d.append(draw.Text("fast", FS_TICK, cx0 - 92, cy0 + gy, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill=MUTED, font_style="italic"))
    d.append(draw.Text("flip observed", FS_TICK, cx0 + gx / 2, cy0 - 92, font_family=FONT,
                       text_anchor="middle", fill=MUTED, font_style="italic"))
    d.append(draw.Text("nodes = observed 0/1;  arrow width = rate", FS_TICK, cx0 + gx / 2,
                       cy0 + gy + 96, font_family=FONT, text_anchor="middle", fill=MUTED))


# --------------------------------------------------------------------------- panel B: a realization
def panel_realization(d, ox, oy, pw, ph):
    ztree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N_TIPS, age=AGE,
                                  direction="backward", seed=TREE_SEED)
    Qslow = [[-SLOW, SLOW], [SLOW, -SLOW]]
    Qfast = [[-FAST, FAST], [FAST, -FAST]]
    model = HiddenStateMk(observed_rates=[Qslow, Qfast], hidden_rate=HIDDEN_RATE,
                          observed_states=OBS_STATES, hidden_states=HID_STATES)
    res = simulate_traits(ztree, model, seed=TRAIT_SEED)
    full = {n.name: res.full_label(i) for n, i in res.node_values.items()}   # (obs, hid)
    obs_tip = {n.name: full[n.name][0] for n in res.node_values if n.is_leaf()}
    hid_chg, obs_chg = {}, {}
    for node, t, frm, to in res.changes():
        if frm[1] != to[1]:
            hid_chg.setdefault(node.name, []).append((t, to[1]))
        if frm[0] != to[0]:
            obs_chg.setdefault(node.name, []).append(t)

    tree = zombi_to_ete3(ztree)
    tfo, present, ys, nleaf = _layout(tree)
    x_at = lambda t: ox + 60 + (t / present) * (pw - 210)      # noqa: E731
    y_at = lambda k: oy + 40 + (k / max(1, nleaf - 1)) * (ph - 90)

    # P#1: panel letter top-left; title centred over the panel's content (tree + chips)
    d.append(draw.Text("B", FS_LABEL, ox, oy - 24, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))
    d.append(draw.Text("a simulated realization", FS_LABEL, ox + pw / 2, oy - 24,
                       font_family=FONT, text_anchor="middle", fill=INK, font_weight="bold"))

    def seg(x1, x2, y, fast):
        d.append(draw.Line(x1, y, x2, y, stroke=FAST_INK if fast else GREY,
                           stroke_width=5.2 if fast else 2.4, stroke_linecap="butt"))

    for n in tree.traverse():
        if n.is_root():
            continue
        y = y_at(ys[n.name])
        t0, cur = tfo[n.up.name], full[n.up.name][1]
        for tt, to in sorted(hid_chg.get(n.name, [])):
            seg(x_at(t0), x_at(tt), y, cur == "fast")
            t0, cur = tt, to
        seg(x_at(t0), x_at(tfo[n.name]), y, cur == "fast")
        for tt in obs_chg.get(n.name, []):                 # mark each observed-state change
            d.append(draw.Circle(x_at(tt), y, 4.6, fill="white", stroke=INK, stroke_width=1.8))
    for n in tree.traverse("postorder"):
        if not n.is_leaf():
            x = x_at(tfo[n.name])
            yy = [y_at(ys[c.name]) for c in n.children]
            d.append(draw.Line(x, min(yy), x, max(yy), stroke=INK, stroke_width=2.4))

    colc = ox + pw - 96
    d.append(draw.Text("obs", FS_TICK, colc, oy + 14, font_family=FONT, text_anchor="middle",
                       fill=INK, font_weight="bold"))
    for n in tree.get_leaves():
        y = y_at(ys[n.name])
        chip(d, colc, y, obs_tip[n.name] == "1", fill_on=FAST_INK)

    base = oy + ph - 14
    d.append(draw.Line(x_at(0), base, x_at(present), base, stroke=INK, stroke_width=1.6))
    for k in range(5):
        tv = present * k / 4
        xx = x_at(tv)
        d.append(draw.Line(xx, base, xx, base + 6, stroke=INK, stroke_width=1.6))
        d.append(draw.Text(f"{tv:.2f}", FS_TICK, xx, base + 22, font_family=FONT,
                           text_anchor="middle", fill=INK))
    d.append(draw.Text("time (root to present)", FS_TICK, (x_at(0) + x_at(present)) / 2, base + 44,
                       font_family=FONT, text_anchor="middle", fill=MUTED))


# --------------------------------------------------------------------------- render
def render_one(name):
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Hidden rate classes (corHMM)", FS_TITLE, W / 2, 46, font_family=FONT,
                       text_anchor="middle", font_weight="bold", fill=INK))
    # legend: hidden-class branch shading + observed chips.
    # Laid out left-to-right with generous gaps (no overlaps) and centred as a block;
    # L#1: text baseline set explicitly so each label is vertically centred on its marker.
    ly = 82
    ty = ly + 0.34 * FS_TICK
    MARK, GAPT, GAPE = 32, 12, 46      # marker length, marker->text gap, entry->entry gap
    entries = [("hidden: fast", 118, "fastline"),
               ("hidden: slow", 126, "slowline"),
               ("observed change", 173, "circle")]
    total = sum(MARK + GAPT + w for _, w, _ in entries) + GAPE * (len(entries) - 1)
    x = W / 2 - total / 2
    for label, tw, kind in entries:
        if kind == "fastline":
            d.append(draw.Line(x, ly, x + MARK, ly, stroke=FAST_INK, stroke_width=5.2))
        elif kind == "slowline":
            d.append(draw.Line(x, ly, x + MARK, ly, stroke=GREY, stroke_width=2.4))
        else:
            d.append(draw.Circle(x + MARK / 2, ly, 4.6, fill="white", stroke=INK, stroke_width=1.8))
        d.append(draw.Text(label, FS_TICK, x + MARK + GAPT, ty, font_family=FONT,
                           text_anchor="start", fill=INK))
        x += MARK + GAPT + tw + GAPE

    panel_model(d, 212, 272)
    panel_realization(d, 520, 150, 620, 470)

    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png")


def render():
    """Render both palettes: colour -> trait_hiddenmk, B&W -> trait_hiddenmk_bw."""
    global FAST_INK
    FAST_INK = FAST_COLOR
    render_one("trait_hiddenmk")
    FAST_INK = FAST_BW
    render_one("trait_hiddenmk_bw")


if __name__ == "__main__":
    render()
