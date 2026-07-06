"""Figure (Ch16): the whole relaxed-clock family rescaling ONE tree.

The flagship comparison. A single time-calibrated tree (top-left, black, in TIME) is
rescaled into substitutions by every clock in the family. Each small tree is a
*phylogram*: branch lengths are expected substitutions per site, and each branch is
painted by the rate multiplier the clock drew for it (viridis, a shared log scale, so a
colour means the same rate in every panel). Every panel shares the SAME horizontal
substitutions-per-pixel scale, so branch lengths are directly comparable across models.

Read it as the chapter's summary table made visual: strict leaves the tree undistorted;
the uncorrelated clocks scatter branch rates independently (tips no longer line up); the
autocorrelated clocks vary the rate smoothly down the tree.

House style: painted trees, one centered title, ASCII text, shared viridis rate bar.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_clock_family.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw

import clock_common as C
from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

# --- panel spec: (key, title, tag, draw-fn) ---------------------------------
# each draw-fn(tree) -> (rate, subst_len); chronogram is special-cased.
SEED = 7
PANELS = [
    ("chrono", "Chronogram", "input (time)", None),
    ("strict", "Strict", "one rate", lambda t: C.strict(t)),
    ("uln", "Uncorrelated lognormal", "uncorrelated", lambda t: C.uncorrelated_lognormal(t, 0.55, SEED)),
    ("ugam", "Uncorrelated gamma", "uncorrelated", lambda t: C.uncorrelated_gamma(t, 2.5, SEED)),
    ("wn", "White noise", "uncorrelated", lambda t: C.white_noise(t, 0.32, SEED)),
    ("aln", "Autocorrelated lognormal", "autocorrelated", lambda t: C.autocorrelated_lognormal(t, 0.5, SEED + 2)),
    ("cir", "Cox-Ingersoll-Ross", "autocorrelated", lambda t: C.cir(t, theta=2.0, sigma=0.6, mean=1.0, seed=SEED)[:2]),
    ("bin", "Discrete-bin (GTDB)", "autocorrelated", lambda t: C.discrete_bin(t, switch_rate=2.0, seed=SEED + 1)[:2]),
]

NCOL, NROW = 4, 2
W, H = 1520, 1030
MARGIN_X = 40
TOP = 150                       # title + rate bar band
ROW_GAP = 34
PANEL_W = (W - 2 * MARGIN_X) / NCOL
PANEL_H = (H - TOP - 24 - ROW_GAP) / NROW
PAD_L, PAD_R, PAD_T = 20, 26, 70              # inside a panel cell
INNER_W = PANEL_W - PAD_L - PAD_R
INNER_H = PANEL_H - PAD_T - 92                # 92 px reserved below the tree for the axis
BRANCH_W = 3.4


def main():
    tree = C.build_tree()
    tfo, present = C.node_times(tree)
    ys, nleaf = C.leaf_ys(tree)

    # compute every panel's phylogram, then a SHARED subs-per-pixel scale
    data = {}
    max_depth = present            # chronogram depth
    for key, _, _, fn in PANELS:
        if key == "chrono":
            data[key] = (None, tfo)          # colour=None (black), x = time
            continue
        rate, sub = fn(tree)
        dr = C.subst_dist_to_root(tree, sub)
        data[key] = (rate, dr)
        max_depth = max(max_depth, max(dr[l.name] for l in tree.get_leaves()))
    import math as _m
    axis_max = _m.ceil(max_depth / 0.5) * 0.5     # round up to a nice tick
    scale = INNER_W / axis_max                    # px per substitution (== per time unit)

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("One tree, every clock: time rescaled into substitutions", FS_TITLE,
                       W / 2, 50, font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))

    # shared viridis rate bar, centered under the title
    _rate_bar(d, W / 2 - 190, 92, 380, 18)

    for idx, (key, title, tag, _) in enumerate(PANELS):
        r, c = divmod(idx, NCOL)
        ox = MARGIN_X + c * PANEL_W
        oy = TOP + r * (PANEL_H + ROW_GAP)
        _panel(d, ox, oy, tree, ys, nleaf, data[key][0], data[key][1], scale, axis_max,
               title, tag, is_chrono=(key == "chrono"))

    name = "clock_family"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png  (shared scale: {scale:.1f}px per subst; max depth {max_depth:.2f})")


def _panel(d, ox, oy, tree, ys, nleaf, rate, dist, scale, axis_max, title, tag, is_chrono):
    x_at = lambda v: ox + PAD_L + v * scale                       # noqa: E731
    y_at = lambda k: oy + PAD_T + (k / max(1, nleaf - 1)) * INNER_H  # noqa: E731

    # panel title + tag
    d.append(draw.Text(title, FS_LABEL, ox + PANEL_W / 2, oy + 22, font_family=FONT,
                       text_anchor="middle", font_weight="bold", fill=INK))
    d.append(draw.Text(tag, FS_TICK, ox + PANEL_W / 2, oy + 44, font_family=FONT,
                       text_anchor="middle", fill=MUTED, font_style="italic"))

    # branches, each coloured by its rate (chronogram = solid black)
    for n in tree.traverse():
        if n.is_root():
            continue
        y = y_at(ys[n.name])
        x0, x1 = x_at(dist[n.up.name]), x_at(dist[n.name])
        col = INK if is_chrono else C.rate_hex(rate[n.name])
        d.append(draw.Line(x0, y, x1, y, stroke=col, stroke_width=BRANCH_W, stroke_linecap="round"))
    # vertical connectors (thin, neutral)
    for n in tree.traverse("postorder"):
        if not n.is_leaf():
            x = x_at(dist[n.name])
            yy = [y_at(ys[ch.name]) for ch in n.children]
            d.append(draw.Line(x, min(yy), x, max(yy), stroke=INK, stroke_width=1.5))
    # tip ticks (short marks so non-ultrametric tips are legible)
    for lf in tree.get_leaves():
        y = y_at(ys[lf.name])
        x = x_at(dist[lf.name])
        d.append(draw.Circle(x, y, 2.4, fill=INK if is_chrono else C.rate_hex(rate[lf.name]),
                             stroke="white", stroke_width=0.5))

    # a light "present" reference line for the chronogram (tips aligned) so the
    # contrast with the ragged phylogram tips is explicit
    if is_chrono:
        xp = x_at(max(dist[l.name] for l in tree.get_leaves()))
        d.append(draw.Line(xp, y_at(0) - 6, xp, y_at(nleaf - 1) + 6, stroke=MUTED,
                           stroke_width=1.0, stroke_dasharray="3,3"))

    # per-panel axis baseline: time for the chronogram, substitutions otherwise
    base = oy + PAD_T + INNER_H + 30
    x0 = x_at(0.0)
    xmax = ox + PAD_L + INNER_W
    d.append(draw.Line(x0, base, xmax, base, stroke=INK, stroke_width=1.4))
    unit = "time" if is_chrono else "subs/site"
    for k in range(3):
        v = axis_max * k / 2
        xx = x_at(v)
        d.append(draw.Line(xx, base, xx, base + 5, stroke=INK, stroke_width=1.4))
        d.append(draw.Text(f"{v:g}", FS_TICK, xx, base + 20, font_family=FONT,
                           text_anchor="middle", fill=MUTED))
    d.append(draw.Text(unit, FS_TICK, (x0 + xmax) / 2, base + 40, font_family=FONT,
                       text_anchor="middle", fill=MUTED, font_style="italic"))


def _rate_bar(d, x, y, w, h):
    grad = draw.LinearGradient(x, y, x + w, y)
    for i in range(21):
        t = i / 20.0
        r = C.RATE_LO * (C.RATE_HI / C.RATE_LO) ** t     # log-spaced
        grad.add_stop(t, C.rate_hex(r))
    d.append(grad)
    d.append(draw.Text("branch rate multiplier", FS_LABEL, x + w / 2, y - 12, font_family=FONT,
                       text_anchor="middle", font_weight="bold", fill=INK))
    d.append(draw.Rectangle(x, y, w, h, fill=grad, stroke=INK, stroke_width=0.8))
    for tx, lab, anc in ((x, "0.33x (slow)", "start"),
                         (x + w / 2, "1x", "middle"),
                         (x + w, "3x (fast)", "end")):
        d.append(draw.Text(lab, FS_TICK, tx, y + h + 18, font_family=FONT, text_anchor=anc, fill=MUTED))


if __name__ == "__main__":
    main()
