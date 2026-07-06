"""Figure (Ch15): the core transformation -- chronogram into phylogram.

The single idea sequence evolution rests on. The SAME tree, drawn twice.

  * Left (chronogram): branch lengths are TIME. Every tip lines up at the present -- the
    tree is ultrametric -- because every lineage has had the same amount of time to evolve.
  * Right (phylogram): branch lengths are SUBSTITUTIONS per site. Each branch has been
    multiplied by an evolutionary rate (a relaxed clock), so fast branches stretch and slow
    branches shrink. The tips no longer line up: what a sequence records is substitutions,
    not time.

Both trees are painted by the SAME per-branch rate, so you can watch each branch move: the
yellow (fast) branches that were short in time become long in substitutions, the purple
(slow) ones collapse. This is what turns a timetree into something you could infer from an
alignment.

House style: painted trees, viridis rate scale, one centered title, ASCII text.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_seq_chrono_phylo.py
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

W, H = 1320, 660
N_TIPS, TREE_SEED = 16, 4
SIGMA, CLOCK_SEED = 0.6, 5
BRANCH_W = 3.4


def draw_panel(d, ox, oy, pw, ph, tree, ys, nleaf, rate, dist, present, title, tag, unit, ragged):
    x_at = lambda v: ox + (v / present) * pw              # noqa: E731
    y_at = lambda k: oy + (k / max(1, nleaf - 1)) * ph    # noqa: E731

    d.append(draw.Text(title, FS_LABEL, ox + pw / 2, oy - 44, font_family=FONT,
                       text_anchor="middle", font_weight="bold", fill=INK))
    d.append(draw.Text(tag, FS_TICK, ox + pw / 2, oy - 22, font_family=FONT,
                       text_anchor="middle", fill=MUTED, font_style="italic"))
    for n in tree.traverse():
        if n.is_root():
            continue
        y = y_at(ys[n.name])
        d.append(draw.Line(x_at(dist[n.up.name]), y, x_at(dist[n.name]), y,
                           stroke=C.rate_hex(rate[n.name]), stroke_width=BRANCH_W, stroke_linecap="round"))
    for n in tree.traverse("postorder"):
        if not n.is_leaf():
            x = x_at(dist[n.name])
            yy = [y_at(ys[c.name]) for c in n.children]
            d.append(draw.Line(x, min(yy), x, max(yy), stroke=INK, stroke_width=1.4))
    for lf in tree.get_leaves():
        d.append(draw.Circle(x_at(dist[lf.name]), y_at(ys[lf.name]), 3.0,
                             fill=C.rate_hex(rate[lf.name]), stroke="white", stroke_width=0.6))

    # tips-aligned reference (chronogram) vs "ragged" note (phylogram)
    if not ragged:
        xp = x_at(max(dist[l.name] for l in tree.get_leaves()))
        d.append(draw.Line(xp, y_at(0) - 8, xp, y_at(nleaf - 1) + 8, stroke=MUTED,
                           stroke_width=1.2, stroke_dasharray="4,4"))
        d.append(draw.Text("tips aligned", FS_TICK, xp, y_at(0) - 16, font_family=FONT,
                           text_anchor="middle", fill=MUTED, font_style="italic"))
    else:
        d.append(draw.Text("tips no longer aligned", FS_TICK, ox + pw, oy + 6,
                           font_family=FONT, text_anchor="end", fill=MUTED, font_style="italic"))

    # axis
    base = oy + ph + 34
    d.append(draw.Line(x_at(0), base, ox + pw, base, stroke=INK, stroke_width=1.5))
    for k in range(3):
        v = present * k / 2
        xx = x_at(v)
        d.append(draw.Line(xx, base, xx, base + 5, stroke=INK, stroke_width=1.5))
        d.append(draw.Text(f"{v:.1f}", FS_TICK, xx, base + 20, font_family=FONT,
                           text_anchor="middle", fill=MUTED))
    d.append(draw.Text(unit, FS_TICK, ox + pw / 2, base + 42, font_family=FONT,
                       text_anchor="middle", fill=INK, font_style="italic"))


def main():
    tree = C.build_tree(n_tips=N_TIPS, seed=TREE_SEED)
    tfo, present = C.node_times(tree)
    ys, nleaf = C.leaf_ys(tree)
    rate, sub = C.autocorrelated_lognormal(tree, SIGMA, seed=CLOCK_SEED)
    dist = C.subst_dist_to_root(tree, sub)
    subpresent = max(dist[l.name] for l in tree.get_leaves())

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("From time to substitutions: chronogram into phylogram", FS_TITLE,
                       W / 2, 44, font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))
    _rate_bar(d, W / 2 - 165, 76, 330, 15)

    pw, ph = 470, 400
    oyt = 178
    lx, rx = 70, W - 70 - pw
    draw_panel(d, lx, oyt, pw, ph, tree, ys, nleaf, rate,
               {k: v for k, v in tfo.items()}, present,
               "Chronogram", "branch length = time", "time", ragged=False)
    draw_panel(d, rx, oyt, pw, ph, tree, ys, nleaf, rate, dist, subpresent,
               "Phylogram", "branch length = substitutions", "substitutions / site", ragged=True)

    # transformation arrow between the panels
    ay = oyt + ph / 2
    ax0, ax1 = lx + pw + 20, rx - 20
    d.append(draw.Line(ax0, ay, ax1 - 10, ay, stroke=INK, stroke_width=3.0))
    d.append(draw.Lines(ax1, ay, ax1 - 16, ay - 9, ax1 - 16, ay + 9, close=True, fill=INK))
    d.append(draw.Text("apply a", FS_TICK, (ax0 + ax1) / 2, ay - 40, font_family=FONT,
                       text_anchor="middle", fill=INK))
    d.append(draw.Text("relaxed clock", FS_ANNOT, (ax0 + ax1) / 2, ay - 20, font_family=FONT,
                       text_anchor="middle", fill=INK, font_weight="bold"))
    d.append(draw.Text("multiply each", FS_TICK, (ax0 + ax1) / 2, ay + 24, font_family=FONT,
                       text_anchor="middle", fill=MUTED))
    d.append(draw.Text("branch by its rate", FS_TICK, (ax0 + ax1) / 2, ay + 42, font_family=FONT,
                       text_anchor="middle", fill=MUTED))

    name = "seq_chrono_phylo"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png  (subst depth {subpresent:.2f})")


def _rate_bar(d, x, y, w, h):
    grad = draw.LinearGradient(x, y, x + w, y)
    for i in range(21):
        t = i / 20.0
        grad.add_stop(t, C.rate_hex(C.RATE_LO * (C.RATE_HI / C.RATE_LO) ** t))
    d.append(grad)
    d.append(draw.Text("branch rate", FS_TICK, x + w / 2, y - 8, font_family=FONT,
                       text_anchor="middle", font_weight="bold", fill=INK))
    d.append(draw.Rectangle(x, y, w, h, fill=grad, stroke=INK, stroke_width=0.8))
    for tx, lab, anc in ((x, "slow", "start"), (x + w / 2, "1x", "middle"), (x + w, "fast", "end")):
        d.append(draw.Text(lab, FS_TICK, tx, y + h + 16, font_family=FONT, text_anchor=anc, fill=MUTED))


if __name__ == "__main__":
    main()
