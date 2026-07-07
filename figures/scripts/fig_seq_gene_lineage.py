"""Figure (Ch15): the gene x lineage clock, rate = R_b * s_g.

Sequence evolution factorises a branch's substitution rate into two independent pieces:

    rate(family g, species branch b)  =  R_b  x  s_g

  * R_b -- the shared lineage clock. One rate is drawn per species branch and SHARED by
    every gene family, so if a clade evolves fast, it is fast for all its genes. (Left: the
    species tree painted by R_b.)
  * s_g -- the per-family speed. Each family draws one constant, a global multiplier on all
    its branches, so some genes are intrinsically fast and others slow.

The three gene-family phylograms (right) are the same tree scaled by three family speeds.
They share R_b, so their rate PATTERN -- the colours, which clades are fast or slow -- is
identical; the family speed only sets each tree's overall SIZE. Length carries s_g, colour
carries R_b.

House style: painted trees, viridis rate scale, ASCII text (R_b, s_g written out).

Run:  python figures/scripts/fig_seq_gene_lineage.py
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

W, H = 1340, 780
N_TIPS, TREE_SEED = 12, 3
SIGMA, CLOCK_SEED = 0.6, 8
FAMILIES = [("g1", 0.6, "slow-evolving gene"),
            ("g2", 1.0, "average gene"),
            ("g3", 1.7, "fast-evolving gene")]


def paint_tree(d, ox, oy, pw, ph, tree, ys, nleaf, rate, dist, maxd, bw, tips=True):
    x_at = lambda v: ox + (v / maxd) * pw                 # noqa: E731
    y_at = lambda k: oy + (k / max(1, nleaf - 1)) * ph    # noqa: E731
    for n in tree.traverse():
        if n.is_root():
            continue
        y = y_at(ys[n.name])
        d.append(draw.Line(x_at(dist[n.up.name]), y, x_at(dist[n.name]), y,
                           stroke=C.rate_hex(rate[n.name]), stroke_width=bw, stroke_linecap="round"))
    for n in tree.traverse("postorder"):
        if not n.is_leaf():
            x = x_at(dist[n.name])
            yy = [y_at(ys[c.name]) for c in n.children]
            d.append(draw.Line(x, min(yy), x, max(yy), stroke=INK, stroke_width=1.3))
    if tips:
        for lf in tree.get_leaves():
            d.append(draw.Circle(x_at(dist[lf.name]), y_at(ys[lf.name]), 2.6,
                                 fill=C.rate_hex(rate[lf.name]), stroke="white", stroke_width=0.5))
    return x_at


def main():
    tree = C.build_tree(n_tips=N_TIPS, seed=TREE_SEED)
    tfo, present = C.node_times(tree)
    ys, nleaf = C.leaf_ys(tree)
    Rb, Sbase = C.autocorrelated_lognormal(tree, SIGMA, seed=CLOCK_SEED)   # shared lineage clock
    base_dist = C.subst_dist_to_root(tree, Sbase)
    base_depth = max(base_dist[l.name] for l in tree.get_leaves())
    smax = max(s for _, s, _ in FAMILIES)
    shared_maxd = base_depth * smax * 1.02

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("The gene x lineage clock: one shared pattern, one speed per family",
                       FS_TITLE, W / 2, 44, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))

    # --- header equation ----------------------------------------------------
    ey = 96
    cx = W / 2
    d.append(draw.Text("branch substitution rate", FS_LABEL, cx - 250, ey, font_family=FONT,
                       text_anchor="middle", fill=INK))
    d.append(draw.Text("=", FS_TITLE, cx - 90, ey + 2, font_family=FONT, text_anchor="middle", fill=INK))
    d.append(draw.Text("R_b", FS_TITLE, cx - 30, ey + 3, font_family=FONT, text_anchor="middle",
                       fill=C.rate_hex(1.9), font_weight="bold"))
    d.append(draw.Text("x", FS_LABEL, cx + 80, ey, font_family=FONT, text_anchor="middle", fill=INK))
    d.append(draw.Text("s_g", FS_TITLE, cx + 190, ey + 3, font_family=FONT, text_anchor="middle",
                       fill=INK, font_weight="bold"))
    d.append(draw.Text("one per tree (shared)", FS_TICK, cx - 30, ey + 26,
                       font_family=FONT, text_anchor="middle", fill=MUTED))
    d.append(draw.Text("one per family", FS_TICK, cx + 190, ey + 26,
                       font_family=FONT, text_anchor="middle", fill=MUTED))

    # --- left: species tree painted by R_b ----------------------------------
    lx, loy, lpw, lph = 70, 210, 430, 430
    d.append(draw.Text("Shared lineage clock  R_b", FS_LABEL, lx + lpw / 2, loy - 26,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))
    d.append(draw.Text("drawn once on the species tree", FS_TICK, lx + lpw / 2, loy - 6,
                       font_family=FONT, text_anchor="middle", fill=MUTED, font_style="italic"))
    paint_tree(d, lx, loy + 14, lpw, lph - 40, tree, ys, nleaf, Rb, base_dist, base_depth, 4.4)
    _rate_bar(d, lx + 30, loy + lph + 6, lpw - 120, 15)

    # --- right: three family phylograms, shared subst scale ------------------
    rx, rpw = 620, 640
    rows_top, row_h = 200, 150
    for i, (gid, sg, desc) in enumerate(FAMILIES):
        roy = rows_top + i * row_h
        fam_dist = {k: v * sg for k, v in base_dist.items()}
        d.append(draw.Text(f"Family {gid}:  s_g = {sg:g}", FS_LABEL, rx, roy - 20, font_family=FONT,
                           text_anchor="start", font_weight="bold", fill=INK))
        d.append(draw.Text(f"({desc})", FS_TICK, rx + 260, roy - 20, font_family=FONT,
                           text_anchor="start", fill=MUTED, font_style="italic"))
        x_at = paint_tree(d, rx, roy, rpw, row_h - 62, tree, ys, nleaf, Rb, fam_dist,
                          shared_maxd, 3.0)
        # size bracket: show how far this family's deepest tip reaches
        tipx = x_at(max(fam_dist[l.name] for l in tree.get_leaves()))
        d.append(draw.Line(rx, roy + row_h - 54, tipx, roy + row_h - 54, stroke=MUTED,
                           stroke_width=1.0, stroke_dasharray="2,3"))

    # shared substitution axis under the three trees
    base = rows_top + len(FAMILIES) * row_h - 40
    xa = lambda v: rx + (v / shared_maxd) * rpw           # noqa: E731
    d.append(draw.Line(xa(0), base, rx + rpw, base, stroke=INK, stroke_width=1.5))
    k = 0
    while k * 0.5 <= shared_maxd:
        v = k * 0.5
        d.append(draw.Line(xa(v), base, xa(v), base + 5, stroke=INK, stroke_width=1.5))
        d.append(draw.Text(f"{v:g}", FS_TICK, xa(v), base + 20, font_family=FONT,
                           text_anchor="middle", fill=MUTED))
        k += 1
    d.append(draw.Text("substitutions / site  (shared scale)", FS_TICK, rx + rpw / 2, base + 42,
                       font_family=FONT, text_anchor="middle", fill=INK, font_style="italic"))

    d.append(draw.Text("same colour pattern (R_b shared) -- different overall length (s_g)",
                       FS_ANNOT, rx + rpw / 2, base + 66, font_family=FONT, text_anchor="middle",
                       fill=INK, font_weight="bold"))

    name = "seq_gene_lineage"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png")


def _rate_bar(d, x, y, w, h):
    # discrete solid-colour cells, NOT a LinearGradient: macOS Preview flattens an axial-shading
    # gradient to a single colour, so a gradient colour bar reads as solid there. Cells render
    # identically everywhere.
    n = 64
    for i in range(n):
        t = i / (n - 1)
        d.append(draw.Rectangle(x + (i / n) * w, y, w / n + 0.7, h,
                                fill=C.rate_hex(C.RATE_LO * (C.RATE_HI / C.RATE_LO) ** t), stroke="none"))
    d.append(draw.Rectangle(x, y, w, h, fill="none", stroke=INK, stroke_width=0.8))
    for tx, lab, anc in ((x, "slow", "start"), (x + w / 2, "R_b = 1", "middle"), (x + w, "fast", "end")):
        d.append(draw.Text(lab, FS_TICK, tx, y + h + 23, font_family=FONT, text_anchor=anc, fill=MUTED))


if __name__ == "__main__":
    main()
