"""Figure: nucleotide events on a 3-leaf tree -> segments -> a tree per segment.

A 12-base root genome descends a species tree ((A,B),C). One event falls on each leaf
branch: a duplication on A, an inversion on B, a loss on C. Read straight from a real
``simulate_nucleotide_genomes`` run (seed 2706); each event happens to act on exactly
one block, so the story is clean.

Three panels, top to bottom:

  A  SEGMENTS ARE CREATED. The initial (root) genome is one smooth gradient with no
     divisions. Each leaf keeps only the segment its own event carved out.
  B  SEGMENTS PROPAGATE. The same tree again, but now the genome is cut into its blocks
     (each a gradient slice of the initial genome) and every leaf is a mosaic of them —
     duplicated on A, reversed on B, one missing on C.
  C  A TREE PER SEGMENT. Each block carries its own genealogy: duplication adds a tip,
     loss prunes one, inversion leaves the tree (only the strand flips), the rest are the
     species tree.

Rendered gradient (each segment a gradient slice), colour and B&W (solid per-block hues).
Run:  python figures/scripts/fig_nucleotide_tree.py
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

import cairosvg
import drawsvg as draw

from zombi2 import simulate_nucleotide_genomes
from zombi2.tree import Tree, TreeNode

from fig_nucleotide import GREY, pos_color, ramp  # noqa: F401
from zombi_style import FONT, INK

OUT_DIR = Path(__file__).resolve().parent.parent

# Panels A and B stack on the LEFT (x < 780); panel C is a vertical column on the RIGHT.
W, H = 1180, 680
TITLE_X = 40
SPLIT_X = 788                                  # divider between the A/B column and panel C

BLOCK_COL = {0: "#6b8fb0", 1: "#7ba884", 2: "#c8a75c", 3: "#b47d6e", 4: "#8f7bb0"}
BLOCK_BW = {0: "#3a3a3a", 1: "#5f5f5f", 2: "#838383", 3: "#a6a6a6", 4: "#c8c8c8"}

BLOCK: dict = {}
BLOCK_LEN: dict = {}
ROOT_L = 12


# --------------------------------------------------------------------------- data
def build_tree() -> Tree:
    root = TreeNode("n1", 0.0)
    n2 = TreeNode("n2", 0.4); root.add_child(n2)
    C = TreeNode("C", 1.0); root.add_child(C)
    A = TreeNode("A", 1.0); n2.add_child(A)
    B = TreeNode("B", 1.0); n2.add_child(B)
    return Tree(root, 1.0)


def classify(ext, inverted):
    cnt = Counter(re.findall(r"([ABC])_g\d+", ext or ""))
    if any(v > 1 for v in cnt.values()):
        return "dup"
    if set(cnt) != {"A", "B", "C"}:
        return "loss"
    return "inv" if inverted else "species"


def simulate():
    r = simulate_nucleotide_genomes(build_tree(), inversion=0.02, duplication=0.02,
                                    loss=0.02, root_length=12, extension=0.6, seed=2706)
    blocks = {a.block_id: a for a in r.blocks}
    gt, hist = r.block_gene_trees(), r.block_histories()
    kind = {aid: classify(gt[aid][1], bool(hist[aid])) for aid in blocks}

    L = max(a.end for a in blocks.values())
    prov = r.registry.provenance
    block_leaf = {aid: None for aid in blocks}
    leaf_event, leaf_breaks = {}, {}
    for rec in r.event_log:
        if rec.event.name not in ("DUPLICATION", "LOSS", "INVERSION"):
            continue
        segs = [prov[op.gid] for op in rec.genes if op.gid in prov]
        ss, se = min(s[1] for s in segs), max(s[2] for s in segs)
        aff = next(aid for aid, a in blocks.items() if ss <= a.start and a.end <= se)
        block_leaf[aff] = rec.branch
        leaf_event[rec.branch] = (kind[aff], aff)
        leaf_breaks[rec.branch] = sorted(b for b in (ss, se) if 0 < b < L)
    mosaics = {leaf.name: r.leaf_mosaic(leaf) for leaf in r.leaf_genomes}
    return blocks, kind, block_leaf, mosaics, leaf_event, leaf_breaks, L


# --------------------------------------------------------------------------- colour
def col(aid, mode):
    if mode == "color":
        return BLOCK_COL[aid]
    if mode == "bw":
        return BLOCK_BW[aid]
    a = BLOCK[aid]
    return pos_color(((a.start + a.end) / 2) / ROOT_L, True)


def lum(h):
    return 0.299 * int(h[1:3], 16) + 0.587 * int(h[3:5], 16) + 0.114 * int(h[5:7], 16)


# --------------------------------------------------------------------------- glyphs
def make_hatch(d, spacing=4.5, sw=1.1):
    p = draw.Pattern(spacing, spacing, id="ev_hatch", patternUnits="userSpaceOnUse")
    p.append(draw.Rectangle(0, 0, spacing, spacing, fill="white"))
    for off in (-spacing, 0, spacing):
        p.append(draw.Line(off, spacing, off + spacing, 0, stroke=INK, stroke_width=sw))
    d.append(p)
    return p


def gl_dup(d, cx, cy, hatch, r=6):
    d.append(draw.Rectangle(cx - r, cy - r, 2 * r, 2 * r, fill=hatch, stroke=INK, stroke_width=1.5))


def gl_loss(d, cx, cy, hatch, r=6):
    d.append(draw.Circle(cx, cy, r, fill=hatch, stroke=INK, stroke_width=1.5))


def gl_inv(d, cx, cy, r=7):
    p = draw.Path(fill="none", stroke=INK, stroke_width=2.2)
    p.M(cx + 6, cy + 3).C(cx + 8, cy - 6, cx - 8, cy - 6, cx - 6, cy + 3)
    d.append(p)
    d.append(draw.Lines(cx - 6, cy + 3, cx - 9, cy - 2, cx - 2, cy - 2, close=True, fill=INK))


def draw_glyph(d, kind, cx, cy, hatch, r=6):
    if kind == "dup":
        gl_dup(d, cx, cy, hatch, r)
    elif kind == "loss":
        gl_loss(d, cx, cy, hatch, r)
    else:
        gl_inv(d, cx, cy, r + 1)


def flip_mark(d, cx, cy, w=11):
    d.append(draw.Line(cx + w, cy, cx - w, cy, stroke="white", stroke_width=4.5, stroke_linecap="round"))
    d.append(draw.Line(cx + w, cy, cx - w + 5, cy, stroke=INK, stroke_width=1.8))
    d.append(draw.Lines(cx - w, cy, cx - w + 7, cy - 4, cx - w + 7, cy + 4, close=True,
                        fill=INK, stroke="white", stroke_width=0.7))


# --------------------------------------------------------------------------- segment / genome strips
def draw_slice(d, x, y, w, h, aid, mode, rev=False, faded=False):
    """One block's content: a gradient slice of the initial genome (or a solid hue)."""
    a = BLOCK[aid]
    if mode == "gradient":
        sub = w / a.length
        for k in range(a.length):
            src = a.end - 1 - k if rev else a.start + k
            g = ramp(GREY, (src + 0.5) / ROOT_L)
            if faded:
                g = tuple(int(v + (255 - v) * 0.72) for v in g)
            d.append(draw.Rectangle(x + k * sub, y, sub + 0.7, h, fill="#%02x%02x%02x" % g))
    else:
        c = col(aid, mode)
        if faded:
            r_, g_, b_ = (int(c[i:i + 2], 16) for i in (1, 3, 5))
            c = "#%02x%02x%02x" % tuple(int(v + (255 - v) * 0.72) for v in (r_, g_, b_))
        d.append(draw.Rectangle(x, y, w, h, fill=c))
    d.append(draw.Rectangle(x, y, w, h, fill="none", stroke="white", stroke_width=1.0))


def _smooth_gradient(d, x0, y, w, h):
    sub = w / ROOT_L
    for k in range(ROOT_L):
        g = ramp(GREY, (k + 0.5) / ROOT_L)
        d.append(draw.Rectangle(x0 + k * sub, y, sub + 0.7, h, fill="#%02x%02x%02x" % g))


def div_line(d, x, y, h):
    """A division line that reads on any shade: white core with thin dark edges."""
    d.append(draw.Line(x, y, x, y + h, stroke="#4a4a4a", stroke_width=2.6))
    d.append(draw.Line(x, y, x, y + h, stroke="white", stroke_width=1.3))


def _block_lines(d, x0, y, w, h, breaks):
    for bp in breaks:
        div_line(d, x0 + (bp / ROOT_L) * w, y, h)


def ancestral_genome(d, x0, y, w, mode, divided, h=18):
    """The initial genome: smooth (divided=False) or cut into its block squares."""
    if not divided:
        if mode == "gradient":
            _smooth_gradient(d, x0, y, w, h)
        else:
            d.append(draw.Rectangle(x0, y, w, h, fill="#d8cebd" if mode == "color" else "#dcdcdc"))
    else:
        for aid in sorted(BLOCK):
            a = BLOCK[aid]
            draw_slice(d, x0 + (a.start / ROOT_L) * w, y, (a.length / ROOT_L) * w, h, aid, mode)
        _block_lines(d, x0, y, w, h, [a.start for a in BLOCK.values() if a.start > 0])
    d.append(draw.Rectangle(x0, y, w, h, fill="none", stroke=INK, stroke_width=1.1))


def native_genome(d, cx, y, mode, breaks, seg_aid, h=18):
    """Panel A: the full genome (correct gradient colours) with only this leaf's own cuts."""
    w = ROOT_L * 9.0
    x0 = cx - w / 2
    if mode == "gradient":
        _smooth_gradient(d, x0, y, w, h)
    else:
        d.append(draw.Rectangle(x0, y, w, h, fill="#d8cebd" if mode == "color" else "#dcdcdc"))
        a = BLOCK[seg_aid]                                        # colour just the carved segment
        draw_slice(d, x0 + (a.start / ROOT_L) * w, y, (a.length / ROOT_L) * w, h, seg_aid, mode)
    _block_lines(d, x0, y, w, h, breaks)
    a = BLOCK[seg_aid]                                            # bracket marks "their" segment
    ux0, ux1, yb = x0 + (a.start / ROOT_L) * w, x0 + (a.end / ROOT_L) * w, y + h + 5
    d.append(draw.Line(ux0, yb, ux1, yb, stroke=INK, stroke_width=2))
    d.append(draw.Line(ux0, yb, ux0, yb - 4, stroke=INK, stroke_width=2))
    d.append(draw.Line(ux1, yb, ux1, yb - 4, stroke=INK, stroke_width=2))
    d.append(draw.Rectangle(x0, y, w, h, fill="none", stroke=INK, stroke_width=1))


def mosaic_bar(d, cx, y, cells, mode, h=22, unit=9.0):
    """A leaf genome: its blocks in order, each a gradient slice (reversed if inverted).

    White lines are drawn at every block limit (the propagated breakpoints)."""
    total = sum(BLOCK_LEN[aid] for aid, _ in cells) * unit
    x = cx - total / 2
    edges = [x]
    for aid, strand in cells:
        w = BLOCK_LEN[aid] * unit
        draw_slice(d, x, y, w, h, aid, mode, rev=strand < 0)
        if strand < 0:
            flip_mark(d, x + w / 2, y + h + 8, w=min(10, w / 2 - 1))
        x += w
        edges.append(x)
    for bx in edges[1:-1]:                                       # every block limit -> division line
        div_line(d, bx, y, h)
    d.append(draw.Rectangle(cx - total / 2, y, total, h, fill="none", stroke=INK, stroke_width=1.0))


# --------------------------------------------------------------------------- species tree
AX, BX, CX = 230, 430, 700
N2X = (AX + BX) / 2
RX = (N2X + CX) / 2


def species_tree(d, y0, th):
    yl, yn2 = y0 + th, y0 + 0.4 * th
    for a, b, c, dd in ((N2X, y0, CX, y0), (N2X, y0, N2X, yn2), (CX, y0, CX, yl),
                        (AX, yn2, BX, yn2), (AX, yn2, AX, yl), (BX, yn2, BX, yl)):
        d.append(draw.Line(a, b, c, dd, stroke=INK, stroke_width=2.6))
    return yl, yn2


# --------------------------------------------------------------------------- panel A
def panel_a(d, blocks, leaf_event, leaf_breaks, mode, hatch):
    y0 = 130
    gw = ROOT_L * 9.0
    d.append(draw.Line(RX, 84, RX, y0, stroke=INK, stroke_width=1.1, stroke_dasharray="3,3"))
    ancestral_genome(d, RX - gw / 2, 64, gw, mode, divided=False)
    d.append(draw.Text("initial genome — one piece, a smooth gradient (no segments yet)", 15, RX, 56,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))

    yl, yn2 = species_tree(d, y0, 112)
    for name, lx in (("A", AX), ("B", BX), ("C", CX)):
        kind, aid = leaf_event[name]
        sy = yn2 + 0.5 * (yl - yn2) if name in ("A", "B") else y0 + 0.5 * (yl - y0)
        draw_glyph(d, kind, lx, sy, hatch)
        native_genome(d, lx, yl + 13, mode, leaf_breaks[name], aid)
        d.append(draw.Text(name, 20, lx, yl + 50, font_family=FONT, text_anchor="middle",
                           font_weight="bold", fill=INK))
    d.append(draw.Text("A — each event carves out one segment (only its own cut shown)", 18, TITLE_X, 34,
                       font_family=FONT, text_anchor="start", font_weight="bold", fill=INK))


# --------------------------------------------------------------------------- panel B
def panel_b(d, blocks, mosaics, mode):
    d.append(draw.Line(TITLE_X, 306, 760, 306, stroke="#dcdcdc", stroke_width=1.2))
    d.append(draw.Text("B — the segments propagate: every genome is a mosaic of the same blocks", 18,
                       TITLE_X, 326, font_family=FONT, text_anchor="start", font_weight="bold", fill=INK))
    y0 = 434
    gw = ROOT_L * 9.0
    d.append(draw.Line(RX, 378, RX, y0, stroke=INK, stroke_width=1.1, stroke_dasharray="3,3"))
    ancestral_genome(d, RX - gw / 2, 356, gw, mode, divided=True)
    d.append(draw.Text("initial genome — now cut into its blocks (gradient squares)", 15, RX, 346,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))

    yl, yn2 = species_tree(d, y0, 112)
    for name, lx in (("A", AX), ("B", BX), ("C", CX)):
        d.append(draw.Line(lx, yl, lx, yl + 14, stroke=INK, stroke_width=1.2))
        mosaic_bar(d, lx, yl + 14, mosaics[name], mode)
        d.append(draw.Text(name, 20, lx, yl + 54, font_family=FONT, text_anchor="middle",
                           font_weight="bold", fill=INK))


# --------------------------------------------------------------------------- panel C
# Panel C is a vertical column on the right: one mini gene-tree per block, stacked.
CX_C = 992                                      # tree centre in the right column
SLOT0, SLOT_H, MH = 88, 118, 64                 # first header y, row pitch, tree height
LEAF_X = {"A": -48, "B": -15, "C": 48}


def mini_tree(d, cx, top, aid, kind, leaf, mode, hatch):
    c = col(aid, mode)
    mleaf, mn2 = top + MH, top + 0.4 * MH
    xA, xB, xC = cx + LEAF_X["A"], cx + LEAF_X["B"], cx + LEAF_X["C"]
    n2 = (xA + xB) / 2
    d.append(draw.Line(n2, top, xC, top, stroke=c, stroke_width=3))
    d.append(draw.Line(n2, top, n2, mn2, stroke=c, stroke_width=3))
    d.append(draw.Line(xA, mn2, xB, mn2, stroke=c, stroke_width=3))
    starts = {"A": (xA, mn2), "B": (xB, mn2), "C": (xC, top)}
    tips = []
    for nm, (x, sy) in starts.items():
        if nm == leaf and kind == "loss":
            yl = sy + 0.34 * (mleaf - sy)
            d.append(draw.Line(x, sy, x, yl - 6, stroke=INK, stroke_width=2, stroke_dasharray="4,3"))
            gl_loss(d, x, yl, hatch, r=6)
        elif nm == leaf and kind == "dup":
            yd = sy + 0.30 * (mleaf - sy)
            x1, x2 = x - 9, x + 9
            d.append(draw.Line(x, sy, x, yd, stroke=c, stroke_width=3))
            d.append(draw.Line(x1, yd, x2, yd, stroke=c, stroke_width=3))
            d.append(draw.Line(x1, yd, x1, mleaf, stroke=c, stroke_width=3))
            d.append(draw.Line(x2, yd, x2, mleaf, stroke=c, stroke_width=3))
            gl_dup(d, x, yd, hatch, r=6)
            tips += [(x1, nm), (x2, nm)]
        else:
            d.append(draw.Line(x, sy, x, mleaf, stroke=c, stroke_width=3))
            tips.append((x, nm))
            # inversions do not change the genealogy -> no mark on the gene tree
    for tx, tn in tips:
        d.append(draw.Text(tn, 15, tx, mleaf + 16, font_family=FONT, text_anchor="middle",
                           font_weight="bold", fill=INK))


def panel_c(d, blocks, kind, block_leaf, mode, hatch):
    d.append(draw.Text("C — reconstruct a tree for every block", 17, SPLIT_X + 20, 34,
                       font_family=FONT, text_anchor="start", font_weight="bold", fill=INK))
    for i, aid in enumerate(sorted(blocks)):
        a = blocks[aid]
        hy = SLOT0 + i * SLOT_H                                     # header baseline for this block
        d.append(draw.Rectangle(CX_C - 96, hy - 9, 22, 13, fill="none", stroke=INK, stroke_width=0.8))
        draw_slice(d, CX_C - 96, hy - 9, 22, 13, aid, mode)         # the block's own gradient chip
        d.append(draw.Text(f"block {aid}", 14, CX_C - 66, hy, font_family=FONT, text_anchor="start",
                           dominant_baseline="central", font_weight="bold", fill=INK))
        d.append(draw.Text(f"[{a.start},{a.end})", 12, CX_C + 18, hy, font_family=FONT,
                           text_anchor="start", dominant_baseline="central", fill="#888"))
        mini_tree(d, CX_C, hy + 16, aid, kind[aid], block_leaf[aid], mode, hatch)


# --------------------------------------------------------------------------- render
def render(blocks, kind, block_leaf, mosaics, leaf_event, leaf_breaks, mode):
    global BLOCK, BLOCK_LEN
    BLOCK = blocks
    BLOCK_LEN = {aid: a.length for aid, a in blocks.items()}
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Line(SPLIT_X, 22, SPLIT_X, H - 16, stroke="#dcdcdc", stroke_width=1.2))
    hatch = make_hatch(d)
    panel_a(d, blocks, leaf_event, leaf_breaks, mode, hatch)
    panel_b(d, blocks, mosaics, mode)
    panel_c(d, blocks, kind, block_leaf, mode, hatch)

    name = {"gradient": "nucleotide_tree_gradient", "color": "nucleotide_tree",
            "bw": "nucleotide_tree_bw"}[mode]
    stem = OUT_DIR / name
    stem.mkdir(parents=True, exist_ok=True)
    stem = stem / name
    stem.with_suffix(".svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(stem.with_suffix(".png")), scale=300 / 72.0)


def main():
    blocks, kind, block_leaf, mosaics, leaf_event, leaf_breaks, L = simulate()
    global ROOT_L
    ROOT_L = L
    for mode in ("gradient", "color", "bw"):
        render(blocks, kind, block_leaf, mosaics, leaf_event, leaf_breaks, mode)
    print("wrote nucleotide_tree {gradient,color,bw}  leaf_event=", leaf_event,
          " leaf_breaks=", leaf_breaks)


if __name__ == "__main__":
    main()
