"""Figure: the couplings of the coevolve mode with SEQUENCES as a fourth node -- colour version.

Extends the S / T / G triangle (species, trait, gene content) with a fourth process, sequences (Q).
The three-node triangle is unchanged: six directed edges, three joint (double-headed) models, and the
two edges *into* S drawn heavy because they make the tree an OUTPUT.

Sequences (Σ) join as a TARGET-ONLY node: a trait drives selection (dN/dS) and substitution speed on
sequences (T->Σ), and gene content drives post-duplication relaxed selection (G->Σ). A sequence does
not drive another level, so Σ has arrows coming *in* and none going out, and no joint double-arrow.
Per the coevolve-grammar design (docs/design/coevolve-grammar.md), sequences ride GENE trees, so the
species-sequence edge (S-Σ) is FORBIDDEN -- Σ connects only to its tier neighbours T and G.

Run:  python figures/scripts/fig_coevolve_modes4.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw

from zombi_style import FONT, INK, FS_TITLE, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1240, 1055
R = 60                                        # node radius

# palette for the four processes -- teal / sage / terracotta, plus a muted violet for sequences.
NODE_COL = {
    "S": "#2f7d84",   # species    -- ZOMBI teal
    "T": "#5f7d34",   # traits     -- sage / olive green
    "G": "#b5654a",   # genes      -- terracotta
    "Σ": "#5a4e86",   # sequences  -- muted violet (the new target-only node)
}

# Each S/T/G node-PAIR gets one hue as a (dark, light) pair; the arrow is drawn in its shade, the LABEL
# always the dark shade so it reads on white.
PAIR = {
    "ST": ("#446b28", "#9ec27c"),   # species <-> traits : green      -- echoes T
    "SG": ("#1f5c63", "#86bdc1"),   # species <-> genes  : teal       -- echoes S
    "TG": ("#93492e", "#dc9f80"),   # traits  <-> genes  : terracotta -- echoes G
}
SEQ_DARK, SEQ_LIGHT = "#4a3f70", "#a79bcf"   # sequence-target edges, in Q's violet
DBL = "#6b6b6b"                       # bidirectional (joint) double-headed arrows -- neutral grey
LW_INTO_S, LW_OVERLAY = 5.2, 3.2      # weight encodes into-S (tree is an output) vs overlay
FS_MODEL = FS_TICK - 1
FS_NODE = FS_TITLE + 6
FS_NAME = FS_TICK - 2


def edge_point(cx, cy, tx, ty, r):
    a = math.atan2(ty - cy, tx - cx)
    return cx + r * math.cos(a), cy + r * math.sin(a)


def _head(d, hx, hy, ax, ay, col, ah=28.0):
    ang = math.atan2(hy - ay, hx - ax)
    d.append(draw.Lines(hx, hy, hx - ah * math.cos(ang - 0.42), hy - ah * math.sin(ang - 0.42),
                        hx - ah * math.cos(ang + 0.42), hy - ah * math.sin(ang + 0.42),
                        close=True, fill=col))


def arrow(d, c1, c2, side, bend, label, col, label_col, lw, label_gap=34, label_dx=0.0, label_dy=0.0):
    x1, y1 = c1
    x2, y2 = c2
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    dx, dy = x2 - x1, y2 - y1
    L = math.hypot(dx, dy)
    px, py = -dy / L, dx / L
    if py > 0 or (abs(py) < 1e-9 and px > 0):
        px, py = -px, -py
    ox, oy = px * side, py * side
    cx, cy = mx + ox * bend, my + oy * bend
    sx, sy = edge_point(x1, y1, cx, cy, R)
    ex, ey = edge_point(x2, y2, cx, cy, R)
    p = draw.Path(fill="none", stroke=col, stroke_width=lw, stroke_linecap="round")
    p.M(sx, sy).Q(cx, cy, ex, ey)
    d.append(p)
    _head(d, ex, ey, cx, cy, col)
    lx, ly = cx + ox * label_gap + label_dx, cy + oy * label_gap + label_dy
    d.append(draw.Text(label, FS_MODEL, lx, ly, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill=label_col, font_weight="bold"))


def biarrow(d, c1, c2, label, col=DBL, lw=3.4):
    x1, y1 = c1
    x2, y2 = c2
    sx, sy = edge_point(x1, y1, x2, y2, R)
    ex, ey = edge_point(x2, y2, x1, y1, R)
    d.append(draw.Line(sx, sy, ex, ey, stroke=col, stroke_width=lw, stroke_linecap="round"))
    _head(d, sx, sy, ex, ey, col, ah=24.0)
    _head(d, ex, ey, sx, sy, col, ah=24.0)
    mx, my = (sx + ex) / 2, (sy + ey) / 2
    tw = 0.55 * FS_MODEL * len(label) + 16
    th = FS_MODEL + 12
    d.append(draw.Rectangle(mx - tw / 2, my - th / 2, tw, th, fill="white", stroke="none"))
    d.append(draw.Text(label, FS_MODEL, mx, my, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill=col, font_weight="bold"))


def node(d, cx, cy, letter, name, dashed=False):
    col = NODE_COL[letter]
    d.append(draw.Circle(cx, cy, R, fill=col, stroke="white", stroke_width=3.2))
    ring = dict(stroke=INK, stroke_width=1.6)
    if dashed:
        ring["stroke_dasharray"] = "5 5"
    d.append(draw.Circle(cx, cy, R + 6, fill="none", **ring))     # outer ring; dashed = target only
    d.append(draw.Text(letter, FS_NODE, cx, cy - 7, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", font_weight="bold", fill="white"))
    d.append(draw.Text(name, FS_NAME, cx, cy + 27, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill="white"))


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("The coevolution diamond: species, traits, genomes, sequences", FS_TITLE, W / 2, 46,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))
    d.append(draw.Text("sequences ride gene trees, so there is no species-sequence edge", FS_TICK,
                       W / 2, 80, font_family=FONT, text_anchor="middle", font_style="italic",
                       fill="#777"))

    # a diamond: S top, T left, G right, Q bottom (downstream, target only)
    S = (W / 2, 292)
    T = (W / 2 - 322, 602)
    G = (W / 2 + 322, 602)
    Q = (W / 2, 934)

    sd, sl = PAIR["ST"]
    gd, gl = PAIR["SG"]
    td, tl = PAIR["TG"]

    # --- the S/T/G triangle: six directed edges (dark = one way, light = the other) ---
    arrow(d, T, S, +1, 60, "trait-driven",     sd, sd, LW_INTO_S)                            # T->S into S
    arrow(d, S, T, -1, 60, "cladogenetic",     sl, sd, LW_OVERLAY, label_dy=-24)             # S->T overlay
    arrow(d, G, S, -1, 60, "key innovation",   gd, gd, LW_INTO_S, label_dx=-26, label_dy=10) # G->S into S
    arrow(d, S, G, +1, 60, "punctuational",    gl, gd, LW_OVERLAY)                           # S->G overlay
    arrow(d, T, G, +1, 46, "trait-linked",     td, td, LW_OVERLAY, label_gap=6, label_dy=-6) # T->G
    arrow(d, G, T, -1, 46, "gene-conditioned", tl, td, LW_OVERLAY, label_gap=6, label_dy=-6) # G->T

    biarrow(d, T, S, "ClaSSE")
    biarrow(d, G, S, "co-diversification")
    biarrow(d, T, G, "trait-gene feedback")

    # --- sequences (Σ): TARGET ONLY. Arrows come IN from T and G; none go out; no joint. ---
    arrow(d, T, Q, -1, 40, "selection (dN/dS)", SEQ_LIGHT, SEQ_DARK, LW_OVERLAY, label_gap=26, label_dx=-8)
    arrow(d, G, Q, +1, 40, "relaxed selection", SEQ_LIGHT, SEQ_DARK, LW_OVERLAY, label_gap=26, label_dx=8)

    node(d, *S, "S", "species")
    node(d, *T, "T", "traits")
    node(d, *G, "G", "genomes")
    node(d, *Q, "Σ", "sequences", dashed=True)

    # legend, top-left
    lx, ly = 60, 108
    d.append(draw.Line(lx, ly, lx + 46, ly, stroke=INK, stroke_width=LW_INTO_S, stroke_linecap="round"))
    d.append(draw.Text("into S: the tree is grown as an OUTPUT", FS_TICK, lx + 60, ly,
                       font_family=FONT, text_anchor="start", dominant_baseline="central", fill=INK))
    d.append(draw.Line(lx, ly + 32, lx + 46, ly + 32, stroke=INK, stroke_width=LW_OVERLAY,
                       stroke_linecap="round"))
    d.append(draw.Text("overlay: runs on a given tree", FS_TICK, lx + 60, ly + 32, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))
    d.append(draw.Line(lx + 7, ly + 64, lx + 39, ly + 64, stroke=DBL, stroke_width=3.4,
                       stroke_linecap="round"))
    _head(d, lx + 7, ly + 64, lx + 39, ly + 64, DBL, ah=9.0)
    _head(d, lx + 39, ly + 64, lx + 7, ly + 64, DBL, ah=9.0)
    d.append(draw.Text("both edges at once: the joint (bidirectional) model", FS_TICK, lx + 60,
                       ly + 64, font_family=FONT, text_anchor="start", dominant_baseline="central",
                       fill=DBL))
    d.append(draw.Circle(lx + 23, ly + 96, 11, fill="none", stroke=INK, stroke_width=1.6,
                         stroke_dasharray="4 4"))
    d.append(draw.Text("sequences (Σ): a coupling target only (never drives)", FS_TICK,
                       lx + 60, ly + 96, font_family=FONT, text_anchor="start",
                       dominant_baseline="central", fill=NODE_COL["Σ"]))

    name = "coevolve_modes4"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png")


if __name__ == "__main__":
    render()
