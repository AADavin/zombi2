"""Figure: the nucleotide (sequence-level) model — root-anchored segments.

A nucleotide genome is a *circular* string of bases (drawn linearly, cut at the
origin); events cut it at breakpoints and rearrange arcs. Every present-day base
still knows the ancestral coordinate it came from (``to_cells()``), so a
present-day genome is a MOSAIC of ancestral intervals, some inverted.

Top bar = ancestral genome with a position gradient. Bottom bar = one leaf genome
painted BY ancestral position, so collinear stretches keep the gradient and an
inversion shows a reversed gradient (with a backward strand arrow). Ribbons join
each segment to the ancestral interval it came from (twisted for inversions).

Rendered in colour and in black & white.  Run:  python figures/scripts/fig_nucleotide.py
"""

from __future__ import annotations

from pathlib import Path

import cairosvg
import drawsvg as draw

from zombi2 import BirthDeath, simulate_nucleotide_genomes, simulate_species_tree

from zombi_style import FONT, INK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1040, 470
XL, XR = 96, 944
Y_ANC, Y_PRES, BAR_H = 108, 336, 32
ROOT_LEN, TREE_SEED, SIM_SEED = 100, 2, 4

RAINBOW = [(0.0, (60, 90, 190)), (0.2, (45, 170, 200)), (0.42, (75, 185, 95)),
           (0.6, (228, 208, 60)), (0.78, (238, 146, 45)), (1.0, (210, 55, 55))]
GREY = [(0.0, (38, 38, 38)), (1.0, (230, 230, 230))]


def ramp(stops, t):
    t = max(0.0, min(1.0, t))
    for (t0, c0), (t1, c1) in zip(stops, stops[1:]):
        if t <= t1:
            f = (t - t0) / (t1 - t0) if t1 > t0 else 0.0
            return tuple(int(round(c0[i] + (c1[i] - c0[i]) * f)) for i in range(3))
    return stops[-1][1]


def hexc(rgb):
    return "#%02x%02x%02x" % tuple(rgb)


def pos_color(t, mono=False):
    return hexc(ramp(GREY if mono else RAINBOW, t))


def strand_arrow(d, cx, cy, strand, aw):
    """A haloed arrow (white outline + black core) — visible on any background."""
    x0, x1 = cx - aw * strand, cx + aw * strand
    d.append(draw.Line(x0, cy, x1, cy, stroke="white", stroke_width=5, stroke_linecap="round"))
    d.append(draw.Line(x0, cy, x1 - 6 * strand, cy, stroke="black", stroke_width=2))
    d.append(draw.Lines(x1, cy, x1 - 9 * strand, cy - 5, x1 - 9 * strand, cy + 5,
                        close=True, fill="black", stroke="white", stroke_width=0.8))


def pick_leaf(r):
    best = None
    for node, g in r.leaf_genomes.items():
        inv = sum(1 for s in g._segments if s.strand < 0)
        score = (2 <= inv <= 3, 3 <= len(g._segments) <= 6, inv, -len(g._segments))
        if best is None or score > best[0]:
            best = (score, g)
    return best[1]


def render(mono):
    tree = simulate_species_tree(BirthDeath(1.0, 0.0), n_tips=6, age=1.0, seed=TREE_SEED)
    r = simulate_nucleotide_genomes(tree, inversion=0.025, root_length=ROOT_LEN, extension=0.93, seed=SIM_SEED)
    g = pick_leaf(r)
    L = r.root_length
    segs = g._segments

    def x_of(p):
        return XL + (p / L) * (XR - XL)

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    dx = (XR - XL) / L

    for p in range(L):                                     # ancestral bar: position gradient
        d.append(draw.Rectangle(x_of(p), Y_ANC, dx + 0.6, BAR_H, fill=pos_color((p + 0.5) / L, mono)))
    d.append(draw.Rectangle(XL, Y_ANC, XR - XL, BAR_H, fill="none", stroke=INK, stroke_width=1.4))

    cells = g.to_cells()                                   # present bar: painted by ancestral pos
    for q in range(L):
        d.append(draw.Rectangle(x_of(q), Y_PRES, dx + 0.6, BAR_H, fill=pos_color((cells[q][1] + 0.5) / L, mono)))
    d.append(draw.Rectangle(XL, Y_PRES, XR - XL, BAR_H, fill="none", stroke=INK, stroke_width=1.4))

    yab, ypt = Y_ANC + BAR_H, Y_PRES
    q = 0
    for s in segs:
        a, b, q0, q1 = s.src_start, s.src_end, q, q + s.length
        q = q1
        xa0, xa1 = x_of(a), x_of(b)
        xp0, xp1 = (x_of(q0), x_of(q1)) if s.strand > 0 else (x_of(q1), x_of(q0))
        dy = ypt - yab
        col = pos_color(((a + b) / 2) / L, mono)
        p = draw.Path(fill=col, fill_opacity=0.30, stroke="none")
        p.M(xa0, yab).C(xa0, yab + dy * 0.5, xp0, ypt - dy * 0.5, xp0, ypt)
        p.L(xp1, ypt).C(xp1, ypt - dy * 0.5, xa1, yab + dy * 0.5, xa1, yab).Z()
        d.append(p)
        d.append(draw.Line(x_of(q0), Y_PRES - 4, x_of(q0), Y_PRES + BAR_H + 4, stroke=INK, stroke_width=1.2))
        cxp = (x_of(q0) + x_of(q1)) / 2
        strand_arrow(d, cxp, Y_PRES + BAR_H / 2, s.strand, min(26, (x_of(q1) - x_of(q0)) * 0.5))

    for p in [0, L // 4, L // 2, 3 * L // 4, L]:
        d.append(draw.Line(x_of(p), Y_ANC - 4, x_of(p), Y_ANC, stroke=INK, stroke_width=1.2))
        d.append(draw.Text(str(p), 12, x_of(p), Y_ANC - 9, font_family=FONT, text_anchor="middle", fill=INK))
    d.append(draw.Text("ancestral (root) genome — base coordinate", 15, XL, Y_ANC - 26,
                       font_family=FONT, text_anchor="start", font_weight="bold", fill=INK))
    d.append(draw.Text("present-day genome — a mosaic of ancestral segments", 15, XL,
                       Y_PRES + BAR_H + 30, font_family=FONT, text_anchor="start", font_weight="bold", fill=INK))
    d.append(draw.Text("each base keeps its ancestral coordinate;  an inversion reverses a segment "
                       "(twisted ribbon) and flips its strand (arrow)", 13, XL, H - 22,
                       font_family=FONT, text_anchor="start", fill="#555"))

    stem = OUT_DIR / ("nucleotide_segments_bw" if mono else "nucleotide_segments")
    stem.mkdir(parents=True, exist_ok=True)
    stem = stem / stem.name
    stem.with_suffix(".svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(stem.with_suffix(".png")), scale=300 / 72.0)
    return len(segs), sum(1 for s in segs if s.strand < 0)


def main():
    for mono in (False, True):
        n, inv = render(mono)
    print(f"wrote nucleotide_segments (+_bw)  (L={ROOT_LEN}, {n} segments, {inv} inverted)")


if __name__ == "__main__":
    main()
