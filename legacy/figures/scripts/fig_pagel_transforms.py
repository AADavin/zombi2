"""Figure: Pagel's (1999) tree transforms -- lambda, kappa, delta.

One small tree, shown under each transform at three parameter values (a 3x3 grid of rectangular
cladograms). The transforms only rescale branch/node lengths, so the topology is fixed and the tip
order is shared across every cell -- what moves is where the internal nodes (and, for kappa, the
tips) sit in time:

  * lambda -- scales internal (shared) depths, tips pinned: 1 = original, 0 = star (independent tips)
  * delta  -- node depths ^delta: <1 pushes change early (toward the root), >1 late (toward the tips)
  * kappa  -- branch lengths ^kappa: 0 = speciational (equal branches, change per speciation event)

House style: near-black branches, ASCII-only labels (no Greek), title centred at the top.

Run:  python figures/scripts/fig_pagel_transforms.py
"""

from __future__ import annotations

from pathlib import Path

import cairosvg

import zombi2 as z
from zombi2.traits import pagel_delta, pagel_kappa, pagel_lambda
from zombi_style import INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT

FIG_DIR = Path(__file__).resolve().parent.parent
OUT_DIR = FIG_DIR / "pagel_transforms"
OUT_STEM = OUT_DIR / "pagel_transforms"

W, H = 1360, 1030
LEFT, TOP, RIGHT, BOT = 210, 120, 26, 54
BRANCH_W = 2.6

# the three transforms, each with three parameter values + a short per-cell tag
ROWS = [
    ("lambda", "phylogenetic signal", pagel_lambda,
     [(1.0, "1 (original)"), (0.5, "0.5"), (0.0, "0 (star tree)")]),
    ("delta", "early vs late change", pagel_delta,
     [(0.4, "0.4 (early)"), (1.0, "1 (original)"), (2.5, "2.5 (late)")]),
    ("kappa", "gradual vs speciational", pagel_kappa,
     [(1.0, "1 (original)"), (0.5, "0.5"), (0.0, "0 (speciational)")]),
]


def _t(x, y, s, size, *, anchor="start", weight="normal", color=INK, italic=False):
    st = ' font-style="italic"' if italic else ""
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" font-family="Helvetica" '
            f'text-anchor="{anchor}" dominant-baseline="middle" font-weight="{weight}"{st} '
            f'fill="{color}">{s}</text>')


def _line(x1, y1, x2, y2, *, w=BRANCH_W, color=INK):
    return (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{color}" '
            f'stroke-width="{w}" stroke-linecap="round" />')


def _leaves(tree):
    return [n for n in tree.nodes_preorder() if not n.children]


def draw_cladogram(o, tree, leaf_row, x0, y0, w, h):
    """Draw ``tree`` as a rectangular cladogram inside the box (x0, y0, w, h). ``leaf_row`` maps a
    leaf name to its (fixed) row index so every cell shares the same tip order. x = depth (node.time),
    scaled per-cell to fill the width; y = tip row."""
    n = len(leaf_row)
    padx, padtop, padbot = 22, 20, 20
    max_depth = max((nd.time for nd in tree.nodes_preorder()), default=1.0) or 1.0

    def xof(depth):
        return x0 + padx + (depth / max_depth) * (w - 2 * padx)

    def yof(row):
        return y0 + padtop + (row + 0.5) / n * (h - padtop - padbot)

    pos = {}

    def rec(node):
        if not node.children:
            y = yof(leaf_row[node.name])
        else:
            y = sum(rec(c) for c in node.children) / len(node.children)
        pos[id(node)] = (xof(node.time), y)
        return y

    rec(tree.root)

    rx, ry = pos[id(tree.root)]
    o.append(_line(x0 + 6, ry, rx, ry))                       # short root stub
    for node in tree.nodes_preorder():
        if not node.children:
            continue
        nx, ny = pos[id(node)]
        cys = [pos[id(c)][1] for c in node.children]
        o.append(_line(nx, min(cys), nx, max(cys)))            # vertical connector
        for c in node.children:
            cx, cy = pos[id(c)]
            o.append(_line(nx, cy, cx, cy))                    # horizontal to child
            o.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="3.0" fill="{INK}" />'
                     if not c.children else "")


def main():
    base = z.simulate_species_tree(z.BirthDeath(1.0, 0.2), n_tips=6, age=1.0,
                                   direction="backward", seed=3)
    leaf_row = {lf.name: i for i, lf in enumerate(_leaves(base))}

    cell_w = (W - LEFT - RIGHT) / 3
    cell_h = (H - TOP - BOT) / 3

    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
         f'viewBox="0 0 {W} {H}">',
         f'<rect x="0" y="0" width="{W}" height="{H}" fill="white" />']

    o.append(_t(W / 2, 52, "Pagel's tree transforms", FS_TITLE, anchor="middle", weight="bold"))

    for r, (name, blurb, fn, values) in enumerate(ROWS):
        ry0 = TOP + r * cell_h
        # row label (left margin): transform name + what it does
        o.append(_t(28, ry0 + cell_h / 2 - 12, name, FS_LABEL, weight="bold"))
        o.append(_t(28, ry0 + cell_h / 2 + 16, blurb, FS_ANNOT - 4, color=MUTED))
        for c, (val, tag) in enumerate(values):
            cx0 = LEFT + c * cell_w
            tree = fn(base, val)
            o.append(_t(cx0 + cell_w / 2, ry0 + 20, f"{name} = {tag}", FS_ANNOT - 2,
                        anchor="middle", color=INK))
            draw_cladogram(o, tree, leaf_row, cx0, ry0 + 38, cell_w, cell_h - 44)

    o.append(_t(W / 2, H - 22,
                "The transform reshapes the tree; the trait model is then run on the result.",
                FS_ANNOT - 4, anchor="middle", color=MUTED))
    o.append("</svg>")
    svg = "\n".join(x for x in o if x)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_STEM.with_suffix(".svg").write_text(svg, encoding="utf-8")
    cairosvg.svg2png(bytestring=svg.encode("utf-8"),
                     write_to=str(OUT_STEM.with_suffix(".png")), scale=300 / 72.0)
    print(f"wrote {OUT_STEM}.svg / .png")


if __name__ == "__main__":
    main()
