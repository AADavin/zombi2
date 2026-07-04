"""Figure: a complete species tree with extinct lineages (forward ZOMBI2 sim).

The forward birth-death process keeps the lineages that died out. We draw the
*complete* tree so the surviving skeleton and the extinct history read as one
picture — the point that motivates ZOMBI2's forward engine (backward simulation
only ever sees the solid tree; the dashes are what really happened).

Encoding
--------
Any segment pointing at a node ``x`` is drawn **dashed** iff ``x`` has no surviving
descendant. Applied to horizontal branches and to each half of a fork's vertical
connector, this gives exactly:

  * both children survive        -> fully solid fork
  * one survives, one dies out   -> half-solid / half-dashed fork
  * speciation inside a dead clade -> fully dashed fork

Extinct tips get a Greek-letter name and their dashed branch simply stops at the
death time; survivors get Latin letters and reach the present.

Extinct tips are found geometrically (a leaf whose root-distance falls short of
the present), so this works on any complete-tree Newick.

Run:  python figures/scripts/fig_species_tree_extinct.py
"""

from __future__ import annotations

import string
from pathlib import Path

import drawsvg as draw
import phylustrator as ph
from phylustrator.io import read_newick

from zombi_style import INK, species_style

FIG_DIR = Path(__file__).resolve().parent.parent
TREE_NWK = FIG_DIR / "species_tree_extinct" / "species_tree.nwk"
OUT_STEM = FIG_DIR / "species_tree_extinct" / "species_tree"

GREEK = list("αβγδεζηθικλμνξοπρστυφχψω")
DASH = "6,5"                 # dash pattern for extinct lineages
PRESENT = "#c9c9c9"          # faint "today" reference line


def annotate_depths(tree) -> float:
    """Set ``node.depth`` = root-to-node distance; return the max (the present)."""
    max_depth = 0.0
    for n in tree.traverse("preorder"):
        n.depth = 0.0 if n.is_root() else n.up.depth + n.dist
        max_depth = max(max_depth, n.depth)
    return max_depth


def mark_survival(tree, present: float, tol: float) -> None:
    """Flag ``is_extant`` (leaf reaches the present) and ``has_survivor``
    (subtree contains an extant leaf)."""
    for n in tree.traverse("postorder"):
        if n.is_leaf():
            n.is_extant = abs(n.depth - present) <= tol
            n.has_survivor = n.is_extant
        else:
            n.has_survivor = any(c.has_survivor for c in n.children)


def _seg(d, x1, y1, x2, y2, dashed: bool) -> None:
    """Draw one branch segment; dashed (butt caps) if extinct, else solid (round)."""
    sw = d.style.branch_stroke_width
    if dashed:
        d.drawing.append(draw.Line(x1, y1, x2, y2, stroke=INK, stroke_width=sw,
                                   stroke_dasharray=DASH, stroke_linecap="butt"))
    else:
        d.drawing.append(draw.Line(x1, y1, x2, y2, stroke=INK, stroke_width=sw,
                                   stroke_linecap="round"))


def draw_skeleton(d, tree) -> None:
    """Draw the tree by hand so each segment can be solid or dashed independently.

    A segment is dashed iff the node it points at has no surviving descendant.
    """
    stub = d.style.root_stub_length
    for n in tree.traverse("postorder"):
        x, y = n.coordinates
        if n.is_root():
            _seg(d, x - stub, y, x, y, dashed=not n.has_survivor)
        else:
            px, _ = n.up.coordinates
            _seg(d, px, y, x, y, dashed=not n.has_survivor)          # branch above n
        if not n.is_leaf():
            for c in n.children:                                      # fork: one stub per child
                _seg(d, x, y, x, c.coordinates[1], dashed=not c.has_survivor)


def add_legend(d, x, y, font_size=15) -> None:
    """Line-swatch legend: solid = extant lineage, dashed = extinct lineage."""
    sw = d.style.branch_stroke_width
    L = 30
    d.drawing.append(draw.Text("Lineages", font_size + 2, x, y, font_weight="bold",
                               font_family=d.style.font_family, text_anchor="start"))
    rows = [("Extant lineage", False), ("Extinct lineage", True)]
    cy = y + font_size * 1.7
    for label, dashed in rows:
        kw = dict(stroke=INK, stroke_width=sw)
        kw.update(stroke_dasharray=DASH, stroke_linecap="butt") if dashed else \
            kw.update(stroke_linecap="round")
        d.drawing.append(draw.Line(x, cy, x + L, cy, **kw))
        d.drawing.append(draw.Text(label, font_size, x + L + 10, cy,
                                   font_family=d.style.font_family, text_anchor="start",
                                   dominant_baseline="middle"))
        cy += font_size * 1.7


def main() -> None:
    tree = read_newick(TREE_NWK)
    present = annotate_depths(tree)
    mark_survival(tree, present, tol=1e-6 * present)

    leaves = tree.get_leaves()                     # top-to-bottom order
    extant = [l for l in leaves if l.is_extant]
    extinct = [l for l in leaves if not l.is_extant]

    # Names: survivors get Latin letters, extinct tips get Greek letters.
    for leaf, letter in zip(extant, string.ascii_uppercase):
        leaf.name = letter
    for leaf, g in zip(extinct, GREEK):
        leaf.name = g

    style = species_style(height=800)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()

    # faint "present" reference (behind everything)
    ys = [l.y_coord for l in leaves]
    present_x = d.root_x + present * d.sf
    d.drawing.append(draw.Line(present_x, min(ys) - 14, present_x, max(ys) + 14,
                               stroke=PRESENT, stroke_width=1.0, stroke_dasharray="2,4"))

    draw_skeleton(d, tree)
    d.add_leaf_names(color=INK, padding=12)

    ticks = [round(present * i / 4, 6) for i in range(5)]
    d.add_time_axis(ticks=ticks, tick_labels=[f"{t:.2f}" for t in ticks],
                    label="Time (root to present)", tick_size=6.0, padding=14.0,
                    stroke_width=1.6)

    # lower-left is the open quadrant here; keeps clear of the basal extinct tip α
    add_legend(d, x=-style.width / 2 + 34, y=style.height / 2 - 150)

    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    print(f"wrote {OUT_STEM}.svg / .png  ({len(extant)} extant + {len(extinct)} extinct tips)")


if __name__ == "__main__":
    main()
