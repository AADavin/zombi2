"""Figure: a gene tree produced by ZOMBI2, with its reconciliation events marked.

This is the gene tree of one family (family 9 from the genome run in
``fig_species_tree_events``) evolving along the 10-species tree. The tips are gene
*copies*, each labelled by the species it lives in, so the event signatures are
visible directly:

  * duplication  -> two sister copies in the same species (I, I), hatched SQUARE at the node
  * transfer     -> a copy sitting next to its DONOR's lineage, not its own species
                    (a "G" copy is sister to F), hatched TRIANGLE at the node
  * loss         -> a dead-end lineage, drawn dashed and capped with a hatched CIRCLE
  * speciation   -> the remaining bifurcations (unmarked)

Same house style, dashing and hatched glyphs as the species-tree figures, so the
set reads as one story: the events of the previous figure produce this gene tree.

Run:  python figures/scripts/fig_gene_tree.py
"""

from __future__ import annotations

from pathlib import Path

import drawsvg as draw
import phylustrator as ph
from phylustrator.io import read_newick

from fig_species_tree_events import make_hatch
from fig_species_tree_extinct import annotate_depths, draw_skeleton, mark_survival
from zombi_style import INK, MUTED, species_style

FIG_DIR = Path(__file__).resolve().parent.parent
TREE_NWK = FIG_DIR / "gene_tree" / "gene_tree.nwk"
OUT_STEM = FIG_DIR / "gene_tree" / "gene_tree"

# species (ZOMBI node) -> display letter, matching the species-tree figures
SPECIES_LETTER = {"n9": "A", "n4": "B", "n10": "C", "n2": "D", "n5": "E",
                  "n1": "F", "n6": "G", "n3": "H", "n8": "I", "n7": "J"}

# reconciliation annotation for family 9 (from 9_events.tsv)
DUP_NODES = {"g477"}                                # duplication in species I
TRANSFER_INFO = {"g314": {"remaining": "n1_g554",   # transfer F -> G: donor copy stays (F),
                          "leaving": "n6_g555"}}     #   transferred copy leaves (to G)

MARKER_R = 8.0


def species_of(tip_name: str) -> str:
    """'n8_g686' -> 'I';  'LOSS_g137' -> '' (handled separately)."""
    if tip_name.startswith("LOSS_"):
        return ""
    return SPECIES_LETTER.get(tip_name.split("_")[0], tip_name.split("_")[0])


def main():
    tree = read_newick(TREE_NWK)
    present = annotate_depths(tree)
    mark_survival(tree, present, tol=1e-6 * present)
    name2node = {n.name: n for n in tree.traverse()}

    # Name tips like genes: <species>_<copy>, so duplicates read as I_1 / I_2, G_1 / G_2.
    # Blank the loss tip (marked only by the dead-end + circle glyph).
    loss_tips, copies = [], {}
    for lf in tree.get_leaves():
        if lf.name.startswith("LOSS_"):
            loss_tips.append(lf)
            lf.name = ""
        else:
            sp = species_of(lf.name)
            copies[sp] = copies.get(sp, 0) + 1
            lf.name = f"{sp}_{copies[sp]}"

    style = species_style(height=820)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()
    draw_skeleton(d, tree)               # solid survivors, dashed lost lineage
    d.add_leaf_names(color=INK, padding=12)

    hatch = make_hatch(d)
    glyph = dict(r=MARKER_R, stroke=INK, stroke_width=2.0)
    for gid in DUP_NODES:
        d._draw_shape_at(*name2node[gid].coordinates, "square", hatch, **glyph)
    for tnode, info in TRANSFER_INFO.items():
        node = name2node[tnode]
        d._draw_shape_at(*node.coordinates, "triangle", hatch, **glyph)
        xn = node.coordinates[0]
        for role, child in (("remaining branch", info["remaining"]),
                            ("leaving branch", info["leaving"])):
            cy = name2node[child].coordinates[1]
            d.add_text(role, xn + 15, cy - 9, font_size=13, color=MUTED)
    for lf in loss_tips:                 # hatched circle at the dead-end, no label
        d._draw_shape_at(*lf.coordinates, "circle", hatch, **glyph)

    ticks = [round(present * i / 4, 6) for i in range(5)]
    d.add_time_axis(ticks=ticks, tick_labels=[f"{t:.2f}" for t in ticks],
                    label="Time (root to present)", tick_size=6.0, padding=14.0,
                    stroke_width=1.6)
    add_legend(d, hatch, x=-style.width / 2 + 34, y=-style.height / 2 + 40)

    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    print(f"wrote {OUT_STEM}.svg / .png  ({len(tree.get_leaves())} gene copies)")


def add_legend(d, hatch, x, y, font_size=15):
    r = MARKER_R
    fam = d.style.font_family
    label_x = x + 2 * r + 14
    d.drawing.append(draw.Text("Gene-tree events", font_size + 2, x, y, font_weight="bold",
                               font_family=fam, text_anchor="start"))
    rows = [("square", "Duplication"), ("triangle", "Transfer"), ("circle", "Loss")]
    cy = y + font_size * 1.9
    for shape, label in rows:
        d._draw_shape_at(x + r, cy, shape, hatch, r=r, stroke=INK, stroke_width=2.0)
        d.drawing.append(draw.Text(label, font_size, label_x, cy, font_family=fam,
                                   text_anchor="start", dominant_baseline="middle"))
        cy += font_size * 1.9


if __name__ == "__main__":
    main()
