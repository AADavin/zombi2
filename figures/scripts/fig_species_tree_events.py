"""Figure: gene-family events painted on the species tree.

This is the *scenario* companion to the gene-tree figure: the exact events of
ZOMBI2 gene family 9 (genome run along this tree, dup=trans=loss=0.2, seed 7),
placed at their true time on the species tree. The same three events reappear as
the gene tree in ``fig_gene_tree`` — species tree -> events -> gene tree.

  * duplication in species I           -> hatched square
  * transfer F -> G                    -> black arc, donor dot -> arrowhead on recipient
  * loss in species J                  -> hatched circle

Monochrome, print-friendly. Run:  python figures/scripts/fig_species_tree_events.py
"""

from __future__ import annotations

import string
from pathlib import Path

import drawsvg as draw
import phylustrator as ph
from phylustrator.io import read_newick

from zombi_style import INK, PANEL, species_style, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

FIG_DIR = Path(__file__).resolve().parent.parent
TREE_NWK = FIG_DIR / "species_tree_10" / "species_tree.nwk"
OUT_STEM = FIG_DIR / "species_tree_events" / "species_tree_events"

# Family 9's events, in DISPLAY-name space (leaves A-J), at their true simulated
# times (from 9_events.tsv) so this figure matches the gene tree exactly.
DUPLICATIONS = [("I", 0.5965)]                              # duplication in species I
LOSSES = [("J", 0.7128)]                                    # loss in species J
TRANSFERS = [{"from": "F", "to": "G", "time": 0.4233}]      # transfer F -> G

MARKER_R = 9.0


def draw_cross(d, x, y, r, stroke_width=3.0):
    """A solid black cross (two crossing segments) marking a loss."""
    d.drawing.append(draw.Line(x - r, y - r, x + r, y + r, stroke=INK,
                               stroke_width=stroke_width, stroke_linecap="round"))
    d.drawing.append(draw.Line(x - r, y + r, x + r, y - r, stroke=INK,
                               stroke_width=stroke_width, stroke_linecap="round"))


def annotate_times(tree) -> float:
    """Set ``node.time_from_origin`` (root = 0); return the present (max)."""
    present = 0.0
    for n in tree.traverse("preorder"):
        n.time_from_origin = 0.0 if n.is_root() else n.up.time_from_origin + n.dist
        present = max(present, n.time_from_origin)
    return present


def make_hatch(d, spacing=4.5, sw=1.2):
    """A diagonal black-on-white hatch pattern, shared by duplication & loss glyphs."""
    p = draw.Pattern(spacing, spacing, id="event_hatch", patternUnits="userSpaceOnUse")
    p.append(draw.Rectangle(0, 0, spacing, spacing, fill=PANEL))
    for off in (-spacing, 0, spacing):
        p.append(draw.Line(off, spacing, off + spacing, 0, stroke=INK, stroke_width=sw))
    d.drawing.append(p)
    return p


def _branch_point(d, node, t):
    """(x, y) at absolute time ``t`` along the branch above ``node``."""
    p = node.up
    x_c, y_c = node.coordinates
    if p is None:
        return x_c, y_c
    x_p, _ = p.coordinates
    span = node.time_from_origin - p.time_from_origin
    frac = 0.5 if abs(span) < 1e-12 else max(0.0, min(1.0, (t - p.time_from_origin) / span))
    return x_p + (x_c - x_p) * frac, y_c


def draw_events(d, tree):
    """Solid-black event glyphs: filled square = duplication, cross = loss,
    black arc with arrowhead = transfer."""
    name2node = {n.name: n for n in tree.traverse()}
    for br, t in DUPLICATIONS:
        d._draw_shape_at(*_branch_point(d, name2node[br], t), "square", INK, r=MARKER_R)
    for br, t in LOSSES:
        x, y = _branch_point(d, name2node[br], t)
        draw_cross(d, x, y, MARKER_R)
    d.plot_transfers(TRANSFERS, mode="time", use_gradient=False, color=INK,
                     stroke_width=2.6, arc_intensity=38.0, opacity=1.0,
                     arrowhead=True, arrow_size=12.0, donor_dot=True)


def add_legend(d, x, y):
    """Symbol legend (solid-black glyphs), a single vertical column in the top-left:
    one row per event, ordered Duplication, Transfer, Loss."""
    r = MARKER_R
    fam = d.style.font_family
    gap = 24                                     # symbol-to-label gap
    row = 40                                      # vertical spacing between rows

    # duplication: filled square
    cy = y
    d._draw_shape_at(x, cy, "square", INK, r=r)
    d.drawing.append(draw.Text("Duplication", FS_LABEL, x + r + gap, cy, font_family=fam,
                               text_anchor="start", dominant_baseline="central", fill=INK))

    # transfer: donor dot -> black arc -> arrowhead
    cy = y + row
    d.drawing.append(draw.Circle(x - r, cy, 3.2, fill=INK))                             # donor dot
    p = draw.Path(stroke=INK, stroke_width=2.6, fill="none", stroke_linecap="round")
    p.M(x - r, cy).C(x - r, cy - 12, x + r, cy - 12, x + r, cy)
    d.drawing.append(p)
    d.drawing.append(draw.Lines(x + r, cy + 1, x + r - 5, cy - 8, x + r + 5, cy - 8,
                                close=True, fill=INK))                                  # arrowhead
    d.drawing.append(draw.Text("Transfer", FS_LABEL, x + r + gap, cy, font_family=fam,
                               text_anchor="start", dominant_baseline="central", fill=INK))

    # loss: black cross
    cy = y + 2 * row
    draw_cross(d, x, cy, r)
    d.drawing.append(draw.Text("Loss", FS_LABEL, x + r + gap, cy, font_family=fam,
                               text_anchor="start", dominant_baseline="central", fill=INK))


def main():
    tree = read_newick(TREE_NWK)
    present = annotate_times(tree)
    for leaf, letter in zip(tree.get_leaves(), string.ascii_uppercase):
        leaf.name = letter

    # taller canvas + generous top margin gives a clean header band (title + legend row)
    # that sits entirely above the tree.
    style = species_style(width=920, height=760, margin=120, font_size=FS_TICK)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d.draw()
    d.add_leaf_names(color=INK, padding=12)

    draw_events(d, tree)

    ticks = [round(present * i / 4, 6) for i in range(5)]
    d.add_time_axis(ticks=ticks, tick_labels=[f"{t:.2f}" for t in ticks],
                    label="Time (root to present)", tick_size=6.0, padding=16.0,
                    stroke_width=1.6)

    # title centered at the top; symbol legend as a vertical column in the top-left,
    # both clear of the tree
    left = -style.width / 2 + 30
    d.drawing.append(draw.Text("Gene-family events on the species tree", FS_TITLE, 0,
                               -style.height / 2 + 42, font_weight="bold",
                               font_family=style.font_family, text_anchor="middle",
                               dominant_baseline="central", fill=INK))
    add_legend(d, x=left + 12, y=-style.height / 2 + 88)

    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    print(f"wrote {OUT_STEM}.svg / .png")


if __name__ == "__main__":
    main()
