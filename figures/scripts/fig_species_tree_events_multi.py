"""Figure: the events of THREE gene families painted on one species tree.

The single-family companion (``fig_species_tree_events``) shows one family's
duplication/transfer/loss on the tree. This wider figure overlays *three*
families at once, each drawn in its own colour, to make the point that a genome
is a whole population of families threaded through the same species tree by one
shared Gillespie process — every branch carries the events of many families.

The tree and the events are REAL ZOMBI2 output (not hand-placed): a 12-tip
backward birth-death tree (seed 6) run through ``simulate_genomes`` (seed 4),
from which three families with a clean one-duplication / one-transfer / one-loss
signature are selected and painted at their true simulated times.

Encoding
--------
  * colour            -> which gene family (Paul-Tol categorical set; the
                         STYLE.md categorical exception, a small identity set)
  * filled square     -> duplication
  * cross (x)         -> loss
  * arc + arrowhead   -> transfer (donor dot -> arrowhead into recipient)
  * small diamond     -> origination (all three originate at the root here)

Monochrome would collapse the three families together, so colour here carries
information (family identity), exactly the case STYLE.md allows.

Run:  python figures/scripts/fig_species_tree_events_multi.py
"""

from __future__ import annotations

from pathlib import Path

import drawsvg as draw
import phylustrator as ph
from phylustrator.io import read_newick

from zombi_style import (
    INK, PANEL, MODULE_COLORS, species_style,
    FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK,
)

FIG_DIR = Path(__file__).resolve().parent.parent
TREE_NWK = FIG_DIR / "species_tree_events_multi" / "species_tree.nwk"
OUT_STEM = FIG_DIR / "species_tree_events_multi" / "species_tree_events_multi"

# --- The three families and their events -----------------------------------
# Real ZOMBI2 output: BirthDeath(birth=1.0, death=0.3) backward tree, n_tips=12,
# seed 6; simulate_genomes(dup=0.35, trans=0.25, loss=0.35, orig=0.5,
# initial_families=25, seed=4). Families 7, 9, 15 each have exactly one
# duplication, one transfer and one loss (plus the root origination). See
# select_ch8.py / fig_species_tree_events_multi header for reproduction.
FAMILIES = [
    {   # family 7
        "label": "Family 1",
        "color": MODULE_COLORS[0],                       # blue
        "duplications": [("n7", 0.928)],
        "losses":       [("n6", 0.266)],
        "transfers":    [{"from": "n4", "to": "n5", "time": 0.859}],
    },
    {   # family 9
        "label": "Family 2",
        "color": MODULE_COLORS[1],                       # orange
        "duplications": [("n12", 0.939)],
        "losses":       [("n3", 0.630)],
        "transfers":    [{"from": "n4", "to": "n9", "time": 0.878}],
    },
    {   # family 15
        "label": "Family 3",
        "color": MODULE_COLORS[2],                       # purple
        "duplications": [("n1", 0.941)],
        "losses":       [("n4", 0.651)],
        "transfers":    [{"from": "i3", "to": "n6", "time": 0.498}],
    },
]

MARKER_R = 8.5


def draw_cross(d, x, y, r, color, stroke_width=3.0):
    """A cross (two crossing segments) marking a loss, in the family colour."""
    d.drawing.append(draw.Line(x - r, y - r, x + r, y + r, stroke=color,
                               stroke_width=stroke_width, stroke_linecap="round"))
    d.drawing.append(draw.Line(x - r, y + r, x + r, y - r, stroke=color,
                               stroke_width=stroke_width, stroke_linecap="round"))


def draw_origin(d, x, y, r, color):
    """A filled circle marking an origination, in the family colour."""
    d.drawing.append(draw.Circle(x, y, r, fill=color, stroke="none"))


def annotate_times(tree) -> float:
    """Set ``node.time_from_origin`` (root = 0); return the present (max)."""
    present = 0.0
    for n in tree.traverse("preorder"):
        n.time_from_origin = 0.0 if n.is_root() else n.up.time_from_origin + n.dist
        present = max(present, n.time_from_origin)
    return present


def _branch_point(node, t):
    """(x, y) at absolute time ``t`` along the branch above ``node``."""
    p = node.up
    x_c, y_c = node.coordinates
    if p is None:
        return x_c, y_c
    x_p, _ = p.coordinates
    span = node.time_from_origin - p.time_from_origin
    frac = 0.5 if abs(span) < 1e-12 else max(0.0, min(1.0, (t - p.time_from_origin) / span))
    return x_p + (x_c - x_p) * frac, y_c


def draw_family_events(d, tree, fam):
    """Paint one family's duplication/loss/transfer glyphs in its colour."""
    name2node = {n.name: n for n in tree.traverse()}
    c = fam["color"]
    for br, t in fam["duplications"]:
        d._draw_shape_at(*_branch_point(name2node[br], t), "square", c, r=MARKER_R)
    for br, t in fam["losses"]:
        x, y = _branch_point(name2node[br], t)
        draw_cross(d, x, y, MARKER_R, c)
    d.plot_transfers(fam["transfers"], mode="time", use_gradient=False, color=c,
                     stroke_width=2.6, arc_intensity=42.0, opacity=1.0,
                     arrowhead=True, arrow_size=12.0, donor_dot=True)


def draw_originations(d, tree):
    """All three families originate at the root at t=0; stack a small coloured
    circle per family just past the root stub so the colour<->family key is
    anchored on the tree itself."""
    root = tree.get_tree_root()
    x0, y0 = root.coordinates
    x = x0 - d.style.root_stub_length * 0.5
    span = (len(FAMILIES) - 1) * 16.0
    for i, fam in enumerate(FAMILIES):
        draw_origin(d, x, y0 - span / 2 + i * 16.0, MARKER_R * 0.8, fam["color"])


def _legend_transfer_glyph(d, x, cy, key, r):
    """A miniature donor-dot -> arc -> arrowhead, centred on (x, cy)."""
    d.drawing.append(draw.Circle(x - r, cy, 3.4, fill=key))
    p = draw.Path(stroke=key, stroke_width=2.6, fill="none", stroke_linecap="round")
    p.M(x - r, cy).C(x - r, cy - 11, x + r, cy - 11, x + r, cy)
    d.drawing.append(p)
    d.drawing.append(draw.Lines(x + r, cy + 1, x + r - 5, cy - 7, x + r + 5, cy - 7,
                                close=True, fill=key))


def add_legends(d, x, y):
    """Two clean legend columns in the top-left header band, well clear of the
    tree: a colour->family key and a glyph->event key. Each column is a heading
    plus one row per entry; the two columns are widely separated so neither text
    block runs into the other or into the crown.

    Columns (glyph centred on ``gx``; label starts at ``gx + LGAP``):
      * family:  coloured square + family name
      * event:   circle(origination) / square(dup) / arc(transfer) / cross(loss)
    """
    fam_font = d.style.font_family
    key = INK
    R = 10                       # legend glyph radius
    ROW = FS_LABEL * 1.55        # vertical spacing between rows
    LGAP = 30                    # glyph centre -> label start
    COL2 = 300                   # family column -> event column offset

    def heading(hx, text):
        d.drawing.append(draw.Text(text, FS_LABEL, hx, y, font_weight="bold",
                                   font_family=fam_font, text_anchor="start",
                                   dominant_baseline="central", fill=INK))

    def label(lx, ly, text):
        d.drawing.append(draw.Text(text, FS_LABEL, lx, ly, font_family=fam_font,
                                   text_anchor="start", dominant_baseline="central",
                                   fill=INK))

    # --- column 1: gene family (colour key) ---
    heading(x, "Gene family")
    gx = x + R
    cy = y + ROW
    for fam in FAMILIES:
        d._draw_shape_at(gx, cy, "square", fam["color"], r=R)
        label(gx + LGAP, cy, fam["label"])
        cy += ROW

    # --- column 2: event (glyph key, colour-neutral) ---
    ex = x + COL2
    heading(ex, "Event")
    gx = ex + R
    cy = y + ROW
    # origination: filled circle
    draw_origin(d, gx, cy, R, key)
    label(gx + LGAP, cy, "Origination")
    # duplication: filled square
    cy += ROW
    d._draw_shape_at(gx, cy, "square", key, r=R)
    label(gx + LGAP, cy, "Duplication")
    # transfer: donor dot -> arc -> arrowhead
    cy += ROW
    _legend_transfer_glyph(d, gx, cy, key, R)
    label(gx + LGAP, cy, "Transfer")
    # loss: cross
    cy += ROW
    draw_cross(d, gx, cy, R, key)
    label(gx + LGAP, cy, "Loss")


def main():
    tree = read_newick(str(TREE_NWK))
    present = annotate_times(tree)

    # wide landscape canvas. The generous top margin opens a clear header band
    # above the crown that holds the title and both legends, so nothing overlaps
    # a branch; the tree still fills the lower ~two-thirds and stays landscape.
    style = species_style(width=1300, height=1030, margin=285, font_size=FS_TICK)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d.draw()
    d.add_leaf_names(color=INK, padding=12)

    draw_originations(d, tree)
    for fam in FAMILIES:
        draw_family_events(d, tree, fam)

    ticks = [round(present * i / 4, 6) for i in range(5)]
    d.add_time_axis(ticks=ticks, tick_labels=[f"{t:.2f}" for t in ticks],
                    label="Time (root to present)", tick_size=6.0, padding=16.0,
                    stroke_width=1.6)

    left = -style.width / 2 + 46
    d.drawing.append(draw.Text("Three gene families on one species tree", FS_TITLE, 0,
                               -style.height / 2 + 44, font_weight="bold",
                               font_family=style.font_family, text_anchor="middle",
                               dominant_baseline="central", fill=INK))
    add_legends(d, x=left, y=-style.height / 2 + 100)

    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    print(f"wrote {OUT_STEM}.svg / .png")


if __name__ == "__main__":
    main()
