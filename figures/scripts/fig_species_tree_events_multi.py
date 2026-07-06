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


def draw_diamond(d, x, y, r, color):
    """A filled diamond marking an origination, in the family colour."""
    d.drawing.append(draw.Lines(x, y - r, x + r, y, x, y + r, x - r, y,
                                close=True, fill=color, stroke="none"))


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
    diamond per family just past the root stub so the colour<->family key is
    anchored on the tree itself."""
    root = tree.get_tree_root()
    x0, y0 = root.coordinates
    x = x0 - d.style.root_stub_length * 0.5
    span = (len(FAMILIES) - 1) * 15.0
    for i, fam in enumerate(FAMILIES):
        draw_diamond(d, x, y0 - span / 2 + i * 15.0, MARKER_R * 0.8, fam["color"])


def add_family_legend(d, x, y):
    """Colour -> family key (one row per family), a vertical column, top-left."""
    fam_font = d.style.font_family
    d.drawing.append(draw.Text("Gene family", FS_LABEL, x, y, font_weight="bold",
                               font_family=fam_font, text_anchor="start",
                               dominant_baseline="central", fill=INK))
    cy = y + FS_LABEL * 1.7
    for fam in FAMILIES:
        d._draw_shape_at(x + 9, cy, "square", fam["color"], r=9)
        d.drawing.append(draw.Text(fam["label"], FS_LABEL, x + 9 + 22, cy,
                                   font_family=fam_font, text_anchor="start",
                                   dominant_baseline="central", fill=INK))
        cy += FS_LABEL * 1.7


def add_event_legend(d, x, y):
    """Shape -> event key (grey, colour-neutral glyphs), a vertical column.
    Order: Origination, Duplication, Transfer, Loss."""
    fam_font = d.style.font_family
    key = INK
    r = 9
    gap = 22
    row = 36
    d.drawing.append(draw.Text("Event", FS_LABEL, x, y, font_weight="bold",
                               font_family=fam_font, text_anchor="start",
                               dominant_baseline="central", fill=INK))
    cy = y + FS_LABEL * 1.7

    # origination: diamond
    draw_diamond(d, x + r, cy, r * 0.85, key)
    d.drawing.append(draw.Text("Origination", FS_LABEL, x + r + gap, cy, font_family=fam_font,
                               text_anchor="start", dominant_baseline="central", fill=INK))
    # duplication: filled square
    cy += row
    d._draw_shape_at(x + r, cy, "square", key, r=r)
    d.drawing.append(draw.Text("Duplication", FS_LABEL, x + r + gap, cy, font_family=fam_font,
                               text_anchor="start", dominant_baseline="central", fill=INK))
    # transfer: donor dot -> arc -> arrowhead
    cy += row
    d.drawing.append(draw.Circle(x + r - 9, cy, 3.2, fill=key))
    p = draw.Path(stroke=key, stroke_width=2.6, fill="none", stroke_linecap="round")
    p.M(x + r - 9, cy).C(x + r - 9, cy - 12, x + r + 9, cy - 12, x + r + 9, cy)
    d.drawing.append(p)
    d.drawing.append(draw.Lines(x + r + 9, cy + 1, x + r + 9 - 5, cy - 8, x + r + 9 + 5, cy - 8,
                                close=True, fill=key))
    d.drawing.append(draw.Text("Transfer", FS_LABEL, x + r + gap, cy, font_family=fam_font,
                               text_anchor="start", dominant_baseline="central", fill=INK))
    # loss: cross
    cy += row
    draw_cross(d, x + r, cy, r, key)
    d.drawing.append(draw.Text("Loss", FS_LABEL, x + r + gap, cy, font_family=fam_font,
                               text_anchor="start", dominant_baseline="central", fill=INK))


def main():
    tree = read_newick(str(TREE_NWK))
    present = annotate_times(tree)

    # wide landscape canvas: three families' events need horizontal room and a
    # left header band clear of the crown for the two legends.
    style = species_style(width=1260, height=820, margin=120, font_size=FS_TICK)
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

    left = -style.width / 2 + 40
    d.drawing.append(draw.Text("Three gene families on one species tree", FS_TITLE, 0,
                               -style.height / 2 + 42, font_weight="bold",
                               font_family=style.font_family, text_anchor="middle",
                               dominant_baseline="central", fill=INK))
    add_family_legend(d, x=left, y=-style.height / 2 + 96)
    add_event_legend(d, x=left + 210, y=-style.height / 2 + 96)

    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    print(f"wrote {OUT_STEM}.svg / .png")


if __name__ == "__main__":
    main()
