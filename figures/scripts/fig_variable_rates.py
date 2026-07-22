"""Figure: the three ways a rate can vary, one tree apiece.

Chapter 3 says a birth or death rate can depend on **time**, on **how crowded the tree is**, or on a
lineage's **ancestry**, and that each is written the same way — a base rate times a modifier. Three
trees stacked show what each one does to the shape of a run:

  A  ``OnTime``             the rate changes at set times, so the tree grows fast and then crawls
  B  ``OnTotalDiversity``   the rate slows as the tree fills, so diversity levels off
  C  ``FromParent``         each lineage inherits its parent's rate, so clades run at their own tempo

The trees are **simulated here**, by the engine the chapter documents, rather than drawn by hand —
so the figure cannot drift away from what the code does. Each panel stops at the same 25 surviving
lineages, which is what makes the three comparable: the same amount of tree, reached differently.
Bounding by tip count rather than by time also keeps the drifting rate in panel C from running away.

The dashed-extinct skeleton is shared with ``fig_species_tree_extinct.py``. No legend here: the
figure just above this one in the chapter defines solid-surviving / dashed-extinct, and a legend
box inside a panel this wide has nowhere to sit that is clear of the tree.

House style: B&W, ASCII text. No title inside the figure — the manual captions it.

Run:  python figures/scripts/fig_variable_rates.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import drawsvg as draw
import phylustrator as ph
from phylustrator.io import read_newick

from fig_species_tree_extinct import PRESENT, annotate_depths, draw_skeleton
from zombi2 import species
from zombi2.rates import modifiers as mod
from zombi_style import save, INK, MUTED, species_style, FS_LABEL, FS_TICK, FS_TITLE

N_EXTANT = 25          # every panel stops at the same standing diversity
DEATH = 0.1            # low, so the dashed extinct lineages stay a garnish rather than the picture
SEED = 3               # a seed where the three come out at comparable depths (~4.4-5.0)

#: (panel letter, the rate as it is written in the chapter, one line on what it does)
PANELS = [
    ("A", "birth = 1.2 * OnTime({0: 1.0, 2.0: 0.3})",
     1.2 * mod.OnTime({0: 1.0, 2.0: 0.3}),
     "the rate drops to a third at time 2, so an early burst gives way to a slow tail",
     2.0),                 # the skyline breakpoint, marked on the panel
    ("B", "birth = 1.2 * OnTotalDiversity(cap=30)",
     1.2 * mod.OnTotalDiversity(cap=30),
     "the rate falls as the tree fills toward its cap, so splits thin out near the present",
     None),                # nothing to mark: the rate falls continuously, not at a moment
    ("C", "birth = 0.45 * FromParent(spread=0.5)",
     0.45 * mod.FromParent(spread=0.5),
     "each lineage inherits its parent's rate, so a clade keeps its own tempo",
     None),
]

PANEL_W, PANEL_H = 1180, 420


def panel(letter: str, written: str, birth, gloss: str, mark_time: float | None):
    """One panel as a drawsvg group: a tree grown under ``birth``, drawn to this figure's style."""
    result = species.simulate_species_tree(birth=birth, death=DEATH, n_extant=N_EXTANT, seed=SEED)
    # phylustrator reads a Newick *file*, so the simulated tree goes through one
    with tempfile.NamedTemporaryFile("w", suffix=".nwk", delete=False) as fh:
        fh.write(result.complete_tree.to_newick() + "\n")
    try:
        tree = read_newick(fh.name)
    finally:
        Path(fh.name).unlink()
    present = annotate_depths(tree)
    # Which lineages survived comes from the run, not from comparing depths: to_newick prints six
    # significant figures, and that rounding accumulates along a root-to-tip path to more than the
    # depth test's tolerance, so a survivor at the end of a deep path gets called extinct and drawn
    # dashed all the way to the present. The labels are n<id>, so the answer can just be looked up.
    extant_ids = {f"n{node.id}" for node in result.complete_tree.extant()}
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node.is_extant = node.name in extant_ids
            node.has_survivor = node.is_extant
        else:
            node.has_survivor = any(c.has_survivor for c in node.children)

    for leaf in tree.get_leaves():          # a dense tree reads better with bare tips
        leaf.name = ""

    style = species_style(width=PANEL_W, height=PANEL_H, margin=96, font_size=FS_TICK)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()

    ys = [leaf.y_coord for leaf in tree.get_leaves()]
    present_x = d.root_x + present * d.sf
    d.drawing.append(draw.Line(present_x, min(ys) - 12, present_x, max(ys) + 12,
                               stroke=PRESENT, stroke_width=1.0, stroke_dasharray="2,4"))
    if mark_time is not None:
        # The moment the rate changes, drawn like the faint "present" rule so it reads as a
        # reference rather than as part of the tree. It carries no label of its own — it lands on
        # its own axis tick below, which says the time and cannot collide with anything.
        mx = d.root_x + mark_time * d.sf
        d.drawing.append(draw.Line(mx, min(ys) - 12, mx, max(ys) + 12, stroke=MUTED,
                                   stroke_width=1.4, stroke_dasharray="5,4"))

    draw_skeleton(d, tree)

    ticks = [round(present * i / 4, 6) for i in range(5)]
    if mark_time is not None:                # give the breakpoint a tick of its own, and drop any
        keep = [x for x in ticks if abs(x - mark_time) > present * 0.15]   # tick it would crowd
        ticks = sorted(keep + [mark_time])
    d.add_time_axis(ticks=ticks, tick_labels=[f"{t:.1f}" for t in ticks],
                    label="time (origin to present)", tick_size=6.0, padding=12.0,
                    stroke_width=1.6)

    left, top = -style.width / 2, -style.height / 2
    d.drawing.append(draw.Text(letter, FS_TITLE, left + 16, top + 46, font_weight="bold",
                               font_family=style.font_family, text_anchor="start",
                               dominant_baseline="central", fill=INK))
    d.drawing.append(draw.Text(written, FS_LABEL, left + 62, top + 40, font_weight="bold",
                               font_family=style.font_family, text_anchor="start",
                               dominant_baseline="central", fill=INK))
    d.drawing.append(draw.Text(gloss, FS_LABEL, left + 62, top + 72, font_style="italic",
                               font_family=style.font_family, text_anchor="start",
                               dominant_baseline="central", fill=MUTED))
    # the drawer works in centred coordinates; shift into this panel's slot in the stack
    group = draw.Group(transform=f"translate({PANEL_W / 2},{PANEL_H / 2})")
    for element in d.drawing.elements:
        group.append(element)
    n_extinct = len([leaf for leaf in tree.get_leaves() if not leaf.is_extant])
    return group, n_extinct


def render() -> None:
    d = draw.Drawing(PANEL_W, PANEL_H * len(PANELS), origin=(0, 0))
    d.append(draw.Rectangle(0, 0, PANEL_W, PANEL_H * len(PANELS), fill="white"))

    for row, (letter, written, birth, gloss, mark) in enumerate(PANELS):
        group, n_extinct = panel(letter, written, birth, gloss, mark)
        slot = draw.Group(transform=f"translate(0,{row * PANEL_H})")
        slot.append(group)
        d.append(slot)
        print(f"  {letter}: {N_EXTANT} extant + {n_extinct} extinct")

    save(d, "variable_rates")


if __name__ == "__main__":
    render()
