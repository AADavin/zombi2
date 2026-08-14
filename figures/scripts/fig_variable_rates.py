"""Figure: the four ways a rate can vary, one tree apiece.

Chapter 3 says a birth or death rate can depend on **time**, on **how crowded the tree is**, or on a
lineage's **ancestry**, and that each is written the same way — a scope holding the base rate, with a
verb chained onto it. Four trees stacked show what each one does to the shape of a run:

  A  ``changing_at``
       the rate changes at set times: an early burst, then a slow tail
  B  ``scaled_by(TotalDiversity(...))``
       the rate slows as the tree fills, so diversity levels off
  C  ``varying_among('lineages', Drift(...))``
       each lineage inherits its parent's rate, so clades run at their own tempo
  D  ``varying_among('lineages', LogNormal(...))``
       each lineage draws its own rate, with no memory of its parent

One panel per bullet of the chapter's section, which is what makes the figure readable beside it —
C and D are the pair that reads as one bullet and is not: the same distribution, inherited or not.

The trees are **simulated here**, by the engine the chapter documents, rather than drawn by hand —
so the figure cannot drift away from what the code does. Each panel stops at the same 25 surviving
lineages, which is what makes them comparable: the same amount of tree, reached differently.
Bounding by tip count rather than by time also keeps the drifting rate in panel C from running away.

Each panel carries its letter and a plain-English label, the chapter's own name for that kind of
variation. No expression is printed: the section's prose carries no code either, and the written
form of each rate is in the chapter's Literature table, so the figure would be the third place to
say the same thing. The caption stays to one line. Solid-surviving / dashed-extinct is the
convention of the figure just above this one in the chapter, and is likewise stated there.

The panels are stacked as nested `<svg>` elements, which is how one SVG holds several
independent coordinate systems; `phylustrator` draws one figure at a time.

House style: B&W, ASCII text. No title inside the figure — the manual captions it.

Run:  python figures/scripts/fig_variable_rates.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import phylustrator as ph

from fig_species_tree_extinct import dead_ends
from zombi2 import species
from zombi2.params import Drift, LogNormal, PerLineage, TotalDiversity
from zombi_style import save, tree_style, FS_LABEL

N_EXTANT = 25          # every panel stops at the same standing diversity
DEATH = 0.1            # low, so the dashed extinct lineages stay a garnish rather than the picture
SEED = 3               # a seed where the four come out at comparable depths (~3.9-5.0)

#: (panel letter, what the panel says, the rate itself, a time to mark or None)
#:
#: The label is the chapter's own name for that kind of variation, in words. It used to be the rate's
#: `repr`, which could not drift from the code but put an expression in a figure whose section
#: deliberately carries none; the written form of each is in the chapter's Literature table. Keep a
#: label saying the same thing as the rate beside it — nothing checks that they agree.
PANELS = [
    ("A", "on time: the rate drops at a set moment",
     PerLineage(1.2).changing_at({0: 1.0, 2.0: 0.3}),
     2.0),                 # the skyline breakpoint, marked on the panel
    ("B", "on total diversity: the rate slows as the tree fills",
     PerLineage(1.2).scaled_by(TotalDiversity(cap=30)),
     None),                # nothing to mark: the rate falls continuously, not at a moment
    ("C", "on the parent's rate: each lineage inherits and drifts from it",
     PerLineage(0.45).varying_among('lineages', Drift(LogNormal(0.0, 0.5))),
     None),
    # D is C's counterpart without the inheritance: the same distribution, drawn afresh per lineage.
    # Its base rate is lower than A and B's because the multipliers have mean 1 either way, and
    # higher than C's because drift accumulates down a path and C's does not need the help.
    ("D", "by lineage: each draws its own rate, with no memory of its parent",
     PerLineage(0.85).varying_among('lineages', LogNormal(0.0, 0.5)),
     None),
]

PANEL_W, PANEL_H = 1180, 420


def panel(letter: str, says: str, birth, mark_time: float | None) -> tuple[str, int]:
    """One panel's SVG: a tree grown under ``birth``, drawn to this figure's style."""
    result = species.simulate_species_tree(birth=birth, death=DEATH, n_extant=N_EXTANT, seed=SEED)
    tree = ph.trees.loads(result.complete_tree.to_newick() + ";")

    # Which lineages survived comes from the run, not from comparing depths: to_newick prints six
    # significant figures, and that rounding accumulates along a root-to-tip path to more than a
    # depth test's tolerance, so a survivor at the end of a deep path would be called extinct and
    # drawn dashed all the way to the present. The labels are n<id>, so the answer is looked up.
    extant = {f"n{node.id}" for node in (result.complete_tree.nodes[_i] for _i in result.complete_tree.extant_leaves())}
    dashed = {node.name for node in dead_ends(tree, lambda leaf: leaf.name in extant)}
    # The names stay as the run wrote them: this figure adds no `tip_labels` layer, so a dense tree
    # comes out bare anyway, and blanking them would collapse every tip onto one name.

    figure = (ph.trees.plot(tree, dashed=dashed, style=tree_style(PANEL_W, PANEL_H, margin=96))
              + ph.trees.note(f"{letter}    {says}", loc="top-left", size=FS_LABEL, dy=-16)
              + ph.trees.time_axis("time (origin to present)"))
    if mark_time is not None:
        # The moment the rate changes, drawn as a faint rule so it reads as a reference rather than
        # as part of the tree. Unlabelled: a marker label is drawn at the top of the panel, where
        # the rate already sits, and the caption says which time this is.
        figure = figure + ph.trees.time_marker(mark_time)
    return figure.as_svg(), len(dashed)


def stack(panels: list[str], width: int, height: int) -> str:
    """Put each panel's SVG in its own row of one document, as a nested `<svg>` with a y offset."""
    rows = []
    for i, svg in enumerate(panels):
        inner = re.sub(r"^<\?xml[^>]*\?>\s*", "", svg.strip())   # one XML declaration per document
        inner = inner.replace("<svg ", f'<svg x="0" y="{i * height}" ', 1)
        rows.append(inner)
    body = "\n".join(rows)
    return (f'<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
            f'width="{width}" height="{height * len(panels)}" '
            f'viewBox="0 0 {width} {height * len(panels)}">\n'
            f'<rect x="0" y="0" width="{width}" height="{height * len(panels)}" fill="white"/>\n'
            f"{body}\n</svg>\n")


def render() -> None:
    svgs = []
    for letter, says, birth, mark in PANELS:
        svg, n_dead = panel(letter, says, birth, mark)
        svgs.append(svg)
        print(f"  {letter}: {N_EXTANT} extant, {n_dead} nodes with no survivor")
    save(stack(svgs, PANEL_W, PANEL_H), "variable_rates")


if __name__ == "__main__":
    render()
