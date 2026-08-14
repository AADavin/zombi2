"""Figure: a complete species tree with the extinct lineages still on it.

The forward birth-death process keeps the lineages that died out. Drawing the *complete* tree puts
the surviving skeleton and the extinct history in one picture, which is the point that motivates
ZOMBI2's forward engine: a backward simulation only ever sees the solid tree, and the dashes are
what really happened.

A branch is drawn **dashed** when the node it leads to has no surviving descendant. Applied to every
branch that gives exactly the three cases you want — a fork whose daughters both survive is solid, a
fork that loses one daughter is half dashed, and a speciation inside a dead clade is fully dashed.
Which lineages survived is read off the tree geometrically (a leaf whose root-distance falls short of
the present), so this works on any complete-tree Newick.

No legend inside the figure: solid-surviving / dashed-extinct is stated in the chapter's caption,
which is where this book says a figure's words belong.

This module also exports `dead_ends` and `present_of`, for any figure that has to tell a
survivor from a lineage that died.

Run:  python figures/scripts/fig_species_tree_extinct.py
"""

from __future__ import annotations

from pathlib import Path

import phylustrator as ph

from zombi_style import save, FS_TICK, tree_style

FIG_DIR = Path(__file__).resolve().parent.parent
TREE_NWK = FIG_DIR / "species_tree_extinct" / "species_tree.nwk"

# ASCII-only names for the extinct tips; survivors are left unlabelled so a dense tree stays legible.
EXTINCT_NAMES = [f"e{i}" for i in range(1, 25)]


def present_of(tree) -> float:
    """The distance from the origin to the furthest tip — where "today" is."""
    return max(tree.depth(node) for node in tree.walk())


def dead_ends(tree, survivor) -> list:
    """Every node with **no** surviving descendant, as nodes rather than names.

    `plot(dashed=…)` matches by name, so a caller that renames its tips has to take the names
    *after* renaming; handing back nodes is what makes that impossible to get wrong.
    ``survivor(leaf)`` says whether one leaf reached the present.
    """
    alive: dict[int, bool] = {}
    for node in reversed(list(tree.walk())):        # children before parents
        alive[id(node)] = (survivor(node) if not node.children
                           else any(alive[id(c)] for c in node.children))
    return [node for node in tree.walk() if not alive[id(node)]]


def reaches_present(tree, present: float, tol: float | None = None):
    """A `survivor` test for `dead_ends`: a leaf survives when its root-distance reaches the present.

    `to_newick` prints six significant figures and that rounding accumulates along a deep path, so
    the test needs a tolerance rather than equality. A caller that knows which tips survived — a run
    in hand rather than a Newick file — should test that instead.
    """
    tol = present * 1e-6 if tol is None else tol
    return lambda leaf: abs(tree.depth(leaf) - present) <= tol


def render() -> None:
    tree = ph.trees.read(TREE_NWK)
    present = present_of(tree)
    dead = dead_ends(tree, reaches_present(tree, present))

    # `plot(dashed=…)` and `tip_labels()` both work off `node.name`, so every node that has to be
    # dashed needs a name of its own before the surviving tips are blanked — otherwise a nameless
    # internal branch and every unlabelled survivor share the name "" and all of them go dashed.
    for i, node in enumerate(tree.walk()):
        if not node.name:
            node.name = f"_i{i}"

    extinct = [node for node in dead if not node.children]
    for leaf in tree.leaves:                       # only the extinct tips are named: they are the point
        leaf.name = ""
    for leaf, name in zip(extinct, EXTINCT_NAMES):
        leaf.name = name
    dashed = {node.name for node in dead}          # names last: `plot(dashed=…)` matches by name

    figure = (ph.trees.plot(tree, dashed=dashed, style=tree_style(1320, 820))
              + ph.trees.tip_labels(size=FS_TICK)
              + ph.trees.time_axis("time (origin to present)"))
    save(figure, "species_tree")
    print(f"  ({len(tree.leaves) - len(extinct)} extant + {len(extinct)} extinct tips)")


if __name__ == "__main__":
    render()
