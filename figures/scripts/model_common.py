"""Shared helpers for the species-tree *model* figures (FBD, episodic, ghosts, ...).

These models produce tips of several kinds (extant / extinct / fossil sample /
sampled ancestor / unsampled 'ghost'), which the Newick can't carry — so we
simulate with the ZOMBI2 Python API and convert the resulting ``Tree`` (which
keeps ``is_extant`` / ``sampled`` flags) to an ete3 tree for Phylustrator.

Convention shared across the model figures (monochrome):
  * a lineage is drawn SOLID if it leads to something we *observe* (an extant
    sampled tip or a fossil sample), and DASHED if its whole subtree is
    unobserved (extinct / unsampled 'ghost') — reusing the Fig-2 dashing.
  * a fossil sample (sampled, not extant) is a solid black DIAMOND, whether it
    sits at a tip (dated fossil) or on a branch (sampled ancestor).
"""

from __future__ import annotations

import ete3

from fig_species_tree_extinct import annotate_depths, draw_skeleton  # noqa: F401  (re-exported)
from zombi_style import INK, PANEL


def zombi_to_ete3(ztree) -> ete3.Tree:
    """Convert a ZOMBI2 ``Tree`` to an ete3 tree, carrying is_extant / sampled / dist."""

    def rec(zn) -> ete3.TreeNode:
        e = ete3.TreeNode(name=zn.name or "")
        e.dist = 0.0 if zn.parent is None else zn.branch_length()
        e.add_feature("is_extant", bool(zn.is_extant))
        e.add_feature("sampled", bool(getattr(zn, "sampled", False)))
        for c in zn.children:
            e.add_child(rec(c))
        return e

    return rec(ztree.root)


def mark_observed(tree) -> None:
    """Set ``node.observed`` and ``node.has_survivor`` (this node, or something below
    it, is observed) — the latter is what :func:`draw_skeleton` dashes on.

    A node is *itself* an observation when it is a sample: an extant sampled tip, or a
    fossil — whether a dated fossil tip OR a **sampled ancestor** (a fossil on a lineage
    that keeps going). Speciation nodes are not samples. Crucially the node's own sample
    status counts, so the branch *into* a fossil is solid even if everything below it is
    unobserved (fossil-then-extinct): solid up to the diamond, dashed after."""
    for n in tree.traverse("postorder"):
        n.observed = bool(n.sampled or (n.is_leaf() and n.is_extant))
        n.has_survivor = n.observed or any(c.has_survivor for c in n.children)


def fossil_nodes(tree) -> list:
    """Nodes that are fossil samples (sampled, not extant) — tips *and* sampled ancestors."""
    return [n for n in tree.traverse() if n.sampled and not n.is_extant]


def draw_fossils(d, nodes, r: float = 6.5) -> None:
    """Draw a solid black diamond at each fossil node."""
    for n in nodes:
        d._draw_shape_at(*n.coordinates, "square", INK, r=r, stroke=PANEL,
                         stroke_width=1.1, rotation=45.0)
