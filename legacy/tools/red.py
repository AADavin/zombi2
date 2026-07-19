"""Relative Evolutionary Divergence (RED) — Parks et al. 2018.

RED places every node of a rooted tree on a ``[0, 1]`` scale of *relative* divergence: the
root is ``0``, every leaf is ``1``, and an internal node sits at the fraction of the
root-to-tip path (measured in branch length) that has elapsed by the time evolution reaches
it. Because it is defined by branch-length *proportions*, RED is invariant to any global rescaling
of the tree — multiplying every branch by a constant leaves RED unchanged — so it reads a
*rate-distorted* phylogram (a tree whose branch lengths are substitutions, not time) as an
approximate, normalised timeline. That is what makes it useful: it recovers relative node ages
without a molecular clock. It is the quantity the **Genome Taxonomy Database (GTDB)** uses to
normalise taxonomic ranks across lineages that evolve at very different rates (Parks et al.,
*Nat. Biotechnol.* 2018; Rinke et al., *Nat. Microbiol.* 2021).

The definition. Assign ``RED(root) = 0``. Visiting nodes from the root outward, for a non-root
node *n* with parent *p*::

    RED(n) = RED(p) + a / (a + b) · (1 − RED(p)),

where ``a`` is the length of the branch above *n* and ``b`` is the mean branch-length distance
from *n* to the leaves of its subtree. A leaf has ``b = 0``, so ``a / (a + b) = 1`` and every
leaf lands exactly at ``1``; the interpolation shares the remaining ``1 − RED(p)`` between a
node and its descendants in proportion to how far down the tree it is.

This tool computes a **number with a right answer** — on an ultrametric tree it returns each
node's true relative age (``node.time / total_age``) exactly — so it sits squarely inside the
:mod:`zombi2.tools` bar. Pair it with a :mod:`zombi2.sequences.clocks` clock to measure how well
RED survives rate variation (see ``docs/tools/red.md``).
"""

from __future__ import annotations

from collections.abc import Callable

from zombi2.tree import Tree, TreeNode

__all__ = ["relative_evolutionary_divergence"]


def relative_evolutionary_divergence(
    tree: "Tree",
    *,
    branch_length: Callable[["TreeNode"], float] | None = None,
) -> dict["TreeNode", float]:
    """Relative Evolutionary Divergence of every node (root ``0.0``, leaves ``1.0``).

    Parameters
    ----------
    tree:
        A :class:`~zombi2.tree.Tree` (its branch lengths are used) **or** a
        :class:`~zombi2.sequences.clocks.RateScaledTree` — the phylogram produced by applying a
        clock — in which case its per-branch substitution lengths are used. A tree read from a
        Newick phylogram via :func:`~zombi2.tree.read_newick` is the common input: its
        ``branch_length()`` is the parsed (substitution) length, so RED reads the phylogram.
    branch_length:
        Optional accessor overriding the branch length used for a node, ``node -> length``
        (default: ``node.branch_length()``). Not allowed together with a ``RateScaledTree``.

    Returns
    -------
    dict
        ``{node: RED}`` for every node in ``tree``, with ``RED(root) == 0.0`` and
        ``RED(leaf) == 1.0`` (except a leaf on a zero-length branch, which inherits its parent's
        value). RED increases monotonically from the root to every tip.

    Notes
    -----
    RED is invariant to a global rescaling of branch lengths, so only the *relative* branch
    proportions matter, not the overall rate. A branch of length zero contributes nothing: its
    child takes the parent's RED. Multifurcations are handled; the tree need not be binary.

    Raises
    ------
    ValueError
        If the tree is empty, or any branch length is negative.
    """
    # A RateScaledTree exposes `.tree` (the topology) and `.branch_lengths` (node -> substitutions).
    if hasattr(tree, "branch_lengths") and hasattr(tree, "tree"):
        if branch_length is not None:
            raise ValueError(
                "pass either a RateScaledTree or an explicit branch_length accessor, not both")
        _bl = tree.branch_lengths
        topology = tree.tree
        length = _bl.__getitem__
    else:
        topology = tree
        length = branch_length if branch_length is not None else TreeNode.branch_length

    nodes = list(topology.nodes_preorder())          # parent before child
    if not nodes:
        raise ValueError("empty tree — nothing to compute RED on")

    # Bottom-up: b(n) = mean branch-length distance from n to the leaves of its subtree.
    mean_tip_dist: dict[TreeNode, float] = {}
    n_leaves: dict[TreeNode, int] = {}
    for node in reversed(nodes):                     # child before parent (valid post-order)
        if node.is_leaf():
            mean_tip_dist[node] = 0.0
            n_leaves[node] = 1
            continue
        total = 0.0
        k = 0
        for child in node.children:
            a = length(child)
            if a < 0.0:
                raise ValueError(f"negative branch length ({a}) above node {child.name!r}")
            total += n_leaves[child] * (a + mean_tip_dist[child])
            k += n_leaves[child]
        mean_tip_dist[node] = total / k
        n_leaves[node] = k

    # Top-down: interpolate RED from the root (0) toward the leaves (1).
    red: dict[TreeNode, float] = {}
    for node in nodes:                               # parent before child
        if node.parent is None:
            red[node] = 0.0
            continue
        a = length(node)
        b = mean_tip_dist[node]
        parent_red = red[node.parent]
        red[node] = parent_red + (a / (a + b)) * (1.0 - parent_red) if (a + b) > 0.0 else parent_red
    return red
