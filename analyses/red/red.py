"""Relative Evolutionary Divergence (RED) — Parks et al. 2018 — on the ZOMBI2 clean-core tree.

RED places every node of a rooted tree on a ``[0, 1]`` scale of *relative* divergence: the root
is ``0``, every leaf is ``1``, and an internal node sits at the fraction of the root-to-tip path
(measured in branch length) that has elapsed by the time evolution reaches it. Because it is
defined by branch-length *proportions*, RED is invariant to any global rescaling of the tree, so
it reads a *rate-distorted* phylogram (branch lengths in substitutions, not time) as an
approximate, normalised timeline — a relative molecular clock with no clock assumed.

Definition. ``RED(root) = 0``. Visiting nodes root-outward, for a non-root node *n* with parent
*p*::

    RED(n) = RED(p) + a / (a + b) · (1 − RED(p))

where ``a`` is the branch above *n* and ``b`` is the mean branch-length distance from *n* to the
leaves of its subtree. A leaf has ``b = 0`` so it lands exactly at ``1``.

On an ultrametric tree (branch lengths proportional to time) RED returns each node's true relative
age exactly — so RED of the *dated* tree is the ground truth, and RED of the *ragged phylogram* is
the estimate we grade against it.

This is a faithful port of the retired ``zombi2.tools.relative_evolutionary_divergence`` to the
clean-core integer-keyed :class:`zombi2.species.Tree` (``Node(id, parent, birth_time, end_time,
children, fate)``; branch length ``end_time − birth_time``).
"""
from __future__ import annotations

from zombi2.species import Tree


def _preorder(tree: Tree) -> list[int]:
    """Node ids, parent before child, from the root."""
    order: list[int] = []
    stack = [tree.root]
    while stack:
        i = stack.pop()
        order.append(i)
        kids = tree.nodes[i].children
        if kids is not None:
            stack.extend(kids)
    return order


def relative_evolutionary_divergence(tree: Tree) -> dict[int, float]:
    """RED of every node (root ``0.0``, leaves ``1.0``), keyed by node id.

    Branch length above node *n* is ``end_time − birth_time``. RED is invariant to a global
    rescaling, so only relative branch proportions matter. A zero-length branch passes the parent's
    value straight down. Multifurcations are handled.
    """
    nodes = tree.nodes
    order = _preorder(tree)
    if not order:
        raise ValueError("empty tree — nothing to compute RED on")

    def length(i: int) -> float:
        nd = nodes[i]
        a = nd.end_time - nd.birth_time
        if a < 0.0:
            raise ValueError(f"negative branch length ({a}) above node n{i}")
        return a

    # Bottom-up: b(n) = mean branch-length distance from n to the leaves of its subtree.
    mean_tip_dist: dict[int, float] = {}
    n_leaves: dict[int, int] = {}
    for i in reversed(order):                       # child before parent (valid post-order)
        kids = nodes[i].children
        if kids is None:
            mean_tip_dist[i] = 0.0
            n_leaves[i] = 1
            continue
        total = 0.0
        k = 0
        for c in kids:
            total += n_leaves[c] * (length(c) + mean_tip_dist[c])
            k += n_leaves[c]
        mean_tip_dist[i] = total / k
        n_leaves[i] = k

    # Top-down: interpolate RED from the root (0) toward the leaves (1).
    red: dict[int, float] = {}
    for i in order:                                 # parent before child
        p = nodes[i].parent
        if p is None:
            red[i] = 0.0
            continue
        a = length(i)
        b = mean_tip_dist[i]
        pr = red[p]
        red[i] = pr + (a / (a + b)) * (1.0 - pr) if (a + b) > 0.0 else pr
    return red


def internal_nodes(tree: Tree) -> list[int]:
    """Ids of the nodes RED is graded on: internal, non-root (leaves are trivially 1, root 0)."""
    return [i for i, nd in tree.nodes.items()
            if nd.children is not None and nd.parent is not None]
