"""A compact, likelihood-ready view of a dated species tree.

ALElite keeps its own indexed species-tree structure rather than reaching into
:class:`zombi2.tree.Tree` directly, so the likelihood core depends only on plain integer
arrays (parent/child indices, times) and stays trivially portable — this is the single
module that knows about ZOMBI2's ``Tree``.

Branches are stored in **post-order** (every child index is smaller than its parent's), so a
single forward pass over ``range(n)`` visits children before parents — exactly the order the
extinction and DP recursions need.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class Branch:
    """One species-tree branch, identified by its lower (child) endpoint.

    ``left``/``right`` are the post-order indices of the two daughter branches (``None`` for a
    leaf). ``time`` is the absolute time of the branch's lower endpoint (forward from origin=0;
    leaves sit at ``total_age``); ``parent_time`` is its upper endpoint. ``length = time -
    parent_time``.
    """

    name: str
    is_leaf: bool
    parent: int | None
    left: int | None
    right: int | None
    time: float
    parent_time: float

    @property
    def length(self) -> float:
        return self.time - self.parent_time


class SpeciesTree:
    """A rooted, dated species tree indexed for likelihood evaluation.

    ``branches`` is the post-order list; ``root`` is the index of the root branch;
    ``leaf_index`` maps a leaf's species name to its branch index. Only strictly bifurcating
    trees are supported (degree-two nodes — e.g. fossilized-birth–death sampled ancestors —
    would need suppressing first); a clear error is raised otherwise.
    """

    def __init__(self, branches: list[Branch], root: int, leaf_index: dict[str, int]):
        self.branches = branches
        self.root = root
        self.leaf_index = leaf_index
        self.n = len(branches)

    @classmethod
    def from_tree(cls, tree) -> "SpeciesTree":
        """Build from a :class:`zombi2.tree.Tree` (or anything with the same node API)."""
        branches: list[Branch] = []

        def visit(node) -> int:
            child_idx = [visit(c) for c in node.children]
            if child_idx and len(child_idx) != 2:
                raise ValueError(
                    f"ALElite needs a strictly bifurcating species tree; node {node.name!r} "
                    f"has {len(child_idx)} children (suppress degree-two nodes first)."
                )
            i = len(branches)
            parent_time = 0.0 if node.parent is None else node.parent.time
            branches.append(Branch(
                name=node.name,
                is_leaf=not node.children,
                parent=None,  # filled in below (parent visited after its children)
                left=child_idx[0] if child_idx else None,
                right=child_idx[1] if child_idx else None,
                time=node.time,
                parent_time=parent_time,
            ))
            return i

        root = visit(tree.root)
        _set_parents(branches)  # wire up parent indices from the child links
        leaf_index = {b.name: i for i, b in enumerate(branches) if b.is_leaf}
        return cls(branches, root, leaf_index)

    def __repr__(self) -> str:
        return f"SpeciesTree(n={self.n}, leaves={len(self.leaf_index)}, root={self.root})"


def _set_parents(branches: list[Branch]) -> None:
    """Fill each branch's ``parent`` index from the child links (post-order guarantees the
    parent is encountered after both its children)."""
    for i, b in enumerate(branches):
        if b.left is not None:
            branches[b.left].parent = i
            branches[b.right].parent = i
