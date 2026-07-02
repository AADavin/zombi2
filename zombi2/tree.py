"""A lightweight, time-aware tree.

Chosen over ete3/dendropy because the load-bearing quantity in both ZOMBI2 algorithms
is the *absolute node time* (the backward sampler assembles a ranked tree from times;
the forward gene-family loop asks "which branches are alive at time ``t``"). We store
``time`` as a first-class attribute and *derive* branch lengths, avoiding the
``dist = child - parent`` bookkeeping that is a frequent bug source. Interop with the
wider ecosystem is preserved by Newick I/O.

Time convention: ``time`` increases forward from the origin. The root sits at
``time == 0`` and the extant leaves at ``time == total_age``. A *branch* is identified
by its child endpoint and spans the half-open interval ``(parent.time, node.time]``.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(eq=False)
class TreeNode:
    """A node in a :class:`Tree`.

    ``time`` is absolute (forward from origin=0). ``branch_length`` is derived from the
    parent, never stored.
    """

    name: str
    time: float
    parent: "TreeNode | None" = field(default=None, repr=False)
    children: list["TreeNode"] = field(default_factory=list, repr=False)
    is_extant: bool = True

    def add_child(self, child: "TreeNode") -> None:
        child.parent = self
        self.children.append(child)

    def branch_length(self) -> float:
        """Length of the branch above this node (0 for the root)."""
        return 0.0 if self.parent is None else self.time - self.parent.time

    def is_leaf(self) -> bool:
        return not self.children

    def __repr__(self) -> str:  # concise, non-recursive
        kind = "leaf" if self.is_leaf() else "internal"
        return f"TreeNode({self.name!r}, time={self.time:.6g}, {kind})"


class Tree:
    """A rooted tree with absolute node times."""

    def __init__(self, root: TreeNode, total_age: float):
        self.root = root
        self.total_age = total_age

    # --- traversal ---------------------------------------------------------
    def nodes_preorder(self):
        """Yield nodes root-first (pre-order)."""
        stack = [self.root]
        while stack:
            node = stack.pop()
            yield node
            # push children reversed so they come out in order
            stack.extend(reversed(node.children))

    def nodes(self) -> list[TreeNode]:
        return list(self.nodes_preorder())

    def leaves(self) -> list[TreeNode]:
        return [n for n in self.nodes_preorder() if n.is_leaf()]

    def internal_nodes(self) -> list[TreeNode]:
        return [n for n in self.nodes_preorder() if not n.is_leaf()]

    def extant_leaves(self) -> list[TreeNode]:
        return [n for n in self.leaves() if n.is_extant]

    # --- the primitive the forward loop needs ------------------------------
    def branches_alive_at(self, t: float) -> list[TreeNode]:
        """Branches (identified by child endpoint) spanning time ``t``.

        A branch is alive at ``t`` when ``parent.time < t <= node.time``. The root has
        no branch above it and is never returned.
        """
        out = []
        for n in self.nodes_preorder():
            if n.parent is None:
                continue
            if n.parent.time < t <= n.time:
                out.append(n)
        return out

    # --- I/O ---------------------------------------------------------------
    def to_newick(self, include_internal_names: bool = True) -> str:
        """Serialize to Newick, with branch lengths derived from node times."""

        def rec(node: TreeNode) -> str:
            if node.children:
                inner = ",".join(rec(c) for c in node.children)
                label = node.name if include_internal_names else ""
                s = f"({inner}){label}"
            else:
                s = node.name
            if node.parent is not None:
                s += f":{node.branch_length():.10g}"
            return s

        return rec(self.root) + ";"

    def __repr__(self) -> str:
        return f"Tree(root={self.root.name!r}, n_leaves={len(self.leaves())}, total_age={self.total_age:.6g})"


def read_newick(newick: str) -> Tree:
    """Parse a Newick string into a :class:`Tree`.

    Branch lengths are read as **durations**: each node's absolute ``time`` is set to
    ``parent.time + length`` (root at 0), so ``branch_length()`` returns the parsed length.
    Unnamed leaves/internal nodes are given names. Useful for loading an externally
    supplied species tree, or a reconstructed gene tree, into the ``Tree`` interface.
    """
    s = newick.strip().rstrip(";")
    i = 0

    def parse() -> TreeNode:
        nonlocal i
        children = []
        if i < len(s) and s[i] == "(":
            i += 1
            while True:
                children.append(parse())
                if s[i] == ",":
                    i += 1
                elif s[i] == ")":
                    i += 1
                    break
        start = i
        while i < len(s) and s[i] not in ",():;":
            i += 1
        name = s[start:i]
        length = 0.0
        if i < len(s) and s[i] == ":":
            i += 1
            start = i
            while i < len(s) and s[i] not in ",():;":
                i += 1
            length = float(s[start:i])
        node = TreeNode(name=name, time=0.0)
        node._parsed_length = length  # scratch, consumed below
        for c in children:
            node.add_child(c)
        return node

    root = parse()
    tree = Tree(root, 0.0)
    for node in tree.nodes_preorder():
        if node.parent is not None:
            node.time = node.parent.time + node._parsed_length
        del node._parsed_length

    leaf_counter = internal_counter = 1
    for node in tree.nodes_preorder():
        if node.is_leaf() and not node.name:
            node.name = f"n{leaf_counter}"
            leaf_counter += 1
        elif not node.is_leaf() and not node.name:
            node.name = f"i{internal_counter}"
            internal_counter += 1

    tree.total_age = max((leaf.time for leaf in tree.leaves()), default=0.0)
    return tree
