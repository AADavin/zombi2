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
    #: True if this node is an observation — an extant *sampled* tip or a serially-sampled
    #: fossil (see fossilized birth–death in ``simulate_species_tree(..., direction="forward")``).
    #: Default False.
    sampled: bool = False

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


def prune(tree: Tree, keep: str = "extant") -> Tree | None:
    """Prune a complete tree to a reconstructed one, keeping only lineages leading to a kept
    leaf, suppressing degree-two nodes and preserving node times.

    ``keep="extant"`` (default) keeps sampled **extant** tips — the survivors-only reconstructed
    tree (the inverse of forward simulation / ghost-grafting). ``keep="sampled"`` keeps every
    **sample** (serially-sampled fossils *and* sampled extant tips), and preserves *sampled
    ancestors* as degree-two nodes — the fossilized-birth–death sampled tree of dated tips.
    Returns ``None`` if nothing is kept.
    """
    if keep == "extant":
        predicate = lambda n: n.is_extant  # noqa: E731
        keep_ancestors = False
    elif keep == "sampled":
        predicate = lambda n: n.sampled  # noqa: E731
        keep_ancestors = True
    else:
        raise ValueError(f"keep must be 'extant' or 'sampled', got {keep!r}")

    def copy(node: TreeNode) -> TreeNode:
        return TreeNode(name=node.name, time=node.time,
                        is_extant=node.is_extant, sampled=node.sampled)

    def rec(node: TreeNode) -> TreeNode | None:
        if node.is_leaf():
            return copy(node) if predicate(node) else None
        kept = [k for k in (rec(c) for c in node.children) if k is not None]
        if not kept:
            # a kept ancestor with no kept descendants becomes a terminal tip
            return copy(node) if (keep_ancestors and predicate(node)) else None
        if len(kept) == 1:
            if keep_ancestors and predicate(node):  # sampled ancestor -> degree-two node
                sa = copy(node)
                sa.add_child(kept[0])
                return sa
            return kept[0]  # suppress a plain degree-two node
        new = copy(node)
        for k in kept:
            new.add_child(k)
        return new

    root = rec(tree.root)
    return Tree(root, tree.total_age) if root is not None else None


def read_newick(newick: str) -> Tree:
    """Parse a Newick string into a :class:`Tree`.

    Branch lengths are read as **durations**: each node's absolute ``time`` is set to
    ``parent.time + length`` (root at 0), so ``branch_length()`` returns the parsed length.
    Unnamed leaves/internal nodes are given names. Useful for loading an externally
    supplied species tree, or a reconstructed gene tree, into the ``Tree`` interface.
    """
    s = newick.strip().rstrip(";")
    if not s:
        raise ValueError("empty Newick string — is the tree file empty?")
    i = 0

    def skip_ws() -> None:
        # unquoted whitespace (spaces, newlines from a pretty-printed / line-wrapped file) is
        # insignificant in Newick — skip it between tokens so it never leaks into a node name
        nonlocal i
        while i < len(s) and s[i].isspace():
            i += 1

    def read_name() -> str:
        # a quoted label ('...' or "...") is read verbatim — embedded whitespace preserved, a
        # doubled quote ('') unwrapped to one — so e.g. 'Homo sapiens' stays a single name. An
        # unquoted label runs to the next whitespace or structural char, so insignificant
        # whitespace never leaks in.
        nonlocal i
        if i < len(s) and s[i] in "'\"":
            quote = s[i]
            i += 1
            chars: list[str] = []
            while i < len(s):
                if s[i] == quote:
                    if i + 1 < len(s) and s[i + 1] == quote:   # doubled quote = a literal quote
                        chars.append(quote)
                        i += 2
                        continue
                    i += 1
                    break
                chars.append(s[i])
                i += 1
            return "".join(chars)
        start = i
        while i < len(s) and s[i] not in ",():;" and not s[i].isspace():
            i += 1
        return s[start:i]

    def parse() -> TreeNode:
        nonlocal i
        children = []
        skip_ws()
        if i < len(s) and s[i] == "(":
            i += 1
            while True:
                children.append(parse())
                skip_ws()
                if s[i] == ",":
                    i += 1
                elif s[i] == ")":
                    i += 1
                    break
        skip_ws()
        name = read_name()
        length = 0.0
        skip_ws()
        if i < len(s) and s[i] == ":":
            i += 1
            skip_ws()
            start = i
            while i < len(s) and s[i] not in ",():;" and not s[i].isspace():
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
