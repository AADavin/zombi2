"""Trees — the shared dated tree datatype and its toolkit.

``Tree`` is the object every level rides on: the species engine and :func:`read_newick` produce it;
genomes, sequences and traits consume it. It lives here (not in ``zombi2.species``) so the datatype
and everything you do to a tree share one home — one import, ``from zombi2 import tree``.

``Tree`` stays a lean dataclass: its **methods** are only structural self-queries (``leaves``,
``extant``, ``extinct``, ``unsampled``, ``to_newick``). Everything that transforms a tree into a new
tree, or analyses one, is a **free function** in this module (``prune``, ``read_newick``, …), so the
toolkit grows by adding functions, never by growing the class.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass


@dataclass
class Node:
    """One lineage segment: born at ``birth_time``, ended at ``end_time`` by a split, a
    death, or reaching the present. A split has two ``children``; a leaf has none."""

    id: int
    parent: int | None
    birth_time: float
    end_time: float | None = None
    children: tuple[int, int] | None = None
    fate: str = "alive"  # alive → "extant" | "extinct" | "unsampled"; internal splits are "speciation"



@dataclass
class Tree:
    """The complete tree: every lineage that ever lived, keyed by id, rooted at ``root``."""

    nodes: dict[int, Node]
    root: int

    def leaves(self) -> list[Node]:
        """Every lineage with no descendants — extant **and** extinct."""
        return [n for n in self.nodes.values() if n.children is None]

    def extant(self) -> list[Node]:
        """The lineages alive at the present."""
        return [n for n in self.nodes.values() if n.fate == "extant"]

    def extinct(self) -> list[Node]:
        """The lineages that died before the present."""
        return [n for n in self.nodes.values() if n.fate == "extinct"]

    def unsampled(self) -> list[Node]:
        """Survivors not observed under incomplete ``sampling`` — kept in the complete tree (told
        apart by their fate) but pruned from the extant tree."""
        return [n for n in self.nodes.values() if n.fate == "unsampled"]

    def to_newick(self) -> str:
        """Serialise to Newick (matching ``tree.to_newick()`` elsewhere in the codebase). Each
        branch length is ``end_time - birth_time`` and every node — leaves and internals — is named
        ``n<id>``.

        The root carries a branch length like any other node: its **stem**, the time from the origin
        to the first split. A forward birth–death run starts from one lineage, so that stem is real
        simulated time in which events happen, and writing ``)n0;`` would silently discard it — for a
        tree whose crown comes late, a large fraction of its history. It is emitted as ``)n0:<stem>;``
        and :func:`read_newick` reads it back."""

        def emit(i: int) -> str:
            node = self.nodes[i]
            bl = node.end_time - node.birth_time
            if node.children is None:
                return f"n{i}:{bl:.6g}"
            inner = ",".join(emit(c) for c in node.children)
            return f"({inner})n{i}:{bl:.6g}"

        root = self.nodes[self.root]
        stem = root.end_time - root.birth_time
        if root.children is None:
            return f"n{self.root}:{stem:.6g};"
        return f"({','.join(emit(c) for c in root.children)})n{self.root}:{stem:.6g};"


def prune(tree: Tree, keep: str = "extant") -> Tree | None:
    """Prune the complete tree to a kept set (matching ``prune(tree, keep=...)`` in the codebase):
    drop the pruned subtrees and suppress the unifurcations they leave behind, giving a dated,
    bifurcating tree. Branch lengths merge across suppressed nodes; ``None`` if nothing is kept.

    ``keep="extant"`` (default) keeps the survivors — the extant tree. (``"sampled"``, the
    fossil/serially-sampled tree, arrives with the sampling and fossils slices; it raises for now.)"""
    if keep != "extant":
        raise ValueError(f"keep must be 'extant' until the sampling/fossils slices land, got {keep!r}")
    nodes = tree.nodes
    surviving: dict[int, bool] = {}
    for i in sorted(nodes, reverse=True):  # children have higher ids → processed before parents
        nd = nodes[i]
        surviving[i] = nd.fate == "extant" if nd.children is None else any(surviving[c] for c in nd.children)
    if not any(surviving.values()):
        return None

    def surv_children(i: int) -> list[int]:
        nd = nodes[i]
        return [] if nd.children is None else [c for c in nd.children if surviving[c]]

    # keep the extant leaves and the genuine bifurcations (≥2 surviving children)
    kept = {i for i in nodes
            if (nodes[i].children is None and nodes[i].fate == "extant") or len(surv_children(i)) >= 2}

    new: dict[int, Node] = {}
    ext_root: int | None = None
    for i in kept:
        p = nodes[i].parent  # walk up to the nearest kept ancestor
        while p is not None and p not in kept:
            p = nodes[p].parent
        branch_start = nodes[p].end_time if p is not None else 0.0  # merge the suppressed edges
        new[i] = Node(i, p, branch_start, nodes[i].end_time, None, nodes[i].fate)
        if p is None:
            ext_root = i
    for i in sorted(kept):  # rebuild children from parents, in id order for a stable Newick
        p = new[i].parent
        if p is not None:
            existing = new[p].children
            new[p].children = (i,) if existing is None else existing + (i,)

    return Tree(new, ext_root)


_ZOMBI_LABEL = re.compile(r"^n(\d+)$")  # the id-bearing label to_newick writes on every node


def _assign_fates_from_map(leaves: list[Node], labels: dict[int, str],
                           tip_fates: dict[str, str], *, source: str) -> None:
    """Set each leaf's fate from ``tip_fates`` (``{label: fate}``), joined to the leaves through
    ``labels`` (``{leaf id: label}``). Every leaf must be uniquely labelled and covered, and every
    value must be ``extant`` / ``extinct`` / ``unsampled``; anything off raises, naming ``source`` (the
    map it came from) so the message points at the right file. Used for both a ZOMBI tree keyed by its
    ``n<id>`` labels and an external tree keyed by the user's labels."""
    labelled = [labels.get(n.id) for n in leaves]  # every tip must be uniquely named to map a fate
    if any(lbl is None for lbl in labelled):
        raise ValueError(f"every tip must be named to declare its fate with {source}, but the "
                         "tree has unlabelled tips")
    if len(set(labelled)) != len(labelled):
        dups = sorted({lbl for lbl in labelled if labelled.count(lbl) > 1})
        raise ValueError(f"tip labels must be unique to map fates; repeated: {', '.join(dups)}")
    fates = {k: str(v).lower() for k, v in tip_fates.items()}
    missing = sorted(set(labelled) - set(fates))
    unknown = sorted(set(fates) - set(labelled))
    if missing:
        raise ValueError(f"{source} is missing a fate for: {', '.join(missing)}")
    if unknown:
        raise ValueError(f"{source} names tips not in the tree: {', '.join(unknown)}")
    bad = sorted({v for v in fates.values() if v not in ("extant", "extinct", "unsampled")})
    if bad:
        raise ValueError(f"tip fates must be 'extant', 'extinct' or 'unsampled', got: {', '.join(bad)}")
    for n in leaves:
        n.fate = fates[labels[n.id]]


def _assign_external_fates(leaves: list[Node], names: dict[int, str],
                           tip_fates: dict[str, str] | None, gap: float) -> None:
    """Set the fate on the tips of a **non-ultrametric** external tree from a user-supplied
    ``tip_fates`` (``{tip label: "extant" | "extinct" | "unsampled"}`` — the same vocabulary a species
    run writes to ``species_fates.tsv``). ZOMBI will not infer a tip's fate from its depth here (a
    shallow tip could be extinct *or* an early sample), so a missing or mismatched map raises rather
    than guesses."""
    if tip_fates is None:
        raise ValueError(
            f"input tree is not ultrametric (tip depths differ by {gap:.3g}); ZOMBI can't tell "
            "extinct lineages from early samples — declare each tip's fate with --tip-fates FILE "
            "(a 'tip_name<TAB>extant|extinct|unsampled' row per tip)")
    _assign_fates_from_map(leaves, names, tip_fates, source="--tip-fates")


def read_newick(newick: str, *, tip_fates: dict[str, str] | None = None) -> tuple[Tree, dict[int, str]]:
    """Parse a Newick string into a complete :class:`Tree` and a name-map ``{id: user label}``.

    This is how the CLI loads a species tree back for the downstream levels. Branch lengths are read
    as **durations**: the root sits at time 0 and each node's ``birth_time`` is its parent's
    ``end_time``, so ``end_time - birth_time`` is the parsed length. Two kinds of tree are accepted,
    told apart by the labels:

    - a **ZOMBI complete tree** (every node — internal ones too — is ``n<id>``, as ``to_newick``
      writes it): the ids come from the labels, and the name-map is empty (the labels *are* the ids).
      Fate comes from ``tip_fates`` when given — the run's ``species_fates.tsv``, keyed by the same
      ``n<id>`` — which is authoritative; without it a leaf is ``"extinct"`` if it ends before the
      tree's greatest depth, else ``"extant"`` (a fallback that cannot recover an ``"unsampled"`` tip).
    - any **external tree** (leaves named freely, internal nodes usually unlabelled): fresh ids are
      minted in traversal order (root 0, parents before children), the original labels are returned as
      the **name-map** (``{minted id: user label}``), and fates depend on whether the tree is
      **ultrametric** (all root-to-tip depths equal, within ``1e-6 × height``):

      - **ultrametric** → the tips are contemporaneous, so **every tip is ``"extant"``** (observed);
      - **not ultrametric** → the differing tip depths could mean extinct lineages *or* early
        samples, which ZOMBI cannot tell apart, so it **refuses to guess**: pass ``tip_fates`` — a
        ``{tip label: "extant" | "extinct" | "unsampled"}`` map covering every tip — or a
        :class:`ValueError` is raised. (The CLI fills ``tip_fates`` from ``--tip-fates FILE``, which
        reads the same format a species run writes to ``species_fates.tsv``.)

    A root branch length is read when present — ``to_newick`` writes one, so a ZOMBI tree round-trips
    with its stem intact. External trees usually have none, and then the root gets zero duration and
    the tree starts at its crown, which is all the file says.

    The ``.nwk`` records only branch lengths, so ``"unsampled"`` fate cannot be read from it alone —
    pass ``tip_fates`` (the run's ``species_fates.tsv``) to recover it; without one a survivor reads
    back ``"extant"``, which is still fine for evolving genomes/traits along the tree.

    Only bifurcating trees are supported (an internal node with other than two children raises).
    """
    s = newick.strip().rstrip(";").strip()
    if not s:
        raise ValueError("empty Newick string — is the tree file empty?")
    i = 0

    def skip_ws() -> None:
        # whitespace (incl. the newlines of a line-wrapped file) is insignificant between tokens
        nonlocal i
        while i < len(s) and s[i].isspace():
            i += 1

    def read_name() -> str:
        # a quoted label ('...' / "...") is taken verbatim (a doubled quote unwrapped to one); an
        # unquoted label runs to the next whitespace or structural char, so stray whitespace never leaks
        nonlocal i
        if i < len(s) and s[i] in "'\"":
            quote = s[i]
            i += 1
            chars: list[str] = []
            while i < len(s):
                if s[i] == quote:
                    if i + 1 < len(s) and s[i + 1] == quote:
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

    # parse into a lightweight tree of (name, length, children) so we can decide ids/fates in a
    # second pass, once we know whether every node is ``n<id>``-labelled and the tree's depth
    class _P:
        __slots__ = ("name", "length", "children")

        def __init__(self, name, length, children):
            self.name, self.length, self.children = name, length, children

    def parse() -> _P:
        nonlocal i
        children: list[_P] = []
        skip_ws()
        if i < len(s) and s[i] == "(":
            i += 1
            while True:
                children.append(parse())
                skip_ws()
                if i >= len(s):
                    raise ValueError("malformed Newick: unbalanced parentheses — the string ended "
                                     "while a clade was still open (is the tree file truncated?)")
                if s[i] == ",":
                    i += 1
                elif s[i] == ")":
                    i += 1
                    break
                else:
                    raise ValueError(f"malformed Newick: expected ',' or ')' after a clade at "
                                     f"position {i}, got {s[i]!r}")
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
            try:
                length = float(s[start:i])
            except ValueError:
                raise ValueError(f"malformed Newick: expected a branch length after ':' at position "
                                 f"{start}, got {s[start:i]!r}") from None
        if children and len(children) != 2:
            raise ValueError(f"only bifurcating trees are supported: node {name or '(unnamed)'!r} has "
                             f"{len(children)} children (collapse polytomies / unifurcations first)")
        return _P(name, length, children)

    root_p = parse()
    skip_ws()
    if i < len(s):
        raise ValueError(f"malformed Newick: unexpected trailing text at position {i}: {s[i:]!r}")

    # ZOMBI complete tree ⟺ every node carries an ``n<id>`` label; then ids come from the labels,
    # otherwise we mint them ourselves and call every tip extant.
    all_labelled = True

    def _scan(p: _P) -> None:
        nonlocal all_labelled
        if not _ZOMBI_LABEL.match(p.name):
            all_labelled = False
        for c in p.children:
            _scan(c)

    _scan(root_p)

    nodes: dict[int, Node] = {}
    names: dict[int, str] = {}  # {minted id: user label} — for external trees; empty for ZOMBI ones
    counter = 0

    def _mint(p: _P) -> int:
        nonlocal counter
        if all_labelled:
            return int(_ZOMBI_LABEL.match(p.name).group(1))
        i_ = counter
        counter += 1
        return i_

    # first pass: assign ids (parents before children) and absolute times from durations
    def _build(p: _P, parent: int | None, birth: float) -> int:
        nid = _mint(p)
        if nid in nodes:
            raise ValueError(f"duplicate node id n{nid} in the Newick (labels must be unique)")
        end = birth + p.length
        child_ids = tuple(_build(c, nid, end) for c in p.children) or None
        nodes[nid] = Node(nid, parent, birth, end, child_ids)  # fate filled in below
        if not all_labelled and p.name:  # external tree: keep the user's label for the name-map
            names[nid] = p.name
        return nid

    root_id = _build(root_p, None, 0.0)

    # second pass: fates. Internal nodes are always speciations; the tips depend on the tree kind.
    for n in nodes.values():
        if n.children is not None:
            n.fate = "speciation"
    leaves = [n for n in nodes.values() if n.children is None]

    if all_labelled:
        if tip_fates is not None:
            # the run's own species_fates.tsv (or a --tip-fates file) states each tip's fate directly,
            # keyed by the same n<id> the tree carries. It is authoritative: depth cannot tell an
            # unsampled survivor from an extant one (both sit at the present), nor an extinct tip that
            # died just before the present, so when it is given we use it rather than guess.
            _assign_fates_from_map(leaves, {n.id: f"n{n.id}" for n in leaves}, tip_fates,
                                   source="tip fates")
            return Tree(nodes, root_id), names
        # no fate table: a tip is extinct if it ends before the present (the greatest end_time). The
        # tolerance is depth-relative — ``to_newick`` prints 6 significant figures, whose rounding
        # accumulates to ~5e-6·height along a root-to-tip path, so a tip at the present can fall a few
        # ×1e-5·height short of the max, far below any real extinction gap. This cannot recover an
        # unsampled tip (it sits at the present, so it reads back extant) — pass the fate table for that.
        present = max(n.end_time for n in nodes.values())
        tol = max(1e-9, 1e-4 * present)
        for n in leaves:
            n.fate = "extinct" if n.end_time < present - tol else "extant"
        return Tree(nodes, root_id), names

    # an external tree: ultrametric ⟺ every tip is contemporaneous (all extant); otherwise the
    # differing depths could be extinctions or early samples, which we refuse to guess (SPEC decision).
    depths = [n.end_time for n in leaves]  # root sits at 0, so a tip's depth is its end_time
    height = max(depths)
    if max(depths) - min(depths) <= max(1e-12, 1e-6 * height):  # ultrametric
        for n in leaves:
            n.fate = "extant"
        return Tree(nodes, root_id), names

    _assign_external_fates(leaves, names, tip_fates, max(depths) - min(depths))
    return Tree(nodes, root_id), names




# --------------------------------------------------------------------------------------------------
# The toolkit — free functions (Tree → Tree transforms, and analyses). Grow it here, not on the class.
# --------------------------------------------------------------------------------------------------


def _copy(tree: Tree) -> Tree:
    return Tree({i: Node(n.id, n.parent, n.birth_time, n.end_time, n.children, n.fate)
                 for i, n in tree.nodes.items()}, tree.root)


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


def _depths(tree: Tree) -> dict[int, float]:
    """Root-to-node branch-length depth of every node (root at 0)."""
    depth: dict[int, float] = {}
    for i in _preorder(tree):
        nd = tree.nodes[i]
        depth[i] = 0.0 if nd.parent is None else depth[nd.parent] + (nd.end_time - nd.birth_time)
    return depth


def with_stem(tree: Tree, length: float, *, mode: str = "set") -> Tree:
    """Return a copy whose **stem** — the branch above the crown (root) — is ``length`` (``mode="set"``)
    or is extended by ``length`` (``mode="add"``). Every other branch length is unchanged, so
    ``to_newick`` writes the new stem as ``)n<root>:<stem>;`` and nothing below moves."""
    if not math.isfinite(length):
        raise ValueError(f"stem length must be finite, got {length!r}")
    if mode not in ("set", "add"):
        raise ValueError(f"mode must be 'set' or 'add', got {mode!r}")
    out = _copy(tree)
    root = out.nodes[out.root]
    root.birth_time = (root.end_time - length) if mode == "set" else (root.birth_time - length)
    if root.end_time - root.birth_time < 0:
        raise ValueError("resulting stem is negative")
    return out


def make_ultrametric(tree: Tree, *, tol: float = 1e-3) -> Tree:
    """Return a copy in which every tip sits at the present (exactly ultrametric), by extending the
    terminal branches to a common depth. Snaps only when the tip-depth spread is within ``tol`` of
    the tree height — i.e. rounding; a larger spread raises, because differing tip depths then carry
    real signal (extinct lineages or serial samples) that this must not silently flatten."""
    depth = _depths(tree)
    tips = [i for i, n in tree.nodes.items() if n.children is None]
    lo, hi = min(depth[i] for i in tips), max(depth[i] for i in tips)
    if hi > 0 and (hi - lo) > tol * hi:
        raise ValueError(
            f"tip depths differ by {hi - lo:.3g} (> {tol:g} × height {hi:.3g}); this is more than "
            "rounding — the tips are not contemporaneous (extinct lineages or serial samples), so "
            "there is no ultrametric tree to snap to")
    out = _copy(tree)
    for i in tips:
        nd = out.nodes[i]
        parent_depth = depth[i] - (nd.end_time - nd.birth_time)   # = depth of this tip's parent
        nd.end_time = nd.birth_time + (hi - parent_depth)         # so the tip lands at depth hi
    return out


def rescale(tree: Tree, *, height: float | None = None, factor: float | None = None) -> Tree:
    """Return a copy with every branch length scaled — either so the root-to-tip height equals
    ``height``, or by a raw ``factor``. Exactly one of the two must be given."""
    if (height is None) == (factor is None):
        raise ValueError("pass exactly one of height= or factor=")
    if factor is None:
        depth = _depths(tree)
        current = max(depth[i] for i, n in tree.nodes.items() if n.children is None)
        if current <= 0:
            raise ValueError("tree has zero height; cannot scale it to a target height")
        factor = height / current
    if factor < 0:
        raise ValueError(f"scale factor must be non-negative, got {factor}")
    out = _copy(tree)
    for nd in out.nodes.values():
        nd.birth_time *= factor
        nd.end_time *= factor
    return out


def relative_evolutionary_divergence(tree: Tree) -> dict[int, float]:
    """Relative Evolutionary Divergence (Parks et al. 2018) of every node — root ``0.0``, leaves
    ``1.0``, keyed by node id. Walking root-outward, a node sits at ``RED(parent) + a/(a+b)·(1 −
    RED(parent))`` where ``a`` is its branch and ``b`` the mean branch-length distance from it to the
    leaves of its subtree. RED is invariant to a global rescaling, so a rate-distorted phylogram reads
    as an approximate relative timeline; on an ultrametric tree it returns each node's exact relative
    age. A zero-length branch passes the parent's value straight down."""
    nodes = tree.nodes
    order = _preorder(tree)
    if not order:
        raise ValueError("empty tree — nothing to compute RED on")

    def length(i: int) -> float:
        a = nodes[i].end_time - nodes[i].birth_time
        if a < 0.0:
            raise ValueError(f"negative branch length ({a}) above node n{i}")
        return a

    mean_tip_dist: dict[int, float] = {}
    n_leaves: dict[int, int] = {}
    for i in reversed(order):                       # child before parent
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

    red: dict[int, float] = {}
    for i in order:                                 # parent before child
        p = nodes[i].parent
        if p is None:
            red[i] = 0.0
            continue
        a, b, pr = length(i), mean_tip_dist[i], red[p]
        red[i] = pr + (a / (a + b)) * (1.0 - pr) if (a + b) > 0.0 else pr
    return red


def red_scaled(tree: Tree) -> Tree:
    """Return a copy whose node depths **are** their RED — ultrametric on ``[0, 1]``, root at 0, every
    tip at 1. Branch lengths become RED increments. This is the tree GTDB-style rank normalisation
    reads (:func:`relative_evolutionary_divergence` gives the raw per-node values)."""
    red = relative_evolutionary_divergence(tree)
    out = _copy(tree)
    for i, nd in out.nodes.items():
        nd.birth_time = 0.0 if nd.parent is None else red[nd.parent]
        nd.end_time = red[i]
    return out


def _clades(tree: Tree) -> dict[frozenset, float]:
    """``{frozenset(descendant leaf ids): branch length above the node}`` for every node."""
    leafset: dict[int, frozenset] = {}
    for i in reversed(_preorder(tree)):             # child before parent
        nd = tree.nodes[i]
        leafset[i] = (frozenset((i,)) if nd.children is None
                      else frozenset().union(*(leafset[c] for c in nd.children)))
    return {leafset[i]: (n.end_time - n.birth_time) for i, n in tree.nodes.items()}


def distance(a: Tree, b: Tree, *, metric: str = "rf") -> float:
    """Distance between two **rooted** trees over their shared tips (matched by node id). Raises if the
    two leaf sets differ. ``metric``: ``"rf"`` (Robinson–Foulds — the number of clades in one tree but
    not the other), ``"rf-normalized"`` (that count over the total number of non-trivial clades), or
    ``"branch-score"`` (Kuhner–Felsenstein — √Σ(branch-length difference)² over all clades, terminal
    branches included)."""
    la = frozenset(i for i, n in a.nodes.items() if n.children is None)
    lb = frozenset(i for i, n in b.nodes.items() if n.children is None)
    if la != lb:
        raise ValueError(
            f"the two trees have different leaf sets ({len(la)} vs {len(lb)} tips, {len(la ^ lb)} not "
            "shared) — treedist needs the same taxa, identically labelled, on both trees")
    ca, cb = _clades(a), _clades(b)
    n = len(la)
    if metric in ("rf", "rf-normalized"):
        sa = {c for c in ca if 2 <= len(c) <= n - 1}
        sb = {c for c in cb if 2 <= len(c) <= n - 1}
        rf = len(sa ^ sb)
        if metric == "rf":
            return float(rf)
        denom = len(sa) + len(sb)
        return float(rf / denom) if denom else 0.0
    if metric == "branch-score":
        return float(sum((ca.get(c, 0.0) - cb.get(c, 0.0)) ** 2 for c in set(ca) | set(cb)) ** 0.5)
    raise ValueError(f"unknown metric {metric!r}; choose 'rf', 'rf-normalized', or 'branch-score'")


__all__ = ["Tree", "Node", "prune", "read_newick", "with_stem", "make_ultrametric", "rescale",
           "relative_evolutionary_divergence", "red_scaled", "distance"]
