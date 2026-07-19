"""Topological and branch-length distances between two trees.

The companion to a simulation: ZOMBI2 emits a *true* tree (species tree, or a reconciled gene
tree), an inference method emits its *estimate*, and this tool reports **how far apart they
are** — one number with a right answer, the standard currency for benchmarking phylogenetic
methods against a known ground truth. Every metric here is exact and defined over the shared
leaf-label set of the two trees; a mismatch in the label sets is an error, not a silent
partial score.

Four classical distances:

* **Robinson–Foulds** (Robinson & Foulds, *Math. Biosci.* 1981) — the symmetric difference of
  the two trees' clade (rooted) or bipartition (unrooted) sets: how many groupings appear in
  one tree but not the other. Fast, exact, and the most widely reported topology distance,
  but coarse — a single misplaced leaf can saturate it.
* **Branch score** (Kuhner & Felsenstein, *Mol. Biol. Evol.* 1994) — like RF but weighted by
  branch *length*: the L2 (or L1) norm of the per-clade branch-length differences, matched
  clades contributing ``(ℓ₁ − ℓ₂)²`` and unmatched clades contributing their full length
  against zero. Uses ZOMBI2's time-derived branch lengths, so it reads a dated tree directly.
* **Quartet distance** (Estabrook, McMorris & Meacham, *Syst. Zool.* 1985) — the fraction of
  four-leaf subsets whose induced (unrooted) topology differs between the trees. Finer-grained
  than RF (it degrades gracefully), but the exact count here is ``O(n⁴)``, so it is guarded by
  ``max_leaves`` and meant for the reconstructed trees people actually score, not million-tip
  simulations.
* **Matching distance** (Bogdanowicz & Giaro, *IEEE/ACM TCBB* 2012; Lin, Rajan & Moret 2012) —
  the minimum-cost *matching* of one tree's splits (unrooted) or clusters (rooted) to the
  other's, each pair costing the number of leaves that must move to make them equal and
  unmatched splits paying their weight against the empty split. Unlike RF it does not saturate:
  a slightly displaced clade costs a little, not a whole unit. Solved as the assignment problem
  via SciPy (an *optional* dependency, lazily imported — the rest of the module needs only
  NumPy); if SciPy is absent, this metric is skipped, not fabricated.

Rooted vs unrooted. ZOMBI2 trees are rooted (``node.time`` runs forward from a root at 0), and
the rooted RF over *clades* is the default — it is the honest distance when the root is known,
as it is for a simulated tree. Pass ``rooted=False`` for the unrooted bipartition RF when
comparing against a method that does not infer a root. The quartet distance is intrinsically
unrooted (a quartet has no root).

Usage::

    from zombi2.tools.treedist import compare_trees
    from zombi2.tree import read_newick

    truth = read_newick(open("species_tree.nwk").read())
    est   = read_newick(open("inferred.nwk").read())
    print(compare_trees(truth, est))          # -> TreeComparison(rf=..., branch_score=..., quartet=...)
"""

from __future__ import annotations

import math
from collections import namedtuple
from itertools import combinations

from zombi2.tree import Tree, read_newick

__all__ = [
    "robinson_foulds",
    "branch_score",
    "quartet_distance",
    "matching_distance",
    "compare_trees",
    "RFResult",
    "QuartetResult",
    "MatchingResult",
    "TreeComparison",
]

#: Result of :func:`robinson_foulds`. ``rf`` is the raw symmetric-difference count, ``max_rf``
#: the largest value it could take (``|A| + |B|`` over the two clade/bipartition sets), and
#: ``normalized`` = ``rf / max_rf`` in ``[0, 1]`` (0 when both trees have no informative splits).
RFResult = namedtuple("RFResult", ["rf", "max_rf", "normalized"])

#: Result of :func:`quartet_distance`. ``differing`` quartets out of ``total`` = ``C(n, 4)``;
#: ``normalized`` = ``differing / total`` in ``[0, 1]``.
QuartetResult = namedtuple("QuartetResult", ["differing", "total", "normalized"])

#: Result of :func:`matching_distance`. ``distance`` is the min-cost matching (an integer number
#: of leaf-moves), ``max_distance`` the cost of matching everything to the empty split (an upper
#: bound), and ``normalized`` = ``distance / max_distance`` in ``[0, 1]``.
MatchingResult = namedtuple("MatchingResult", ["distance", "max_distance", "normalized"])

#: Aggregate of every metric between two trees. ``quartet`` fields are ``None`` when the quartet
#: distance was skipped (``quartet=False`` or more than ``max_leaves`` tips); ``matching`` fields
#: are ``None`` when it was skipped (``matching=False``, too many tips, or SciPy is not installed).
TreeComparison = namedtuple(
    "TreeComparison",
    ["n_leaves", "rf", "rf_max", "rf_normalized", "rf_unrooted", "rf_unrooted_max",
     "branch_score", "quartet", "quartet_normalized", "matching", "matching_normalized"],
)


# --------------------------------------------------------------------------- helpers

def _as_tree(t) -> Tree:
    """Accept a :class:`~zombi2.tree.Tree` or a Newick string."""
    if isinstance(t, Tree):
        return t
    if isinstance(t, str):
        return read_newick(t)
    raise TypeError(f"expected a Tree or a Newick string, got {type(t).__name__}")


def _leaf_names(tree: Tree) -> list[str]:
    return [leaf.name for leaf in tree.leaves()]


def _check_same_leaves(t1: Tree, t2: Tree) -> frozenset[str]:
    """Return the shared leaf-name set, or raise if the two trees disagree on their leaves."""
    a, b = set(_leaf_names(t1)), set(_leaf_names(t2))
    if len(a) != len(t1.leaves()):
        raise ValueError("tree 1 has duplicate leaf names; distances need unique labels")
    if len(b) != len(t2.leaves()):
        raise ValueError("tree 2 has duplicate leaf names; distances need unique labels")
    if a != b:
        only1 = sorted(a - b)[:5]
        only2 = sorted(b - a)[:5]
        raise ValueError(
            "the two trees must have the same leaf set; "
            f"{len(a - b)} only in tree 1 (e.g. {only1}), "
            f"{len(b - a)} only in tree 2 (e.g. {only2})"
        )
    return frozenset(a)


def _node_clades(tree: Tree) -> dict:
    """``{node: frozenset(leaf names in its subtree)}`` for every node (post-order fill)."""
    clades: dict = {}
    for node in reversed(tree.nodes()):          # preorder reversed ⇒ children before parents
        if node.is_leaf():
            clades[node] = frozenset((node.name,))
        else:
            acc: set = set()
            for c in node.children:
                acc |= clades[c]
            clades[node] = frozenset(acc)
    return clades


def _clusters(tree: Tree, n: int) -> set:
    """Rooted informative clades: leaf-sets of size ``2 .. n-1`` (drop singletons and the root)."""
    return {cl for cl in _node_clades(tree).values() if 1 < len(cl) < n}


def _bipartitions(tree: Tree, all_leaves: frozenset) -> set:
    """Unrooted informative splits, canonicalised so ``A|B`` and ``B|A`` collapse to one key."""
    ref = min(all_leaves)
    n = len(all_leaves)
    out: set = set()
    for cl in _node_clades(tree).values():
        if min(len(cl), n - len(cl)) < 2:        # trivial split (a leaf vs the rest)
            continue
        side = cl if ref not in cl else (all_leaves - cl)
        out.add(frozenset(side))
    return out


def _clade_lengths(tree: Tree) -> dict:
    """``{clade: summed branch length}`` over every non-root node (terminal branches included)."""
    clades = _node_clades(tree)
    out: dict = {}
    for node, cl in clades.items():
        if node.parent is None:
            continue
        out[cl] = out.get(cl, 0.0) + node.branch_length()
    return out


# --------------------------------------------------------------------------- metrics

def robinson_foulds(t1, t2, *, rooted: bool = True) -> RFResult:
    """Robinson–Foulds distance between two trees over a shared leaf set.

    ``rooted=True`` (default) compares **clade** sets (the honest distance when the root is
    known, as for a simulated tree); ``rooted=False`` compares unrooted **bipartitions**. In
    both cases the distance is the number of groupings present in exactly one tree, and
    ``normalized = rf / (|splits₁| + |splits₂|)`` lands in ``[0, 1]``.
    """
    t1, t2 = _as_tree(t1), _as_tree(t2)
    leaves = _check_same_leaves(t1, t2)
    n = len(leaves)
    if rooted:
        s1, s2 = _clusters(t1, n), _clusters(t2, n)
    else:
        s1, s2 = _bipartitions(t1, leaves), _bipartitions(t2, leaves)
    rf = len(s1 ^ s2)
    max_rf = len(s1) + len(s2)
    return RFResult(rf, max_rf, (rf / max_rf) if max_rf else 0.0)


def branch_score(t1, t2, *, order: int = 2) -> float:
    """Kuhner–Felsenstein branch-score distance using the trees' branch lengths.

    Clades are matched by their leaf-set; a clade in both trees contributes ``|ℓ₁ − ℓ₂|`` (its
    branch-length difference), a clade in only one contributes its length against zero. With
    ``order=2`` (default) the contributions are squared and the result is their square root
    (the L2 / Euclidean branch score); ``order=1`` returns the L1 sum. Terminal branches are
    included. Rooted (clade-keyed): ZOMBI2 trees carry a root, so this is the natural form.
    """
    if order not in (1, 2):
        raise ValueError(f"order must be 1 or 2, got {order!r}")
    t1, t2 = _as_tree(t1), _as_tree(t2)
    _check_same_leaves(t1, t2)
    b1, b2 = _clade_lengths(t1), _clade_lengths(t2)
    keys = set(b1) | set(b2)
    if order == 2:
        return math.sqrt(sum((b1.get(k, 0.0) - b2.get(k, 0.0)) ** 2 for k in keys))
    return sum(abs(b1.get(k, 0.0) - b2.get(k, 0.0)) for k in keys)


def _leaf_root_paths(tree: Tree) -> dict:
    """``{leaf name: (id(root), ..., id(leaf))}`` — the node-id path from the root to each leaf.

    Two leaves' LCA depth is ``len(common prefix) - 1`` (the root sits at depth 0), which is all
    the quartet topology needs.
    """
    paths: dict = {}

    def rec(node, prefix):
        prefix = prefix + (id(node),)
        if node.is_leaf():
            paths[node.name] = prefix
        else:
            for c in node.children:
                rec(c, prefix)

    rec(tree.root, ())
    return paths


def _lca_depth(px: tuple, py: tuple) -> int:
    k = 0
    for a, b in zip(px, py):
        if a != b:
            break
        k += 1
    return k - 1


def _quartet_topology(a, b, c, d, P) -> frozenset | None:
    """Induced unrooted topology of quartet ``{a,b,c,d}``: the canonical split, or ``None`` if
    unresolved. The grouped pair is the one whose LCA is *deepest*; a tie between two
    non-complementary pairs (or an all-way tie, i.e. a star) is unresolved."""
    d_ab, d_cd = _lca_depth(P[a], P[b]), _lca_depth(P[c], P[d])
    d_ac, d_bd = _lca_depth(P[a], P[c]), _lca_depth(P[b], P[d])
    d_ad, d_bc = _lca_depth(P[a], P[d]), _lca_depth(P[b], P[c])
    m = max(d_ab, d_cd, d_ac, d_bd, d_ad, d_bc)
    winners = set()
    if d_ab == m or d_cd == m:
        winners.add(0)
    if d_ac == m or d_bd == m:
        winners.add(1)
    if d_ad == m or d_bc == m:
        winners.add(2)
    if winners not in ({0}, {1}, {2}):
        return None                              # tie across splits ⇒ unresolved
    split = next(iter(winners))
    if split == 0:
        return frozenset((frozenset((a, b)), frozenset((c, d))))
    if split == 1:
        return frozenset((frozenset((a, c)), frozenset((b, d))))
    return frozenset((frozenset((a, d)), frozenset((b, c))))


def quartet_distance(t1, t2, *, max_leaves: int = 100) -> QuartetResult:
    """Quartet distance: the number of four-leaf subsets whose induced unrooted topology
    differs between the two trees (a quartet resolved differently, or resolved in one tree and
    unresolved in the other, counts as differing).

    Exact and ``O(n⁴)`` in the number of leaves ``n``; ``max_leaves`` (default 100) guards
    against accidentally launching it on a large tree — raise it explicitly to override.
    """
    t1, t2 = _as_tree(t1), _as_tree(t2)
    leaves = sorted(_check_same_leaves(t1, t2))
    n = len(leaves)
    if n > max_leaves:
        raise ValueError(
            f"quartet distance is O(n^4); this tree has {n} leaves > max_leaves={max_leaves}. "
            f"Pass max_leaves={n} (or higher) to run it anyway."
        )
    if n < 4:
        return QuartetResult(0, 0, 0.0)
    P1, P2 = _leaf_root_paths(t1), _leaf_root_paths(t2)
    differing = 0
    total = 0
    for a, b, c, d in combinations(leaves, 4):
        total += 1
        if _quartet_topology(a, b, c, d, P1) != _quartet_topology(a, b, c, d, P2):
            differing += 1
    return QuartetResult(differing, total, (differing / total) if total else 0.0)


def _clade_bitmasks(tree: Tree, index: dict) -> dict:
    """``{node: int bitmask of its leaves}`` (post-order), leaves numbered by ``index``."""
    masks: dict = {}
    for node in reversed(tree.nodes()):
        if node.is_leaf():
            masks[node] = 1 << index[node.name]
        else:
            m = 0
            for c in node.children:
                m |= masks[c]
            masks[node] = m
    return masks


def _split_bitmasks(tree: Tree, index: dict, n: int, full: int, rooted: bool) -> list:
    """Informative splits (unrooted, canonicalised to the ref-free side) or clusters (rooted),
    as a de-duplicated list of bitmasks."""
    out: set = set()
    for m in _clade_bitmasks(tree, index).values():
        size = m.bit_count()
        if rooted:
            if 1 < size < n:                         # informative cluster
                out.add(m)
        elif min(size, n - size) >= 2:               # informative bipartition
            out.add(m if not (m & 1) else full ^ m)  # drop bit 0 → canonical side
    return list(out)


def matching_distance(t1, t2, *, rooted: bool = True, max_leaves: int = 2000) -> MatchingResult:
    """Matching split (unrooted) / cluster (rooted) distance between two trees.

    Each split/cluster of one tree is matched to at most one of the other; a matched pair costs
    the number of leaves that must move to make them identical (``|A₁ △ A₂|`` for clusters, the
    min over both orientations for splits), and an unmatched split pays its cost against the
    empty split. The total is the minimum over all matchings — the assignment problem, solved
    with :func:`scipy.optimize.linear_sum_assignment`. Unlike RF this degrades gracefully: a
    slightly displaced clade costs a little, not a full unit.

    ``rooted=True`` (default) matches clades (the *matching cluster* distance); ``rooted=False``
    matches bipartitions (the *matching split* distance). ``O(n³)``, guarded by ``max_leaves``.

    Raises :class:`ImportError` if SciPy is not installed (it is an optional dependency — the
    rest of :mod:`zombi2.tools.treedist` needs only NumPy).
    """
    t1, t2 = _as_tree(t1), _as_tree(t2)
    leaves = sorted(_check_same_leaves(t1, t2))
    n = len(leaves)
    if n > max_leaves:
        raise ValueError(
            f"matching distance is O(n^3); this tree has {n} leaves > max_leaves={max_leaves}. "
            f"Pass max_leaves={n} (or higher) to run it anyway."
        )
    try:
        from scipy.optimize import linear_sum_assignment
    except ImportError as e:
        raise ImportError(
            "the matching distance needs SciPy (an optional dependency); install it with "
            "'pip install scipy' or the 'zombi2[dev]' extra."
        ) from e
    import numpy as np

    index = {name: i for i, name in enumerate(leaves)}
    full = (1 << n) - 1
    s1 = _split_bitmasks(t1, index, n, full, rooted)
    s2 = _split_bitmasks(t2, index, n, full, rooted)

    def null_weight(m: int) -> int:
        size = m.bit_count()
        return size if rooted else min(size, n - size)

    def split_cost(a: int, b: int) -> int:
        if rooted:
            return (a ^ b).bit_count()
        return min((a ^ b).bit_count(), (a ^ (full ^ b)).bit_count())

    max_distance = sum(null_weight(m) for m in s1) + sum(null_weight(m) for m in s2)
    if not s1 and not s2:
        return MatchingResult(0, 0, 0.0)

    # Pad the smaller split set with empty splits to a square N×N assignment: rows are the
    # larger (all real) set, columns are the smaller set's reals then empty columns; matching a
    # real split to an empty column pays that split's null weight.
    rows, cols_real = (s2, s1) if len(s1) <= len(s2) else (s1, s2)
    N = len(rows)
    cost = np.empty((N, N))
    for i, rm in enumerate(rows):
        rw = null_weight(rm)
        for j in range(N):
            cost[i, j] = split_cost(rm, cols_real[j]) if j < len(cols_real) else rw
    ri, ci = linear_sum_assignment(cost)
    distance = int(round(float(cost[ri, ci].sum())))
    return MatchingResult(distance, max_distance,
                          (distance / max_distance) if max_distance else 0.0)


def compare_trees(t1, t2, *, quartet: bool = True, matching: bool = True, max_leaves: int = 100,
                  matching_max_leaves: int = 2000, branch_score_order: int = 2) -> TreeComparison:
    """Every metric between two trees in one call.

    Computes rooted RF, unrooted RF, and the branch score always. The quartet distance is
    computed only when ``quartet=True`` and the tree has at most ``max_leaves`` tips; the
    (rooted) matching distance only when ``matching=True``, the tree has at most
    ``matching_max_leaves`` tips, and SciPy is installed. A skipped metric's fields are ``None``
    — reported, never silently dropped.
    """
    t1, t2 = _as_tree(t1), _as_tree(t2)
    leaves = _check_same_leaves(t1, t2)
    n = len(leaves)
    rooted = robinson_foulds(t1, t2, rooted=True)
    unrooted = robinson_foulds(t1, t2, rooted=False)
    bs = branch_score(t1, t2, order=branch_score_order)
    if quartet and n <= max_leaves:
        q = quartet_distance(t1, t2, max_leaves=max_leaves)
        qd, qn = q.differing, q.normalized
    else:
        qd = qn = None
    md = mn = None
    if matching and n <= matching_max_leaves:
        try:
            m = matching_distance(t1, t2, rooted=True, max_leaves=matching_max_leaves)
            md, mn = m.distance, m.normalized
        except ImportError:
            pass                                     # SciPy absent → matching skipped (fields None)
    return TreeComparison(
        n_leaves=n,
        rf=rooted.rf, rf_max=rooted.max_rf, rf_normalized=rooted.normalized,
        rf_unrooted=unrooted.rf, rf_unrooted_max=unrooted.max_rf,
        branch_score=bs, quartet=qd, quartet_normalized=qn,
        matching=md, matching_normalized=mn,
    )
