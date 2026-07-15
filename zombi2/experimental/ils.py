"""Incomplete lineage sorting via the multispecies coalescent (**experimental**).

ZOMBI2's core gene-family engine forces every gene lineage to coalesce *exactly* at the species
divergence it passes through: at a speciation the parent copy is cloned into both daughters, so the
reconstructed gene tree always matches the species tree topology. That is the *no-ILS* limit. This
module relaxes it.

Incomplete lineage sorting (ILS) is modelled by the **multispecies coalescent** (MSC). Reading time
backward, the ``k`` gene lineages inside a branch coalesce as a Kingman coalescent at pairwise rate
``1 / N`` per unit branch time (total rate ``k(k-1) / (2N)``); any that fail to coalesce before the
branch's older (parent) end *escape* into the parent branch, where they may coalesce with lineages
from the *sister* branch. When they do, the gene tree disagrees with the species tree — ILS.

The single knob is ``population_size`` ``N``, in the tree's own time units: the amount of ILS is
governed by ``branch_length / N`` -- large ``N`` (or short branches) => slow coalescence => much
discordance; small ``N`` (or long branches) => the no-ILS limit. For the classic rooted triple
``((A,B),C)`` with internal branch length ``T`` the gene tree matches the species tree with
probability ``1 - (2/3)*e^{-T/N}`` (:func:`expected_triple_concordance`), the two discordant
resolutions being equally likely.

**Container seam (why this is v2-ready).** :meth:`MultispeciesCoalescent.sample_gene_tree` runs the
coalescent inside an arbitrary *container* :class:`~zombi2.tree.Tree`. For plain ILS the container is
the **species tree** (this module, v1). Feed the same routine a **locus tree** -- a DTL gene family
reconstructed as a ``Tree`` -- and it yields gene trees under *DTL + ILS*: duplication nodes act as
coalescent bounds exactly like speciations, and losses are simply branches with no sampled
descendant. (Only transfers fall outside the seam, since they make the container a network rather
than a tree; that is the planned v2.)

Pure ``numpy``; no optional dependencies.
"""
from __future__ import annotations

import math

import numpy as np

from zombi2.experimental import warn_experimental
from zombi2.tree import Tree, TreeNode

__all__ = [
    "MultispeciesCoalescent",
    "expected_triple_concordance",
    "rooted_clades",
    "is_concordant",
]


class MultispeciesCoalescent:
    """Draw gene trees under the multispecies coalescent inside a container tree (ILS).

    Parameters
    ----------
    population_size:
        The effective population size ``N`` in the container tree's own time units: the pairwise
        coalescence rate is ``1 / N`` per unit time, so ``k`` co-occurring lineages coalesce at total
        rate ``k(k-1) / (2N)``. Larger ``N`` => more ILS. Constant across the whole tree (a single
        global value; per-branch ``N`` is a planned extension).
    """

    def __init__(self, population_size: float):
        warn_experimental("MultispeciesCoalescent")
        n = float(population_size)
        if not np.isfinite(n) or n <= 0:
            raise ValueError(
                f"population_size must be a positive finite number, got {population_size!r}"
            )
        self.population_size = n

    # -- public API ---------------------------------------------------------
    def sample_gene_tree(self, container: Tree, *, samples=1, rng=None) -> Tree:
        """One gene tree drawn under the MSC within ``container``.

        ``samples`` is the number of gene copies sampled per container leaf: an ``int`` (the same
        count everywhere, default 1 -- one allele per species, the standard single-copy-ortholog
        setting) or a ``dict`` mapping leaf name -> count (leaves absent from the dict get 0).

        Returns a gene-tree :class:`~zombi2.tree.Tree` laid out on the container's own time axis, so a
        coalescence that predates the species root has ``time < 0`` -- the signature of deep
        coalescence. Tips are named ``<leaf>`` (or ``<leaf>_<i>`` when several copies are sampled from
        one leaf); internal coalescences are unnamed.
        """
        rng = _as_rng(rng)
        n = self.population_size
        entering: dict[int, list[TreeNode]] = {}
        gene_root: TreeNode | None = None

        # reversed pre-order visits every node after all of its descendants (a valid post-order), so a
        # branch is coalesced only once the lineages entering it from below are known.
        for node in reversed(container.nodes()):
            if node.is_leaf():
                lineages = _sampled_lineages(node, samples)
            else:
                lineages = []
                for child in node.children:
                    lineages.extend(entering.pop(id(child), []))
            if node.parent is None:
                # the root branch is unbounded: run the coalescent down to the MRCA
                survivors = _coalesce(lineages, node.time, None, n, rng)
                gene_root = survivors[0] if survivors else None
            else:
                duration = node.time - node.parent.time
                entering[id(node)] = _coalesce(lineages, node.time, duration, n, rng)

        if gene_root is None:
            raise ValueError("no lineages were sampled -- `samples` must place at least one copy")
        total_age = max((lf.time for lf in _iter_leaves(gene_root)), default=container.total_age)
        return Tree(gene_root, total_age)

    def sample_gene_trees(self, container: Tree, n: int, *, samples=1, rng=None) -> list[Tree]:
        """``n`` independent gene trees (each as in :meth:`sample_gene_tree`)."""
        rng = _as_rng(rng)
        return [self.sample_gene_tree(container, samples=samples, rng=rng) for _ in range(n)]

    # -- DTL + ILS (v2): the coalescent run within each family's locus tree -----------------
    def sample_family_gene_trees(self, genomes, *, samples=1, rng=None) -> dict:
        """Per gene family in a ZOMBI2 ``Genomes`` result, sample a gene tree under **DTL + ILS**: the
        bounded multispecies coalescent run *within that family's locus tree* (its duplication / loss /
        transfer history). Origination, a duplication's new copy, and a transferred/converted copy are
        single-copy foundings (the bounded coalescent); speciations allow deep coalescence (ILS).

        Returns ``{family: Tree}`` (families with no surviving copy are omitted). ``samples`` is the
        number of alleles sampled per extant gene copy (default 1). Same container seam as
        :meth:`sample_gene_tree` — only the container (a locus tree) and the per-event routing differ.
        """
        fams = self._family_trees(genomes.gene_families, genomes._gid_to_species(),
                                  genomes.species_tree.total_age, samples, 1, _as_rng(rng))
        return {family: trees[0] for family, trees in fams.items()}

    def _family_trees(self, families, gid2species, total_age, samples, replicates, rng) -> dict:
        """Shared engine for :meth:`sample_family_gene_trees` and the CLI: per family, build the locus
        tree once and draw ``replicates`` gene trees under DTL + ILS. Returns ``{family: list[Tree]}``
        (families with no surviving copy omitted). ``families`` is ``{family: list[EventRecord]}`` and
        ``gid2species`` is ``{gid: extant species}`` — exactly what a ``Genomes`` result or a written
        ``events_trace.tsv`` (via ``read_events_trace`` + ``extant_species_from_records``) provides."""
        from zombi2.genomes.reconciliation import _node_tree     # lazy: keep the module import light
        out: dict = {}
        for family, records in families.items():
            locus = _node_tree(records, gid2species, total_age)
            if locus is None:
                continue
            trees = [Tree(root, total_age) for root in
                     (_sample_locus_gene_tree(locus, samples, self.population_size, rng)
                      for _ in range(replicates)) if root is not None]
            if trees:
                out[family] = trees
        return out


# --------------------------------------------------------------------------- #
# the censored coalescent within one branch
# --------------------------------------------------------------------------- #
def _coalesce(lineages, start_time, duration, n, rng):
    """Coalesce ``lineages`` backward from ``start_time`` over ``duration`` (``None`` = unbounded).

    Returns the survivors -- the lineages that had not coalesced when the branch's older end was
    reached (all of them fused into one MRCA when ``duration`` is ``None``). Merges create new internal
    :class:`TreeNode`s in place.
    """
    lineages = list(lineages)
    elapsed = 0.0
    while len(lineages) >= 2:
        k = len(lineages)
        rate = k * (k - 1) / (2.0 * n)
        elapsed += float(rng.exponential(1.0 / rate))
        if duration is not None and elapsed >= duration:
            break
        i, j = (int(x) for x in rng.choice(k, size=2, replace=False))
        parent = TreeNode(name="", time=start_time - elapsed, is_extant=False)
        parent.add_child(lineages[i])
        parent.add_child(lineages[j])
        lineages = [ln for t, ln in enumerate(lineages) if t != i and t != j]
        lineages.append(parent)
    return lineages


def _sampled_lineages(leaf: TreeNode, samples) -> list[TreeNode]:
    count = int(samples.get(leaf.name, 0)) if isinstance(samples, dict) else int(samples)
    if count < 0:
        raise ValueError(f"samples per leaf must be >= 0, got {count} for {leaf.name!r}")
    if count == 1:
        return [TreeNode(name=leaf.name, time=leaf.time, is_extant=True)]
    return [
        TreeNode(name=f"{leaf.name}_{i}", time=leaf.time, is_extant=True)
        for i in range(1, count + 1)
    ]


def _iter_leaves(root: TreeNode):
    stack = [root]
    while stack:
        node = stack.pop()
        if node.children:
            stack.extend(node.children)
        else:
            yield node


def _as_rng(rng) -> np.random.Generator:
    if isinstance(rng, np.random.Generator):
        return rng
    return np.random.default_rng(rng)  # None -> fresh; otherwise treated as a seed


# --------------------------------------------------------------------------- #
# DTL + ILS: the bounded ("founder") coalescent and the locus-tree traversal (v2)
# --------------------------------------------------------------------------- #
# A locus is *founded by a single copy* at three kinds of event: the family origination, the new copy
# of a duplication, and the transferred/converted copy. Its sampled alleles must therefore coalesce
# to that one founder by the event time -- a coalescent *conditioned to reach a single lineage within
# the branch* (the bounded coalescent). A speciation is NOT a founding event (the lineages existed in
# the ancestral population and may coalesce deeper -- ILS). The bounded coalescent is what keeps the
# model causal: no allele is older than the locus it belongs to.
_FOUNDER_PARENT_EVENTS: frozenset = frozenset()


def _founder_events() -> frozenset:
    """The parent event types whose second child (`children[1]`) is a single-copy founding."""
    global _FOUNDER_PARENT_EVENTS
    if not _FOUNDER_PARENT_EVENTS:
        from zombi2.genomes.events import EventType  # lazy: keep module import numpy-only
        _FOUNDER_PARENT_EVENTS = frozenset(
            {EventType.DUPLICATION, EventType.TRANSFER, EventType.CONVERSION})
    return _FOUNDER_PARENT_EVENTS


def _g_lineages_to_one(m: int, t: float, n: float) -> float:
    """P(m lineages coalesce to a single common ancestor within time t) at population size n
    (Tavaré 1984, the j=1 case): monotone increasing in t, 0 at t=0 for m>=2, ->1 as t->inf."""
    if m <= 1:
        return 1.0
    total = 0.0
    ratio = 1.0                       # prod_{a=0}^{i-1} (m-a)/(m+a) == m_falling_i / m_rising_i
    for i in range(1, m + 1):
        ratio *= (m - (i - 1)) / (m + (i - 1))
        total += math.exp(-i * (i - 1) * t / (2.0 * n)) * (2 * i - 1) * (-1) ** (i - 1) * ratio
    return min(1.0, max(0.0, total))


def _truncated_exponential(rate: float, s: float, rng) -> float:
    """A draw from Exp(rate) conditioned to fall in [0, s]."""
    if rate * s > 700.0:              # truncation numerically irrelevant
        return float(rng.exponential(1.0 / rate))
    u = float(rng.random())
    return -math.log1p(-u * (1.0 - math.exp(-rate * s))) / rate


def _sample_bounded_wait(rate: float, s: float, m_after: int, n: float, rng) -> float:
    """The next coalescence wait for the bounded coalescent: density ∝ rate·e^{-rate·w}·g_{m_after,1}(s-w)
    on [0, s], sampled by rejection against the truncated exponential."""
    if s <= 0.0:
        return 0.0                    # no budget left -> force immediate coalescence
    if m_after <= 1:
        return _truncated_exponential(rate, s, rng)     # last merge: g_{1,1} == 1, plain truncated Exp
    gs = _g_lineages_to_one(m_after, s, n)
    if gs <= 0.0:
        return 0.0
    for _ in range(10000):
        w = _truncated_exponential(rate, s, rng)
        if float(rng.random()) * gs <= _g_lineages_to_one(m_after, max(0.0, s - w), n):
            return w
    return _truncated_exponential(rate, s, rng)          # numerical fallback (extreme conditioning)


def _coalesce_bounded(lineages, start_time, duration, n, rng):
    """Coalesce `lineages` backward from `start_time`, conditioned to reach ONE lineage within
    `duration` (the founder bottleneck). Returns [founder] (or the input unchanged if <= 1 lineage)."""
    lineages = list(lineages)
    elapsed = 0.0
    while len(lineages) >= 2:
        m = len(lineages)
        rate = m * (m - 1) / (2.0 * n)
        elapsed += _sample_bounded_wait(rate, duration - elapsed, m - 1, n, rng)
        i, j = (int(x) for x in rng.choice(m, size=2, replace=False))
        parent = TreeNode(name="", time=start_time - elapsed, is_extant=False)
        parent.add_child(lineages[i])
        parent.add_child(lineages[j])
        lineages = [ln for t, ln in enumerate(lineages) if t != i and t != j]
        lineages.append(parent)
    return lineages


def _seed_copies(node, samples) -> list[TreeNode]:
    """Sampled allele lineages at one extant gene copy (locus-tree leaf), named ``<species>_<gid>``."""
    label = f"{node.species}_{node.gid}"
    if isinstance(samples, dict):
        count = int(samples.get(label, samples.get(node.species, 1)))
    else:
        count = int(samples)
    if count <= 1:
        return [TreeNode(name=label, time=node.end, is_extant=True)]
    return [TreeNode(name=f"{label}_{i}", time=node.end, is_extant=True) for i in range(1, count + 1)]


def _sample_locus_gene_tree(root, samples, n, rng):
    """One gene tree under DTL + ILS: the bounded MSC run within one family's locus tree ``root`` (a
    reconciliation ``_Node``). Returns the gene-tree root :class:`TreeNode`, or None if no extant copy."""
    parent_of: dict[int, tuple] = {}
    pre = []
    stack = [root]
    while stack:
        nd = stack.pop()
        pre.append(nd)
        for idx, ch in enumerate(nd.children):
            parent_of[id(ch)] = (nd, idx)
        stack.extend(reversed(nd.children))

    founder_parents = _founder_events()
    out: dict[int, list] = {}
    gene_root = None
    for nd in reversed(pre):          # every node after all its descendants
        if nd.children:
            entering = []
            for ch in nd.children:
                entering.extend(out.pop(id(ch), []))
        elif nd.is_extant and nd.species is not None:
            entering = _seed_copies(nd, samples)
        else:
            entering = []             # a loss / non-extant leaf contributes nothing
        info = parent_of.get(id(nd))
        is_founder = info is None or (info[0].kind in founder_parents and info[1] == 1)
        duration = nd.end - nd.birth
        if is_founder:
            result = _coalesce_bounded(entering, nd.end, duration, n, rng)
        else:
            result = _coalesce(entering, nd.end, duration, n, rng)
        if info is None:
            gene_root = result[0] if result else None
        else:
            out[id(nd)] = result
    return gene_root


# --------------------------------------------------------------------------- #
# analytics / diagnostics (used by the CLI summary and the tests)
# --------------------------------------------------------------------------- #
def expected_triple_concordance(internal_branch_length: float, population_size: float) -> float:
    """P(gene tree == species tree) for a rooted triple ``((A,B),C)`` with internal branch length
    ``T`` and population size ``N``: ``1 - (2/3)*e^{-T/N}`` (Hudson 1983; Nei 1987). The two
    discordant topologies are each ``(1/3)*e^{-T/N}``."""
    return 1.0 - (2.0 / 3.0) * math.exp(-internal_branch_length / population_size)


def rooted_clades(tree: Tree) -> set[frozenset[str]]:
    """The set of rooted clades of ``tree``, each a ``frozenset`` of descendant leaf names. Two trees
    on the same leaf set have identical rooted topology iff their clade sets are equal."""
    clades: set[frozenset[str]] = set()

    def rec(node: TreeNode) -> frozenset[str]:
        if node.is_leaf():
            return frozenset((node.name,))
        leaves = frozenset().union(*(rec(c) for c in node.children))
        clades.add(leaves)
        return leaves

    rec(tree.root)
    return clades


def is_concordant(gene_tree: Tree, species_tree: Tree) -> bool:
    """True iff ``gene_tree`` and ``species_tree`` have identical rooted topology. Meaningful only
    when the two share a leaf set (i.e. one copy sampled per species)."""
    return rooted_clades(gene_tree) == rooted_clades(species_tree)
