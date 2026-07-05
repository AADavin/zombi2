"""Family Sequence Evolution — reconciled gene trees from time into substitutions/site.

ZOMBI2 gene trees are **timetrees** (branch lengths are time). Real sequence data would show
branch lengths in *substitutions per site*, which differ from time because lineages and gene
families evolve at different speeds. This module overlays a **gene × lineage** substitution-rate
model on the reconciled gene trees:

.. math:: \\text{rate}(\\text{family } g,\\ \\text{species branch } b) = R_b \\cdot s_g

* **Lineage rate** ``R_b`` — an autocorrelated lognormal relaxed clock on the **species tree**:
  ``R_child = R_parent · exp(N(0, sigma·sqrt(branch_length)))``. There is one rate per species
  branch and it is **shared by every gene family**, so a fast-evolving clade is fast for all the
  genes passing through it. ``sigma = 0`` gives a strict clock (every branch rate ``root_rate``).
* **Family speed** ``s_g`` — each gene family draws **one constant** speed multiplier from a
  distribution (e.g. ``LogNormal(0, sigma)``), so some families are globally fast, others slow.

A gene-tree branch that lives on species branch ``b`` over the time interval ``[t0, t1]`` gets
substitution length ``s_g · R_b · (t1 - t0)``. Because ZOMBI2 reconciliation is exact, each
complete-tree branch lies on a single species branch; after pruning to the extant tree a branch
may span several species branches, and the pieces simply **sum**. With ``sigma = 0`` and unit
family speeds the phylogram is identical to the input chronogram.

The lineage clock has two forms: the per-branch lognormal relaxed clock (``branch_sigma`` /
``--branch-speed``, constant along a branch), or the discrete-bin within-branch model of
:class:`~zombi2.RateVariation` (the GTDB model — pass ``lineage=RateVariation(...)`` /
``--branch-bins``), which can vary the rate *within* a branch. The family speed multiplies
whichever lineage clock you choose.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .distributions import as_distribution
from .profiles import _natkey
from .reconciliation import _Node, _node_tree


@dataclass
class GenePhylograms:
    """Result of :meth:`SequenceEvolution.scale`.

    ``complete`` / ``extant`` map each family id to its phylogram Newick (branch lengths in
    substitutions/site), or ``None`` for a family with no lineages / no survivors.
    ``family_speed`` is the drawn per-family speed ``s_g``; ``branch_rate`` is the shared
    per-species-branch lineage rate ``R_b`` (keyed by species-branch name). The last two make
    a run fully reproducible and inspectable.
    """

    complete: dict
    extant: dict
    family_speed: dict
    branch_rate: dict


class SequenceEvolution:
    """Gene × lineage substitution-rate clock over reconciled gene trees.

    The **lineage** clock has two forms — pick one:

    Parameters
    ----------
    branch_sigma:
        Drift of the autocorrelated **lognormal** lineage clock on the species tree, per
        ``sqrt(time)`` (``>= 0``): ``R_child = R_parent · exp(N(0, branch_sigma·sqrt(L)))``,
        one rate per branch. ``0`` is a strict clock (all branches share ``root_rate``).
    lineage:
        A :class:`~zombi2.RateVariation` — the **discrete-bin** within-branch clock (the GTDB
        model), used as an alternative to ``branch_sigma``. It is run once on the species tree,
        so a branch may carry several rate segments. Mutually exclusive with ``branch_sigma``.
    family_speed:
        The per-**family** speed distribution (``Distribution`` / float / scipy frozen dist /
        callable — see :func:`~zombi2.as_distribution`). Each family draws one constant, clamped
        to ``>= 0``. A bare float is a fixed speed for every family (``1.0`` = no family effect).
    root_rate:
        Lineage rate of the root branch, which the lognormal walk starts from (default ``1.0``).
    """

    def __init__(self, *, branch_sigma: float = 0.0, lineage=None, family_speed=1.0,
                 root_rate: float = 1.0):
        if lineage is not None and branch_sigma:
            raise ValueError("give either branch_sigma (lognormal per-branch clock) OR "
                             "lineage=RateVariation (discrete-bin clock), not both")
        if branch_sigma < 0:
            raise ValueError(f"branch_sigma must be >= 0, got {branch_sigma}")
        if root_rate <= 0:
            raise ValueError(f"root_rate must be > 0, got {root_rate}")
        self.branch_sigma = float(branch_sigma)
        self.lineage = lineage
        self.family_speed = as_distribution(family_speed)
        self.root_rate = float(root_rate)

    def _lineage_segments(self, tree, rng):
        """Resolve the lineage clock to ``({branch: [(rate, lo_abs, hi_abs), ...]}, {branch: avg})``.

        ``segments`` is what the integrator consumes; ``avg`` is the time-averaged rate per
        branch, reported as ``branch_rate`` (for the lognormal clock this is exactly ``R_b``).
        """
        segments: dict = {}
        avg: dict = {}
        if self.lineage is not None:
            scaled = self.lineage.scale(tree, rng=rng)      # RateScaledTree (within-branch bins)
            bins = self.lineage.bins
            for node in tree.nodes_preorder():
                if node.parent is None:
                    segments[node.name], avg[node.name] = [], self.root_rate
                    continue
                start = node.parent.time
                pieces = []
                for b, dur in scaled.segments[node]:
                    pieces.append((bins[b], start, start + dur))
                    start += dur
                segments[node.name] = pieces
                length = node.branch_length()
                avg[node.name] = (scaled.branch_lengths[node] / length if length > 0
                                  else bins[scaled.end_bin[node]])
            return segments, avg

        # lognormal per-branch clock: one constant-rate segment spanning the whole branch
        for node in tree.nodes_preorder():
            if node.parent is None:
                segments[node.name], avg[node.name] = [], self.root_rate
                continue
            scale = self.branch_sigma * math.sqrt(max(node.branch_length(), 0.0))
            drift = math.exp(rng.normal(0.0, scale)) if scale > 0 else 1.0
            r = avg[node.parent.name] * drift
            segments[node.name] = [(r, node.parent.time, node.time)]
            avg[node.name] = r
        return segments, avg

    def scale(self, genomes, *, rng: np.random.Generator | None = None,
              seed: int | None = None) -> GenePhylograms:
        """Overlay the model on a :class:`~zombi2.Genomes` result and return the phylograms.

        The lineage clock is drawn once (shared across families); each family then draws its
        own constant speed. Reproducible given ``seed``.
        """
        return self.scale_families(genomes.species_tree, genomes.gene_families,
                                   genomes._gid_to_species(), rng=rng, seed=seed)

    def scale_families(self, species_tree, families, gid2species, *,
                       rng: np.random.Generator | None = None,
                       seed: int | None = None) -> GenePhylograms:
        """Scale reconciled gene trees given the raw ingredients (not a live ``Genomes``).

        ``families`` is ``{family: [EventRecord]}`` and ``gid2species`` maps each extant gene
        id to its species. This is the entry point used by ``zombi2 sequence`` after replaying a
        written ``Events_trace.tsv``; :meth:`scale` delegates here for the in-process case.
        """
        return self._scale_families(species_tree, families, gid2species, rng=rng, seed=seed)[0]

    def scale_families_trees(self, species_tree, families, gid2species, *,
                             rng: np.random.Generator | None = None,
                             seed: int | None = None):
        """Like :meth:`scale_families` but also return the substitution-scaled node trees.

        Returns ``(phylograms, node_trees)`` where ``node_trees`` maps each family to
        ``{"complete": (root, subst), "extant": (root, subst)}`` — the reconciliation ``_Node``
        trees with a ``subst`` map (branch length in substitutions/site of the branch ending at
        each node). These are exactly what :func:`~zombi2.evolve_on_tree` consumes to simulate a
        sequence alignment along the rescaled tree. A family with no lineages / no survivors maps
        its ``complete`` / ``extant`` entry to ``None``.
        """
        return self._scale_families(species_tree, families, gid2species, rng=rng, seed=seed)

    def _scale_families(self, species_tree, families, gid2species, *, rng, seed):
        if rng is None:
            rng = np.random.default_rng(seed)

        tree = species_tree
        segments, avg = self._lineage_segments(tree, rng)
        total_age = tree.total_age

        complete: dict = {}
        extant: dict = {}
        speeds: dict = {}
        node_trees: dict = {}
        # canonical family order so a family always draws the same speed for a given seed,
        # whether we came from a live Genomes or a replayed trace (dict order may differ)
        for fam in sorted(families, key=_natkey):
            records = families[fam]
            root_node = _node_tree(records, gid2species, total_age)
            if root_node is None:
                complete[fam] = extant[fam] = None
                node_trees[fam] = {"complete": None, "extant": None}
                continue
            s_g = max(0.0, self.family_speed.sample(rng))
            speeds[fam] = s_g
            subst: dict = {}
            _annotate(root_node, segments, s_g, subst)
            complete[fam] = _phylo_newick(root_node, subst) + ";"
            pruned, psubst = _prune_subst(root_node, subst)
            extant[fam] = (_phylo_newick(pruned, psubst) + ";") if pruned is not None else None
            node_trees[fam] = {
                "complete": (root_node, subst),
                "extant": (pruned, psubst) if pruned is not None else None,
            }

        return GenePhylograms(complete, extant, speeds, avg), node_trees


def _integrate(segs, t0: float, t1: float) -> float:
    """Integrate the piecewise-constant lineage rate over ``[t0, t1]`` (∑ rate × overlap)."""
    total = 0.0
    for rate, lo, hi in segs:
        a, b = max(t0, lo), min(t1, hi)
        if b > a:
            total += rate * (b - a)
    return total


def _annotate(root: _Node, segments: dict, s_g: float, subst: dict) -> None:
    """Fill ``subst[node]`` = substitution length of each node's branch = ``s_g · ∫ R_b dt``."""
    def visit(node: _Node) -> None:
        # the species branch this gene lineage lives on over [birth, end]: the branch of its
        # own terminating event, or (for a survivor) its extant leaf species
        name = node.branch if node.branch is not None else node.species
        subst[node] = s_g * _integrate(segments.get(name, ()), node.birth, node.end)
        for child in node.children:
            visit(child)
    visit(root)


def _prune_subst(node: _Node, subst: dict):
    """Prune to extant lineages, summing substitution lengths across suppressed degree-2 nodes.

    Mirrors :func:`~zombi2.reconciliation._prune` but tracks branch lengths in the ``out`` map
    (keyed by the freshly built pruned nodes) instead of relying on ``end - birth``.
    """
    out: dict = {}

    def rec(n: _Node):
        if not n.children:
            if not n.is_extant:
                return None
            leaf = _Node(n.gid, n.birth)
            leaf.end, leaf.species, leaf.is_extant = n.end, n.species, True
            out[leaf] = subst[n]
            return leaf
        kept = [k for k in (rec(c) for c in n.children) if k is not None]
        if not kept:
            return None
        if len(kept) == 1:  # suppress: fold this node's branch into the single survivor
            survivor = kept[0]
            out[survivor] += subst[n]
            survivor.birth = n.birth
            return survivor
        inner = _Node(n.gid, n.birth)
        inner.end, inner.kind, inner.children, inner.branch = n.end, n.kind, kept, n.branch
        out[inner] = subst[n]
        return inner

    return rec(node), out


def _phylo_newick(node: _Node, subst: dict) -> str:
    """Newick with branch lengths from ``subst`` (mirrors reconciliation._to_newick's labels)."""
    if not node.children:
        if node.is_loss:
            name = f"LOSS_{node.gid}"
        elif node.is_extant:
            name = f"{node.species}_{node.gid}"
        else:
            name = node.gid
        return f"{name}:{subst[node]:.6g}"
    inner = ",".join(_phylo_newick(c, subst) for c in node.children)
    return f"({inner}){node.gid}:{subst[node]:.6g}"
