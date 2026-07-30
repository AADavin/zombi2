"""Marker quality — can this family be used as a phylogenetic marker, and would it mislead you?

The question most people are asking when they ask for "the orthologs" is not about a pair of genes.
It is about a **family**: *can I put this one in a concatenation and trust the tree that comes out?*
No amount of pairwise labelling adds up to that answer, which is why the homology table
(`homology`) never quite gives it. This one does, family by family.

Three things decide it, and ZOMBI has all three exactly, because it recorded the history rather than
inferring it:

- **single-copy** — every genome that has the family has exactly one copy of it, so there is no
  choosing which copy to align;
- **universal** — every extant genome has it, the criterion a BUSCO-style marker set is built on;
- **congruent** — the family's true gene tree has the same shape as the species tree over the genomes
  it occupies, so aligning it recovers the species history rather than the gene's own.

The third is the one worth the exercise. A family can be single-copy **and** universal and still give
the wrong tree — a duplication followed by loss of the other copy in each descendant, or a transfer
that replaced the resident gene. That is hidden paralogy, and it is *invisible in real data*: it
passes every filter a phylogenomicist applies and quietly poisons the concatenation. Here it is a
column. Being able to say "these 200 families pass every filter you would apply, and these 12 will
give you the wrong topology" is the thing a simulator can offer that a dataset cannot.

Congruence is measured as the Robinson–Foulds distance between the family's extant gene tree, with
each gene read as the genome it sits in, and the species tree **restricted to those same genomes** —
so a family present in half the tree is judged against the half it occupies, not against the whole.
``rf = 0`` is a perfect marker. It is the same distance `zombi2 tools treedist` reports, on the same
rooted-clade convention, so the two agree.
"""

from __future__ import annotations

import collections
import pathlib

from ..genomes.gene_trees import GeneTree

#: the columns, in order — the header and the row builder read from the one list
_COLS = ("family", "genomes", "copies", "single_copy", "universal",
         "duplications", "transfers", "losses", "rf", "congruent")


def _clades(root, species_of) -> set[frozenset]:
    """The non-trivial rooted clades of a tree, as frozensets of species.

    ``2 <= size <= n-1`` — the same convention `zombi2.tree.distance()` uses for Robinson–Foulds, so a
    number here means what a number there means."""
    below: dict[int, frozenset] = {}
    order, stack = [], [root]
    while stack:
        n = stack.pop()
        order.append(n)
        stack.extend(n.children)
    for n in reversed(order):
        below[id(n)] = (frozenset({species_of(n)}) if n.is_leaf
                        else frozenset().union(*(below[id(c)] for c in n.children)))
    total = len(below[id(root)])
    return {c for c in below.values() if 2 <= len(c) <= total - 1}


def _species_clades(tree, keep: set[int]) -> set[frozenset]:
    """The species tree's clades **induced on** ``keep`` — each clade intersected with the genomes the
    family actually occupies, which is what a family present in part of the tree must be judged
    against. Restricting the clades is the same thing as pruning the tree to those tips, without
    building one."""
    below: dict[int, frozenset] = {}
    for i in sorted(tree.nodes, reverse=True):      # children have higher ids than their parent
        node = tree.nodes[i]
        below[i] = (frozenset({i} & keep) if node.children is None
                    else frozenset().union(*(below[c] for c in node.children)))
    total = len(keep)
    return {c for c in below.values() if 2 <= len(c) <= total - 1}


def marker_row(family: int, gt: GeneTree, tree) -> dict:
    """One family's marker verdict — the fields of ``_COLS`` as a dict.

    Event counts are over the family's **whole** history, dead lineages included: that is what
    happened to the family, and a duplication in a lineage that later died is part of it even though
    it cannot affect the marker. ``rf`` is the congruence check, and is left empty where it has no
    meaning — a family with several copies in one genome (no one-to-one gene→genome map) or fewer
    than three genomes (no non-trivial clade to disagree about)."""
    kinds = collections.Counter()
    leaves, stack = [], [gt.complete]
    while stack:
        n = stack.pop()
        kinds[n.kind] += 1
        if n.kind == "extant":
            leaves.append(n)
        stack.extend(n.children)

    per_genome = collections.Counter(n.species for n in leaves)
    single = bool(per_genome) and max(per_genome.values()) == 1
    n_extant = sum(1 for _ in tree.extant_leaves())
    row = {"family": family,
           "genomes": len(per_genome),
           "copies": len(leaves),
           "single_copy": single,
           "universal": len(per_genome) == n_extant,
           "duplications": kinds["duplication"],
           "transfers": kinds["transfer"],
           "losses": kinds["loss"],
           "rf": "",
           "congruent": ""}
    if single and len(per_genome) >= 3:
        rf = len(_clades(gt.extant, lambda n: n.species) ^ _species_clades(tree, set(per_genome)))
        row["rf"], row["congruent"] = rf, rf == 0
    return row


def markers_tsv(gene_trees: dict[int, GeneTree], tree) -> str:
    """The marker table as TSV — one row per family that left a surviving copy, in family order."""
    rows = ["\t".join(_COLS)]
    for family, gt in sorted(gene_trees.items()):
        if gt.extant is None:                       # no survivor: nothing to use as a marker
            continue
        row = marker_row(family, gt, tree)
        rows.append("\t".join(_fmt(row[c]) for c in _COLS))
    return "\n".join(rows) + "\n"


def _fmt(v) -> str:
    return {True: "yes", False: "no"}.get(v, v) if isinstance(v, bool) else str(v)


def write_markers(gene_trees: dict[int, GeneTree], tree, directory) -> str:
    """Write ``markers.tsv`` — one row per family — into ``directory``. One file, not one per family:
    the whole point is to sort and filter it, which needs every family in one table."""
    d = pathlib.Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    (d / "markers.tsv").write_text(markers_tsv(gene_trees, tree), encoding="utf-8")
    return "1 table"


__all__ = ["marker_row", "markers_tsv", "write_markers"]
