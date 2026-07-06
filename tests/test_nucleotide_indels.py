"""Intergenic indels in the nucleotide model.

An **insertion** lays down a run of novel nucleotides (a fresh source, its own block) inside an
intergene stretch; a **deletion** removes a run from *within a single* intergene stretch. Both are
per-nucleotide rates, off by default. The invariants pinned here:

* insertions only grow a genome, deletions only shrink it (bounded by the min-genome floor);
* in genic mode indels never split, span, or delete a gene (they stay inside one intergene run);
* the run length is geometric with mean ``indel_mean_length`` (a knob separate from ``extension``);
* the feature is deterministic and the Rust profiles fast path refuses it.
"""

import numpy as np
import pytest

from zombi2 import BirthDeath, simulate_species_tree
from zombi2.genomes.genome import IdManager
from zombi2.genomes.events import EventType, TargetParams
from zombi2.genomes.nucleotide_genome import MIN_GENOME_LENGTH, NucleotideGenome, SegmentRegistry
from zombi2.genomes.nucleotide_sim import simulate_nucleotide_genomes


def _tree(n_tips=8, age=2.0, seed=0):
    return simulate_species_tree(BirthDeath(birth=1.0, death=0.25), n_tips=n_tips, age=age, seed=seed)


def _genic(root_length=200, genes=((30, 60, "geneA"),), **kw):
    reg = SegmentRegistry(pending_genes=[tuple(g) for g in genes])
    g = NucleotideGenome(IdManager(), root_length=root_length, extension=kw.pop("ext", 0.9),
                         registry=reg, **kw)
    g.originate(np.random.default_rng(0), TargetParams())
    return g


# --------------------------------------------------------------------------- #
# Growth / shrinkage direction
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("seed", range(6))
def test_insertion_only_never_shrinks(seed):
    tree = _tree(seed=seed)
    res = simulate_nucleotide_genomes(tree, output="genomes", root_length=1000, extension=0.98,
                                      insertion=0.002, indel_mean_length=15, seed=seed * 3 + 1)
    assert any(r.event is EventType.INSERTION for r in res.event_log)
    assert all(g._length >= 1000 for g in res.leaf_genomes.values())


@pytest.mark.parametrize("seed", range(6))
def test_deletion_only_never_grows_and_respects_floor(seed):
    tree = _tree(seed=seed)
    res = simulate_nucleotide_genomes(tree, output="genomes", root_length=1000, extension=0.98,
                                      deletion=0.003, indel_mean_length=20, seed=seed * 5 + 2)
    assert any(r.event is EventType.DELETION for r in res.event_log)
    # deletions only shrink, and never below the floor (checked on every node, not just leaves)
    for g in res.node_genomes.values():
        assert MIN_GENOME_LENGTH <= g._length <= 1000


def test_deletion_cannot_empty_a_tiny_genome():
    """Heavy deletion on a tiny chromosome must stall at the floor, not empty it (which would
    stop the size-proportional process)."""
    tree = _tree(n_tips=10, age=4.0, seed=1)
    res = simulate_nucleotide_genomes(tree, output="genomes", root_length=8, extension=0.5,
                                      deletion=0.5, indel_mean_length=50, seed=9)
    for g in res.node_genomes.values():
        assert g._length >= MIN_GENOME_LENGTH


# --------------------------------------------------------------------------- #
# Genic invariant: indels never touch a gene
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("seed", range(20))
def test_indels_never_split_span_or_delete_a_gene(seed):
    tree = _tree(seed=seed)
    genes = [(50, 90, "a"), (140, 200, "b"), (300, 340, "c"), (600, 660, "d")]
    res = simulate_nucleotide_genomes(
        tree, output="genomes", root_length=1000, extension=0.9,
        insertion=0.01, deletion=0.01, indel_mean_length=25,   # indels ONLY (isolate them)
        gene_intervals=genes, seed=seed * 7 + 3)
    assert any(r.event in (EventType.INSERTION, EventType.DELETION) for r in res.event_log)
    gene_iv = {(gi.source, gi.start, gi.end): gi.gene_id
               for gis in res.registry.genes.values() for gi in gis}
    for a in res.blocks:
        if a.kind == "gene":
            assert (a.source, a.start, a.end) in gene_iv
        for (src, s, e) in gene_iv:                 # no block straddles a gene boundary
            if a.source == src and not (a.end <= s or a.start >= e):
                assert a.start >= s and a.end <= e
    # every surviving gene copy keeps its full length — never internally cut
    for genome in res.leaf_genomes.values():
        for seg in genome._segments:
            if seg.gene_id is None:
                continue
            gi = next(g for gs in res.registry.genes.values() for g in gs
                      if g.gene_id == seg.gene_id)
            assert seg.src_end - seg.src_start == gi.length


def test_intergene_run_bounds_and_gene_interior_returns_none():
    g = _genic(200, [(30, 60, "geneA"), (100, 140, "geneB")])
    # runs: [0,30) | gene | [60,100) | gene | [140,200)
    assert g._intergene_run_at(10) == (0, 30)
    assert g._intergene_run_at(80) == (60, 100)
    assert g._intergene_run_at(170) == (140, 200)
    assert g._intergene_run_at(45) is None       # strictly inside geneA
    assert g._intergene_run_at(120) is None      # strictly inside geneB


# --------------------------------------------------------------------------- #
# Run length: geometric with mean indel_mean_length
# --------------------------------------------------------------------------- #
def test_draw_indel_length_matches_mean():
    g = _genic(200, indel_mean_length=12.0)
    rng = np.random.default_rng(0)
    draws = [g._draw_indel_length(rng) for _ in range(20000)]
    assert all(d >= 1 for d in draws)
    assert abs(np.mean(draws) - 12.0) < 0.6      # geometric mean == indel_mean_length


def test_indel_mean_length_one_gives_single_nucleotides():
    g = _genic(200, indel_mean_length=1.0)
    rng = np.random.default_rng(0)
    assert all(g._draw_indel_length(rng) == 1 for _ in range(1000))


# --------------------------------------------------------------------------- #
# Determinism, default-off, and the fast-path gate
# --------------------------------------------------------------------------- #
def test_indels_are_deterministic():
    tree = _tree(seed=4)
    kw = dict(output="genomes", root_length=1000, extension=0.97,
              insertion=0.002, deletion=0.002, indel_mean_length=14, seed=123)
    a = simulate_nucleotide_genomes(tree, **kw)
    b = simulate_nucleotide_genomes(tree, **kw)
    assert [r.event for r in a.event_log] == [r.event for r in b.event_log]
    assert (sorted(g._length for g in a.leaf_genomes.values())
            == sorted(g._length for g in b.leaf_genomes.values()))


def test_default_off_no_indel_events():
    tree = _tree(seed=2)
    res = simulate_nucleotide_genomes(tree, output="genomes", root_length=500, extension=0.95,
                                      loss=0.01, duplication=0.01, seed=1)
    assert not any(r.event in (EventType.INSERTION, EventType.DELETION) for r in res.event_log)


def test_profiles_fast_path_rejects_indels():
    tree = _tree(seed=1)
    with pytest.raises(ValueError, match="indel"):
        simulate_nucleotide_genomes(tree, output="profiles", insertion=0.01, seed=1)
    with pytest.raises(ValueError, match="indel"):
        simulate_nucleotide_genomes(tree, output="profiles", deletion=0.01, seed=1)


def test_negative_indel_rate_rejected():
    from zombi2.genomes.rates import SharedRates
    with pytest.raises(ValueError, match="insertion"):
        SharedRates(insertion=-1.0)
    with pytest.raises(ValueError, match="deletion"):
        SharedRates(deletion=-1.0)
