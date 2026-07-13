"""Linear (and mixed-topology) nucleotide chromosomes.

A nucleotide chromosome is a **ring** by default: an event arc ``[s, s+ell)`` may wrap the origin.
A *linear* chromosome has two ends and no wrap. The whole feature rests on one invariant — **an
arc never wraps on a linear chromosome** — enforced by clamping the arc at draw time; every event's
existing non-wrapping path then handles it. These tests are the adversarial safety net for that
invariant, and for the topology-aware tier (fission at one breakpoint, same-topology fusion) and
per-sequence seeding (including *mixed* circular + linear genomes, e.g. a linear chromosome with
circular plasmids).

Circular-only runs must stay byte-identical (that guarantee lives in test_chromosome_fingerprints.py);
here every assertion is about linear / mixed behaviour that is off unless a chromosome is made linear.
"""
from collections import Counter

import numpy as np
import pytest

from zombi2 import BirthDeath, simulate_species_tree
from zombi2.genomes.events import EventType, TargetParams
from zombi2.genomes.genome import IdManager
from zombi2.genomes.nucleotide_genome import NucleotideGenome, SegmentRegistry
from zombi2.genomes.nucleotide_sim import simulate_nucleotide_genomes

E = EventType
_ARC_EVENTS = (E.INVERSION, E.LOSS, E.DUPLICATION, E.TRANSPOSITION, E.TRANSLOCATION)


def _linear(root_length=200, genes=None):
    """A freshly-seeded genome whose (single) chromosome has been made linear."""
    reg = SegmentRegistry(pending_genes=[tuple(g) for g in genes]) if genes else None
    g = NucleotideGenome(IdManager(), root_length=root_length, extension=0.9, registry=reg)
    g.originate(np.random.default_rng(1), TargetParams())
    g.chromosomes[0].circular = False
    return g


# --------------------------------------------------------------------------- #
# L1 — the invariant: an arc never wraps on a linear chromosome
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("genes", [None, [(20, 60, "a"), (100, 150, "b"), (170, 195, "c")]])
def test_linear_arc_never_wraps(genes):
    """Every event arc drawn on a linear chromosome stays within [0, L): start + length <= L."""
    g = _linear(genes=genes)
    L = g._chrom_length(g.chromosomes[0])
    rng = np.random.default_rng(0)
    for _ in range(5000):
        for ev in _ARC_EVENTS:
            r = g.draw_target(ev, rng, TargetParams()).region
            assert r.start + r.length <= L, f"{ev}: arc [{r.start},{r.start + r.length}) wraps L={L}"


@pytest.mark.parametrize("genes", [None, [(20, 60, "a"), (100, 150, "b")]])
def test_events_apply_on_a_linear_chromosome(genes):
    """Inversion / transposition applied to a linear chromosome permute the mosaic without wrapping:
    the (source, position) multiset is conserved, ids stay unique, the chromosome stays linear."""
    g = _linear(genes=genes)
    before = Counter((s, p) for (s, p, _st) in g.to_cells())
    rng = np.random.default_rng(3)
    for _ in range(300):
        for ev in (E.INVERSION, E.TRANSPOSITION):
            sel = g.draw_target(ev, rng, TargetParams())
            g.apply(ev, sel, rng, TargetParams())
    after = Counter((s, p) for (s, p, _st) in g.to_cells())
    assert after == before                                 # only permuted — nothing created/lost
    ids = [seg.seg_id for seg in g._iter_segments()]
    assert len(ids) == len(set(ids))
    assert g.chromosomes[0].circular is False and g.size() == before.total()


# --------------------------------------------------------------------------- #
# L2 — the tier respects topology (fission at one breakpoint, same-topology fusion)
# --------------------------------------------------------------------------- #
def test_linear_fission_splits_into_two_linear_pieces():
    g = _linear(genes=[(20, 60, "a"), (100, 150, "b")])
    before = Counter(g.to_cells())
    src, new = g.fission(np.random.default_rng(4), TargetParams())
    c0, c1 = g.chromosomes[src], g.chromosomes[new]
    assert not c0.circular and not c1.circular            # both pieces stay linear
    assert len(c0) >= 1 and len(c1) >= 1                  # one breakpoint, both non-empty
    assert Counter(g.to_cells()) == before               # content conserved


def test_fusion_is_same_topology_only():
    # two linear chromosomes fuse into a linear one, content conserved
    g = _linear(genes=[(20, 60, "a"), (100, 150, "b")])
    before = Counter(g.to_cells())
    g.fission(np.random.default_rng(4), TargetParams())   # -> two linear
    keep, _ = g.fusion(np.random.default_rng(4), TargetParams())
    assert not g.chromosomes[keep].circular and Counter(g.to_cells()) == before
    # a mixed circular+linear genome cannot fuse -> fusion returns None (the sampler skips)
    g2 = NucleotideGenome(IdManager(), root_length=100, extension=0.9)
    g2.originate(np.random.default_rng(1), TargetParams())     # chromosome 0, circular
    lin = g2._new_root_chromosome(circular=False)
    lin.elements.append(g2._new_segment("Z", 0, 50, 1))
    assert g2.fusion(np.random.default_rng(0), TargetParams()) is None


# --------------------------------------------------------------------------- #
# L3 — per-sequence seeding, including MIXED circular + linear genomes
# --------------------------------------------------------------------------- #
def _tree(seed=0, n_tips=8):
    return simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=n_tips, age=2.5, seed=seed)


def test_seed_mixed_topology_from_specs():
    # a linear chromosome + a circular plasmid, each with its own genes and topology
    res = simulate_nucleotide_genomes(
        _tree(seed=3), root_chromosomes=[(400, [(30, 90, "c1"), (150, 230, "c2")], False),
                                         (200, [(20, 80, "p1")], True)],
        inversion=0.004, seed=3, retain_internal=True)
    root = list(res.node_genomes[res.species_tree.root].chromosomes.values())
    assert [c.circular for c in root] == [False, True]    # seeded with the right per-replicon topology


def test_circular_flag_makes_all_chromosomes_linear():
    res = simulate_nucleotide_genomes(_tree(seed=2), initial_chromosomes=2, circular=False,
                                      root_length=200, inversion=0.005, seed=2)
    assert all(not c.circular for g in res.leaf_genomes.values() for c in g.chromosomes.values())


def test_mixed_topology_survives_all_events():
    """The full event zoo on a MIXED circular+linear genome: it must run, keep unique ids, conserve
    each chromosome's topology through the tier, and reconstruct — the adversarial linear/mixed
    counterpart of the circular all-events regression test."""
    roots = [(400, [(30, 90, "c1"), (150, 230, "c2"), (300, 360, "c3")], False),   # linear
             (200, [(20, 80, "p1"), (120, 170, "p2")], True)]                        # circular
    fired, tier = Counter(), Counter()
    for seed in range(6):
        res = simulate_nucleotide_genomes(
            _tree(seed=seed), root_chromosomes=roots,
            inversion=0.004, transposition=0.004, translocation=0.008,
            duplication=0.004, transfer=0.004, loss=0.006,
            insertion=0.004, deletion=0.006, origination=0.2,
            fission=0.25, fusion=0.25, chromosome_origination=0.1, chromosome_loss=0.1,
            pseudogenization=0.4, replacement=0.4, extension=0.9, seed=seed)
        fired += Counter(r.event.value for r in res.event_log.records if r.event.value != "S")
        tier += Counter(r.event.value for r in res.event_log.chromosome_records)
        for g in res.leaf_genomes.values():
            ids = [s.seg_id for s in g._iter_segments()]
            assert len(ids) == len(set(ids)) and len(g.chromosomes) >= 1
            g.to_cells()                                  # assembling a mixed genome must not wrap
        assert len(res.block_gene_trees()) == len(res.blocks)
    assert fired["O"] >= 1 and tier["FI"] >= 1 and tier["FU"] >= 1   # the crash-prone combos fired


def test_all_linear_survives_all_events():
    """The full event zoo on an all-LINEAR multi-chromosome genome: it must run, keep unique ids,
    keep every chromosome linear through the tier, and reconstruct — the pure-linear counterpart of
    the circular all-events regression test."""
    fired, tier = Counter(), Counter()
    for seed in range(6):
        res = simulate_nucleotide_genomes(
            _tree(seed=seed), initial_chromosomes=3, circular=False, root_length=250,
            gene_intervals=None,
            inversion=0.004, transposition=0.004, translocation=0.008,
            duplication=0.004, transfer=0.004, loss=0.006,
            insertion=0.004, deletion=0.006, origination=0.2,
            fission=0.25, fusion=0.25, chromosome_origination=0.1, chromosome_loss=0.1,
            extension=0.9, seed=seed)
        fired += Counter(r.event.value for r in res.event_log.records if r.event.value != "S")
        tier += Counter(r.event.value for r in res.event_log.chromosome_records)
        for g in res.leaf_genomes.values():
            ids = [s.seg_id for s in g._iter_segments()]
            assert len(ids) == len(set(ids)) and len(g.chromosomes) >= 1
            # every non-plasmid chromosome stays linear (plasmid origination is always circular)
            g.to_cells()
        assert len(res.block_gene_trees()) == len(res.blocks)
    assert fired["O"] >= 1 and tier["FI"] >= 1


def test_linear_requires_the_python_engine():
    with pytest.raises(ValueError, match="linear"):
        simulate_nucleotide_genomes(_tree(seed=1), circular=False, output="profiles",
                                    extension=0.9, seed=1)
