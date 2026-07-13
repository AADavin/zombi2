"""Chromosome-tier events for the nucleotide model: fission, fusion, plasmid origination, loss.

These are one tier above the segment events — they act on whole chromosomes. The invariants
they must respect:

* **Content conservation** — fission and fusion only move segments between chromosomes; the
  multiset of trace-back cells (and every segment id) is preserved. Only chromosome *loss*
  removes material.
* **Genes stay whole** — in genic mode a fission breakpoint snaps to a segment boundary, so no
  gene is ever split across the two daughter chromosomes.
* **Reconstruction is chromosome-agnostic** — blocks tile by ancestral coordinate, so per-block
  gene trees build identically whether the surviving segments sit on one chromosome or many.
* **Off by default** — all four rates default to 0, so a genome stays single-chromosome unless
  asked otherwise (the byte-identity guarantee lives in test_chromosome_fingerprints.py).
"""
from collections import Counter

import numpy as np
import pytest

from zombi2 import BirthDeath, simulate_species_tree
from zombi2.genomes.events import EventType, Region, Selection, TargetParams
from zombi2.genomes.genome import IdManager
from zombi2.genomes.nucleotide_genome import NucleotideGenome, SegmentRegistry
from zombi2.genomes.nucleotide_sim import simulate_nucleotide_genomes


def _tree(n_tips=8, age=3.0, seed=0):
    return simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=n_tips, age=age, seed=seed)


def _fresh(length=300, ext=0.9, genes=None):
    reg = SegmentRegistry(pending_genes=[tuple(g) for g in genes]) if genes else None
    g = NucleotideGenome(IdManager(), root_length=length, extension=ext, registry=reg)
    g.originate(np.random.default_rng(0), TargetParams())  # seed the root chromosome
    return g


def _cells(genome):
    return Counter(genome.to_cells())


# --------------------------------------------------------------------------- #
# Unit level: the four methods on a single genome
# --------------------------------------------------------------------------- #
def test_originate_chromosome_adds_empty_circular():
    g = _fresh()
    before = _cells(g)
    cid = g.originate_chromosome(np.random.default_rng(1), TargetParams())
    assert cid in g.chromosomes and cid != 0
    assert len(g.chromosomes[cid]) == 0 and g.chromosomes[cid].circular
    assert len(g.chromosomes) == 2
    assert _cells(g) == before  # a plasmid starts empty — content unchanged


def test_fission_splits_into_two_conserving_content():
    g = _fresh()
    before = _cells(g)
    src, new = g.fission(np.random.default_rng(3), TargetParams())
    assert src == 0 and new != 0 and new in g.chromosomes
    assert len(g.chromosomes) == 2
    assert g.chromosomes[new].circular
    # both daughters carry material and their union is exactly the parent's content
    assert len(g.chromosomes[0]) >= 1 and len(g.chromosomes[new]) >= 1
    assert _cells(g) == before


def test_fission_splits_at_most_two_boundaries():
    # fission cuts the ring at two breakpoints; only the (up to two) segments straddling them are
    # re-minted (split, with provenance links) — everything else keeps its id and its content.
    g = _fresh()
    before_n = g.n_segments()
    before_cells = _cells(g)
    g.fission(np.random.default_rng(5), TargetParams())
    assert g.n_segments() <= before_n + 2
    assert _cells(g) == before_cells  # no material created or destroyed by the cuts


def test_fusion_merges_two_into_one_conserving_content():
    g = _fresh()
    g.fission(np.random.default_rng(7), TargetParams())
    assert len(g.chromosomes) == 2
    before = _cells(g)
    keep, absorbed = g.fusion(np.random.default_rng(7), TargetParams())
    assert absorbed not in g.chromosomes and keep in g.chromosomes
    assert len(g.chromosomes) == 1
    assert _cells(g) == before  # fusion only concatenates — content unchanged


def test_lose_chromosome_drops_it_and_reports_losses():
    g = _fresh()
    _, new = g.fission(np.random.default_rng(9), TargetParams())
    lost_segments = {s.seg_id for s in g.chromosomes[new].elements}
    # bias uniform choice so we lose the freshly split-off chromosome deterministically enough:
    # try seeds until the second chromosome is the one dropped
    for seed in range(50):
        h = _fresh()
        _, n2 = h.fission(np.random.default_rng(9), TargetParams())
        expect = {s.seg_id for s in h.chromosomes[n2].elements}
        cid, groups = h.lose_chromosome(np.random.default_rng(seed))
        if cid == n2:
            assert cid not in h.chromosomes and len(h.chromosomes) == 1
            gone = {op.gid for grp in groups for op in grp}
            assert gone == expect  # every segment on the lost chromosome is reported lost
            assert all(op.role == "lost" for grp in groups for op in grp)
            break
    else:
        pytest.fail("uniform loss never selected the second chromosome in 50 tries")
    assert lost_segments  # sanity: the split-off chromosome had content to lose


def test_fission_keeps_genes_whole():
    genes = [(20, 40, "geneA"), (60, 75, "geneB"), (120, 160, "geneC")]
    for seed in range(40):
        g = _fresh(length=200, genes=genes)
        g.fission(np.random.default_rng(seed), TargetParams())
        per_gene = Counter(s.gene_id for s in g._iter_segments() if s.gene_id is not None)
        # each gene remains exactly one whole segment — never split across the breakpoint
        assert per_gene == {"geneA": 1, "geneB": 1, "geneC": 1}


# --------------------------------------------------------------------------- #
# Initial karyotype: seeding N root chromosomes
# --------------------------------------------------------------------------- #
def test_initial_chromosomes_seeds_independent_copies():
    g = NucleotideGenome(IdManager(), root_length=100, extension=0.9, initial_chromosomes=3)
    g.originate(np.random.default_rng(0), TargetParams())  # one seed call lays down all three
    assert len(g.chromosomes) == 3
    assert all(len(c) == 1 and c.circular for c in g.chromosomes.values())
    assert g.size() == 300                       # three independent root-length copies
    assert sorted(g.families()) == ["1", "2", "3"]  # each under its own source namespace


def test_initial_chromosomes_genic_copies_share_layout_under_own_source():
    reg = SegmentRegistry(pending_genes=[(20, 40, "geneA"), (60, 75, "geneB")])
    g = NucleotideGenome(IdManager(), root_length=100, extension=0.9,
                         registry=reg, initial_chromosomes=3)
    g.originate(np.random.default_rng(0), TargetParams())
    assert len(g.chromosomes) == 3
    for chrom in g.chromosomes.values():
        layout = [(s.src_start, s.src_end, s.gene_id) for s in chrom.elements]
        assert layout == [(0, 20, None), (20, 40, "geneA"), (40, 60, None),
                          (60, 75, "geneB"), (75, 100, None)]
        assert len({s.source for s in chrom.elements}) == 1  # one source per chromosome
    assert len(g.families()) == 3                            # three distinct sources


def test_single_seed_call_regardless_of_count():
    # seeding is idempotent after the first origination: a second call adds a novel gene, not a copy
    g = NucleotideGenome(IdManager(), root_length=100, extension=0.9, initial_chromosomes=2)
    g.originate(np.random.default_rng(0), TargetParams())
    assert len(g.chromosomes) == 2 and g.size() == 200
    g.originate(np.random.default_rng(1), TargetParams())  # novel gene now, no new chromosome
    assert len(g.chromosomes) == 2 and g.size() > 200


def test_initial_chromosomes_reconstruct_end_to_end():
    res = simulate_nucleotide_genomes(
        _tree(seed=9), initial_chromosomes=3, inversion=0.02, duplication=0.005, loss=0.005,
        transposition=0.01, root_length=200, extension=0.9, seed=9)
    assert all(len(g.chromosomes) == 3 for g in res.leaf_genomes.values())  # no fission -> stays 3
    assert len(res.block_gene_trees()) == len(res.blocks)


# --------------------------------------------------------------------------- #
# Integration: driven through the forward simulator
# --------------------------------------------------------------------------- #
def test_zero_tier_rates_stay_single_chromosome():
    res = simulate_nucleotide_genomes(
        _tree(seed=1), inversion=0.01, loss=0.005, duplication=0.005,
        root_length=300, extension=0.9, seed=1)
    assert all(len(g.chromosomes) == 1 for g in res.leaf_genomes.values())


def test_fission_produces_multichromosome_leaves():
    res = simulate_nucleotide_genomes(
        _tree(seed=11), inversion=0.005, loss=0.004, duplication=0.004,
        fission=0.4, root_length=300, extension=0.9, seed=11)
    assert any(len(g.chromosomes) > 1 for g in res.leaf_genomes.values())
    # tier events are recorded in the karyotype log with parent/child chromosome ids
    fi = [r for r in res.event_log.chromosome_records if r.event is EventType.FISSION]
    assert fi and all(len(r.children) == 2 and r.parents[0] == r.children[0] for r in fi)


def test_all_tier_events_fire_and_reconstruct():
    fired = Counter()
    for seed in range(12):
        res = simulate_nucleotide_genomes(
            _tree(seed=seed), inversion=0.004, loss=0.003, duplication=0.003,
            fission=0.5, fusion=0.5, chromosome_origination=0.15, chromosome_loss=0.3,
            root_length=300, extension=0.9, seed=seed, retain_internal=True)
        for r in res.event_log.chromosome_records:
            fired[r.event] += 1
        # reconstruction must succeed and segment ids stay unique per genome
        assert len(res.block_gene_trees()) == len(res.blocks)
        for g in res.leaf_genomes.values():
            ids = [s.seg_id for s in g._iter_segments()]
            assert len(ids) == len(set(ids)) and len(g.chromosomes) >= 1
    # over a dozen seeds every one of the four tier events should have fired at least once
    assert set(fired) == {EventType.FISSION, EventType.FUSION,
                          EventType.CHROMOSOME_ORIGINATION, EventType.CHROMOSOME_LOSS}


def test_content_conserved_under_fission_and_fusion_only():
    # fission + fusion (no loss/deletion) must conserve each surviving family's copy number
    res = simulate_nucleotide_genomes(
        _tree(seed=4), inversion=0.006, fission=0.5, fusion=0.5,
        root_length=300, extension=0.9, seed=4)
    for g in res.leaf_genomes.values():
        # every genome still spells out one contiguous circular history per chromosome
        assert sum(len(c) for c in g.chromosomes.values()) == g.n_segments()


# --------------------------------------------------------------------------- #
# Transfer is the gene event that carries a family to another chromosome
# --------------------------------------------------------------------------- #
def _multi_chrom_recipient(lengths):
    """A genome whose chromosomes have exactly the given bp lengths (chromosome 0 seeded first)."""
    g = NucleotideGenome(IdManager(), root_length=lengths[0], extension=0.9)
    g.originate(np.random.default_rng(0), TargetParams())
    for L in lengths[1:]:
        chrom = g._new_root_chromosome()
        chrom.elements.append(g._new_segment(g.ids.new_family(), 0, L, 1))
    return g


def _donor_segment(rng, ids):
    """A one-segment additive transfer copy from a throwaway donor sharing ``ids`` (so its source
    id does not collide with the recipient's own seed sources)."""
    donor = NucleotideGenome(ids, root_length=30, extension=0.9)
    donor.originate(rng, TargetParams())
    return donor.extract_segment(Selection(genes=(), region=Region(0, 10, 6)), rng)


def test_transfer_choose_returns_chromosome_and_position():
    rng = np.random.default_rng(1)
    rec = _multi_chrom_recipient([100, 50, 20])
    seg = _donor_segment(rng, rec.ids)
    at = rec.choose_insertion_point(seg, rng)
    assert isinstance(at, tuple) and len(at) == 2            # (chrom_id, position)
    cid, pos = at
    assert cid in rec.chromosomes and isinstance(pos, int)


def test_transfer_lands_on_the_chosen_chromosome():
    rng = np.random.default_rng(2)
    rec = _multi_chrom_recipient([100, 50, 20])
    seg = _donor_segment(rng, rec.ids)
    donor_source = seg.genes[0].source
    at = rec.choose_insertion_point(seg, rng)
    cid, _pos = at
    rec.insert_segment(seg, at, rng)
    on = {c for c, chrom in rec.chromosomes.items()
          if any(s.source == donor_source for s in chrom.elements)}
    assert on == {cid}                                       # landed on exactly the chosen chromosome


def test_transfer_recipient_chromosome_is_uniform_not_size_weighted():
    """The recipient chromosome is chosen uniformly: a tiny chromosome receives transfers about as
    often as a large one (unlike the size-weighted choice that targets rearrangements)."""
    rng = np.random.default_rng(3)
    rec = _multi_chrom_recipient([10, 10, 10, 400])          # deliberately very uneven
    seg = _donor_segment(rng, rec.ids)
    hits = [rec.choose_insertion_point(seg, rng)[0] for _ in range(400)]
    assert set(hits) == set(rec.chromosomes)                 # every chromosome reachable
    cids = list(rec.chromosomes)
    # under uniform(4) each is expected ~100; the tiny and the huge chromosome are hit comparably,
    # which a bp-weighted choice (10/430 vs 400/430) never would.
    assert hits.count(cids[0]) > 40 and hits.count(cids[3]) > 40


def test_profiles_path_rejects_tier_events():
    tree = _tree(seed=2)
    for kw in (dict(fission=0.1), dict(fusion=0.1),
               dict(chromosome_origination=0.1), dict(chromosome_loss=0.1)):
        with pytest.raises(ValueError, match="chromosome-tier events"):
            simulate_nucleotide_genomes(tree, inversion=0.01, root_length=200,
                                        extension=0.9, output="profiles", seed=2, **kw)


# --------------------------------------------------------------------------- #
# Translocation: an arc moves to another chromosome of the same genome
# --------------------------------------------------------------------------- #
def _source_chromosomes(genome):
    """source -> set of chromosome ids it occupies in this genome."""
    on = {}
    for cid, chrom in genome.chromosomes.items():
        for s in chrom.elements:
            on.setdefault(s.source, set()).add(cid)
    return on


def test_translocation_moves_material_across_chromosomes():
    # two chromosomes, each seeded under its own source; a translocation carries a source's material
    # onto the other chromosome, so afterwards some source occupies >1 chromosome.
    res = simulate_nucleotide_genomes(
        _tree(seed=5), initial_chromosomes=2, inversion=0.005, translocation=0.02,
        root_length=250, extension=0.9, seed=5)
    crossed = sum(any(len(cids) > 1 for cids in _source_chromosomes(g).values())
                  for g in res.leaf_genomes.values())
    assert crossed >= 1                                    # material reached another chromosome
    # translocation never changes the chromosome count (it moves content, doesn't add/remove)
    assert all(len(g.chromosomes) == 2 for g in res.leaf_genomes.values())


def test_translocation_is_content_preserving_and_reconstructs():
    tree = _tree(seed=6)
    kw = dict(initial_chromosomes=2, root_length=250, extension=0.9, seed=6, retain_internal=True)
    withtl = simulate_nucleotide_genomes(tree, inversion=0.006, translocation=0.03, **kw)
    plain = simulate_nucleotide_genomes(tree, **kw)  # no content-changing events either way
    # translocation only relocates material, so every leaf's total size is unchanged
    assert ({l.name: g.size() for l, g in withtl.leaf_genomes.items()}
            == {l.name: g.size() for l, g in plain.leaf_genomes.items()})
    for g in withtl.leaf_genomes.values():                # segment ids stay unique (no re-mint)
        ids = [s.seg_id for s in g._iter_segments()]
        assert len(ids) == len(set(ids))
    assert len(withtl.block_gene_trees()) == len(withtl.blocks)   # gene trees unaffected


def test_translocation_needs_the_python_engine():
    with pytest.raises(ValueError, match="translocation"):
        simulate_nucleotide_genomes(_tree(seed=2), initial_chromosomes=2, translocation=0.02,
                                    root_length=200, extension=0.9, output="profiles", seed=2)
