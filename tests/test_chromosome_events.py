"""Chromosome-tier events for the OrderedGenome model (Stage 4a): origination (a de-novo
plasmid), whole-chromosome loss, fission (1->2) and fusion (2->1).

These act one tier above the gene events. Fission and fusion only move genes between
chromosomes, so gene ids (lineages) are preserved and gene-tree reconstruction is unaffected;
only chromosome loss ends gene lineages. The four operations are tested here in isolation; the
sampler wiring + gated rates (byte-identical when off) come in the next sub-step.
"""

import numpy as np

from zombi2 import OrderedGenome
from zombi2.genomes.genome import IdManager, OrderedGene


def _seeded(n_chromosomes, circular, per_chrom=4):
    """A genome with ``per_chrom`` distinct-family genes on each of ``n_chromosomes`` chromosomes."""
    ids = IdManager()
    genome = OrderedGenome(ids, n_chromosomes=n_chromosomes, circular=circular)
    for c, chrom in enumerate(genome.chromosomes.values()):
        for k in range(per_chrom):
            chrom.genes.append(OrderedGene(ids.new_gene(), f"c{c}f{k}", 1))
    return genome, ids


def _gids(genome):
    return sorted(g.gid for chrom in genome.chromosomes.values() for g in chrom.genes)


# --- origination (a de-novo plasmid) --------------------------------------------------

def test_originate_chromosome_adds_empty_replicon():
    genome, ids = _seeded(2, circular=True)
    genes_before = ids._gene
    before = set(genome.chromosomes)
    cid = genome.originate_chromosome(np.random.default_rng(0), None)
    assert cid not in before                           # a fresh id
    assert len(genome.chromosomes) == 3
    assert genome.chromosomes[cid].genes == []         # an empty replicon
    assert genome.chromosomes[cid].circular is True    # the genome's topology
    assert ids._gene == genes_before                   # a chromosome id never touches the gene counter


def test_originate_chromosome_matches_genome_topology():
    genome, _ = _seeded(1, circular=False)
    cid = genome.originate_chromosome(np.random.default_rng(0), None)
    assert genome.chromosomes[cid].circular is False


# --- whole-chromosome loss ------------------------------------------------------------

def test_lose_chromosome_removes_it_and_reports_gene_losses():
    genome, _ = _seeded(3, circular=True, per_chrom=4)
    all_before = set(_gids(genome))
    cid, groups = genome.lose_chromosome(np.random.default_rng(1))
    assert cid not in genome.chromosomes               # the chromosome is gone
    assert len(genome.chromosomes) == 2
    lost = {op.gid for grp in groups for op in grp}
    assert len(lost) == 4                              # one LOSS group per gene it held
    assert all(grp[0].role == "lost" for grp in groups)
    assert set(_gids(genome)) == all_before - lost     # exactly those genes are gone; the rest survive


# --- fission (1 -> 2) -----------------------------------------------------------------

def test_fission_linear_one_breakpoint_conserves_content():
    genome, _ = _seeded(1, circular=False, per_chrom=8)
    (only,) = genome.chromosomes.values()
    original = list(only.genes)                        # capture the order before the cut
    src, new = genome.fission(np.random.default_rng(2), None)
    assert len(genome.chromosomes) == 2
    assert genome.chromosomes[new].circular is False   # a linear daughter
    # one cut: prefix (src) then suffix (new) reproduce the original order exactly, ids unchanged
    assert genome.chromosomes[src].genes + genome.chromosomes[new].genes == original


def test_fission_circular_two_breakpoints_conserves_content():
    genome, _ = _seeded(1, circular=True, per_chrom=10)
    before = _gids(genome)
    src, new = genome.fission(np.random.default_rng(3), None)
    assert len(genome.chromosomes) == 2
    assert genome.chromosomes[src].circular and genome.chromosomes[new].circular
    assert _gids(genome) == before                     # gene-id multiset conserved (no re-minting)


def test_fission_targets_a_nonempty_chromosome():
    """Fission is size-weighted, so it never picks an empty chromosome (which could not be split)."""
    genome, ids = _seeded(1, circular=True, per_chrom=6)
    empty = genome.originate_chromosome(np.random.default_rng(0), None)  # an empty plasmid
    for seed in range(20):
        src, _ = genome.fission(np.random.default_rng(seed), None)
        assert src != empty


# --- fusion (2 -> 1) ------------------------------------------------------------------

def test_fusion_merges_two_into_one_conserves_content():
    genome, _ = _seeded(3, circular=True, per_chrom=4)
    before = _gids(genome)
    keep, absorbed = genome.fusion(np.random.default_rng(4), None)
    assert absorbed not in genome.chromosomes and keep in genome.chromosomes
    assert len(genome.chromosomes) == 2
    assert _gids(genome) == before                     # every gene conserved, ids stable


# --- invariant over a random mix ------------------------------------------------------

def test_content_conserved_over_a_mix_of_chromosome_events():
    """Origination / fission / fusion never change the gene multiset; only loss removes genes. Over
    a random mix, the surviving gene ids stay equal to the originals minus everything ever lost."""
    genome, _ = _seeded(2, circular=True, per_chrom=6)
    rng = np.random.default_rng(7)
    alive = set(_gids(genome))
    for _ in range(200):
        r = rng.random()
        if r < 0.3:
            genome.originate_chromosome(rng, None)
        elif r < 0.6 and genome.size() > 0:
            genome.fission(rng, None)
        elif r < 0.8 and len(genome.chromosomes) >= 2:
            genome.fusion(rng, None)
        elif len(genome.chromosomes) >= 2:
            _, groups = genome.lose_chromosome(rng)
            alive -= {op.gid for grp in groups for op in grp}
        assert set(_gids(genome)) == alive             # the invariant holds after every event
    assert len(genome.chromosomes) >= 1                # a genome always keeps at least one chromosome
