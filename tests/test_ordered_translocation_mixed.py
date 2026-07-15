"""Ordered-genome parity with the nucleotide model: translocation + mixed chromosome topology.

The nucleotide chromosome tier already ships two features the ordered model lacked — moving a
segment *between* chromosomes (**translocation**), and seeding a genome whose chromosomes are a
**mix** of circular and linear replicons. This module closes that gap and pins it down.

Two invariants carry the correctness argument and are exercised adversarially here:

* **Translocation is lineage-neutral.** It only changes *which* chromosome holds a gene run; gene
  ids are never re-minted, so the multiset of gids is conserved and the gene trees never see it —
  exactly like transposition, but cross-chromosome. It is a pure move: total gene count conserved,
  the source sheds precisely the run, some other chromosome gains it.

* **Fusion is same-topology only.** Fusing a linear replicon into a circular one is topologically
  ill-defined, so a chromosome fuses only with a same-topology partner (and the event is a no-op —
  ``fusion`` returns ``None``, the caller skips it — when there is none). For a topology-homogeneous
  genome every other chromosome qualifies, so the draws stay byte-identical to the pre-mixed engine
  (guarded by test_chromosome_fingerprints.py).
"""

import numpy as np
import pytest

from zombi2 import (
    BirthDeath,
    OrderedGenome,
    SharedRates,
    simulate_genomes,
    simulate_species_tree,
)
from zombi2.cli import main
from zombi2.genomes.genome import IdManager, OrderedGene
from zombi2.genomes.events import EventType, TargetParams


def _seed_genes(genome, ids, per_chromosome=5):
    """Drop ``per_chromosome`` throwaway genes on every chromosome (families are irrelevant here)."""
    for chrom in genome.chromosomes.values():
        for _ in range(per_chromosome):
            chrom.genes.append(OrderedGene(ids.new_gene(), "f", 1))


def _all_events_rates(**over):
    kw = dict(duplication=0.3, transfer=0.25, loss=0.3, origination=0.15,
              inversion=0.2, transposition=0.2, translocation=0.4,
              chromosome_origination=0.08, chromosome_loss=0.05, fission=0.1, fusion=0.1)
    kw.update(over)
    return SharedRates(**kw)


# --- Part A: translocation is a lineage-neutral pure move -----------------------------------

def test_translocation_conserves_gids_and_gene_counts():
    """400 translocations on a 3-chromosome genome: gids conserved (lineage-neutral), total gene
    count conserved (a move, not a birth/death), and the source chromosome sheds exactly the run."""
    rng = np.random.default_rng(0)
    ids = IdManager()
    g = OrderedGenome(ids, extension=0.7, n_chromosomes=3, circular=True)
    for _ in range(30):
        g.originate(rng, TargetParams())

    def snapshot():
        gids = sorted(gene.gid for c in g.chromosomes.values() for gene in c.genes)
        per = {cid: len(c.genes) for cid, c in g.chromosomes.items()}
        return gids, per

    fires = 0
    for _ in range(400):
        before_gids, before_per = snapshot()
        sel = g.draw_target(EventType.TRANSLOCATION, rng, TargetParams(extension=0.7))
        src, moved = sel.region.chromosome, len(sel.genes)
        groups = g.apply(EventType.TRANSLOCATION, sel, rng, TargetParams())
        after_gids, after_per = snapshot()
        assert after_gids == before_gids, "translocation re-minted a gid — not lineage-neutral"
        if groups:  # >= 2 chromosomes present, so it fired
            fires += 1
            assert sum(after_per.values()) == sum(before_per.values()), "gene count changed"
            assert after_per[src] == before_per[src] - moved, "source did not shed the run"
            assert all(op.role == "translocated" for group in groups for op in group)
    assert fires == 400, fires


def test_translocation_needs_a_second_chromosome():
    """On a single-chromosome genome there is nowhere to move to: ``apply`` is a no-op (returns [])
    and nothing is lost — the sampler's n_chrom>=2 gate is mirrored defensively in the genome."""
    rng = np.random.default_rng(1)
    ids = IdManager()
    g = OrderedGenome(ids, n_chromosomes=1, circular=True)
    for _ in range(10):
        g.originate(rng, TargetParams())
    before = sorted(gene.gid for gene in next(iter(g.chromosomes.values())).genes)
    sel = g.draw_target(EventType.TRANSLOCATION, rng, TargetParams(extension=0.5))
    assert g.apply(EventType.TRANSLOCATION, sel, rng, TargetParams()) == []
    after = sorted(gene.gid for gene in next(iter(g.chromosomes.values())).genes)
    assert after == before  # genome untouched


def test_translocation_is_a_supported_ordered_event():
    g = OrderedGenome(IdManager(), n_chromosomes=2)
    assert EventType.TRANSLOCATION in g.supported_events()


@pytest.mark.parametrize("seed", range(6))
def test_translocation_fires_and_reconstructs(seed):
    """A high-rate translocation run simulates end to end: gene trees and reconciliations both
    reconstruct (translocation records are lineage-neutral, so the genealogy is well-formed)."""
    rates = SharedRates(duplication=0.3, transfer=0.2, loss=0.3, origination=0.1,
                        inversion=0.2, transposition=0.2, translocation=0.6, fission=0.05)
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=seed)
    g = simulate_genomes(tree, rates, initial_families=15, seed=1000 + seed,
                         genome_factory=lambda ids: OrderedGenome(ids, extension=0.7,
                                                                  n_chromosomes=2, circular=True))
    n_x = sum(1 for r in g.event_log.records if r.event is EventType.TRANSLOCATION)
    assert n_x > 0, "translocation never fired at rate 0.6 with two chromosomes"
    assert g.gene_trees()            # reconstructs
    assert g.reconciliations()       # reconstructs
    # every translocation record is lineage-neutral (role 'translocated', gene ids intact)
    for r in g.event_log.records:
        if r.event is EventType.TRANSLOCATION:
            assert all(op.role == "translocated" for op in r.genes)


# --- Part B: mixed circular/linear topology -------------------------------------------------

def test_mixed_topology_seeding():
    g = OrderedGenome(IdManager(), extension=0.6, circular=[True, False, True, False])
    assert [c.circular for c in g.chromosomes.values()] == [True, False, True, False]


def test_circular_sequence_sets_the_chromosome_count():
    """A per-chromosome ``circular`` sequence sets the count on its own; passing a conflicting
    ``n_chromosomes`` is rejected, but the default n_chromosomes=1 is accepted."""
    g = OrderedGenome(IdManager(), circular=[True, False, False])
    assert len(g.chromosomes) == 3
    with pytest.raises(ValueError, match="n_chromosomes must equal len"):
        OrderedGenome(IdManager(), circular=[True, False], n_chromosomes=3)


def test_mixed_topology_survives_clone():
    ids = IdManager()
    g = OrderedGenome(ids, extension=0.6, circular=[True, False, True, False])
    _seed_genes(g, ids, per_chromosome=3)
    child, _ = g.clone_reminting()
    assert [c.circular for c in child.chromosomes.values()] == [True, False, True, False]


def test_fusion_is_same_topology_only():
    """2 circular + 2 linear, fused 2000 times: fusion never crosses topology, and it returns a
    same-topology pair (never ``None`` here, because each topology has a partner)."""
    rng = np.random.default_rng(0)
    crossed = none = 0
    for _ in range(2000):
        ids = IdManager()
        g = OrderedGenome(ids, extension=0.6, circular=[True, True, False, False])
        _seed_genes(g, ids, per_chromosome=3)
        topo = {cid: c.circular for cid, c in g.chromosomes.items()}
        res = g.fusion(rng, TargetParams())
        if res is None:
            none += 1
        else:
            keep, absorb = res
            if topo[keep] != topo[absorb]:
                crossed += 1
    assert crossed == 0, f"fusion crossed topology {crossed} times"
    assert none == 0, "every topology had a partner, so fusion should always have fired"


def test_fusion_returns_none_when_no_same_topology_partner():
    """One circular + one linear: neither has a same-topology partner, so fusion is always a
    no-op (``None``) — which the simulation loop skips."""
    rng = np.random.default_rng(2)
    for _ in range(200):
        ids = IdManager()
        g = OrderedGenome(ids, extension=0.6, circular=[True, False])
        _seed_genes(g, ids, per_chromosome=4)
        assert g.fusion(rng, TargetParams()) is None
        assert len(g.chromosomes) == 2  # untouched


@pytest.mark.parametrize("seed", range(8))
def test_mixed_topology_all_events_reconstruct(seed):
    """The full integration: a root of 2 circular + 2 linear replicons under *every* event
    (including translocation and the whole chromosome tier) simulates and reconstructs, and every
    chromosome on every leaf keeps a well-formed boolean topology."""
    rates = _all_events_rates()
    tree = simulate_species_tree(BirthDeath(1.0, 0.25), n_tips=10, age=3.5, seed=seed)
    g = simulate_genomes(tree, rates, initial_families=18, seed=5000 + seed,
                         genome_factory=lambda ids: OrderedGenome(
                             ids, extension=0.6, circular=[True, False, True, False]))
    assert g.gene_trees()
    assert g.reconciliations()
    for genome in g.leaf_genomes.values():
        for chrom in genome.chromosomes.values():
            assert isinstance(chrom.circular, bool)


# --- CLI wiring -----------------------------------------------------------------------------

def _species(tmp_path):
    sp = tmp_path / "sp"
    main(["species", "--birth", "1", "--death", "0.3", "--tips", "8", "--age", "3",
          "--seed", "1", "-o", str(sp)])
    return str(sp / "species_tree.nwk")


def test_cli_ordered_translocation_runs(tmp_path):
    """`--translocation` is accepted for the ordered model with >= 2 chromosomes and produces a
    layout (gene_order.tsv is auto-surfaced for a multi-chromosome run)."""
    tree = _species(tmp_path)
    out = tmp_path / "gen"
    rc = main(["genomes", "--tree", tree, "--genome-model", "ordered", "--n-chromosomes", "2",
               "--dup", "0.2", "--loss", "0.2", "--orig", "0.3", "--translocation", "0.5",
               "--initial-families", "18", "--seed", "3", "-o", str(out)])
    assert rc == 0
    assert (out / "gene_order.tsv").exists()


def test_cli_translocation_rejected_without_chromosomes(tmp_path):
    """The unordered model has no chromosomes, so `--translocation` there is a usage error."""
    tree = _species(tmp_path)
    with pytest.raises(SystemExit):
        main(["genomes", "--tree", tree, "--genome-model", "unordered", "--translocation", "0.3",
              "--seed", "1", "-o", str(tmp_path / "x")])
