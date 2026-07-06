"""Genes & intergenes in the nucleotide model.

Genes are user-supplied, non-overlapping intervals on the root chromosome that structural
events may never break — so each gene is exactly one block (one genealogy). Losses may
*pseudogenize* a gene (retain the sequence, flip it to intergene); replacement transfers are
*homologous* (the copy replaces the recipient's syntenic locus, found via flanking genes).
These tests pin the never-cut invariant, the pseudogenization state machine, the homologous
transfer, and that both gene and intergene trees are recovered.
"""

import numpy as np
import pytest

from zombi2 import BirthDeath, simulate_species_tree
from zombi2.events import EventType, Region, Selection, TargetParams
from zombi2.genome import IdManager
from zombi2.nucleotide_genome import NucleotideGenome, SegmentRegistry
from zombi2.nucleotide_sim import simulate_nucleotide_genomes


def _tree(n_tips=6, age=1.0, seed=0):
    return simulate_species_tree(BirthDeath(birth=1.0, death=0.25), n_tips=n_tips, age=age, seed=seed)


def _genic(root_length=100, genes=((20, 40, "geneA"),), **kw):
    reg = SegmentRegistry(pending_genes=[tuple(g) for g in genes])
    g = NucleotideGenome(IdManager(), root_length=root_length, extension=kw.pop("ext", 0.9),
                         registry=reg, **kw)
    g.originate(np.random.default_rng(0), TargetParams())
    return g


# --------------------------------------------------------------------------- #
# Seeding: the root chromosome is tiled into gene / intergene blocks
# --------------------------------------------------------------------------- #
def test_seed_tiles_into_gene_and_intergene_segments():
    g = _genic(100, [(20, 40, "geneA"), (60, 75, "geneB")])
    assert g.size() == 100
    tiling = [(s.src_start, s.src_end, s.gene_id) for s in g._segments]
    assert tiling == [(0, 20, None), (20, 40, "geneA"), (40, 60, None),
                      (60, 75, "geneB"), (75, 100, None)]
    assert all(s.is_gene for s in g._segments if s.gene_id is not None)


def test_seed_without_genes_is_one_segment():
    g = NucleotideGenome(IdManager(), root_length=100, extension=0.9)
    g.originate(np.random.default_rng(0), TargetParams())
    assert g.n_segments() == 1 and g._segments[0].gene_id is None


# --------------------------------------------------------------------------- #
# Breakpoint snapping: never cut inside a gene
# --------------------------------------------------------------------------- #
def test_snap_keeps_intergene_moves_out_of_genes():
    g = _genic(100, [(20, 40, "geneA")])
    # intergene interior: unchanged
    assert g._snap(10) == 10
    assert g._snap(50) == 50
    # inside the gene [20,40): snapped to a boundary (nearer / down / up)
    assert g._snap(25) in (20, 40)
    assert g._snap(25, -1) == 20
    assert g._snap(25, +1) == 40
    # boundaries themselves are legal
    assert g._snap(20) == 20 and g._snap(40) == 40


def test_draw_target_never_lands_inside_a_gene():
    g = _genic(200, [(30, 60, "geneA"), (100, 140, "geneB")])
    rng = np.random.default_rng(3)
    genes = [(30, 60), (100, 140)]
    for _ in range(500):
        sel = g.draw_target(EventType.LOSS, rng, TargetParams())
        s, ell = sel.region.start, sel.region.length
        L = g.size()
        e = (s + ell) % L
        for lo, hi in genes:          # neither endpoint strictly inside a gene
            assert not (lo < s < hi)
            assert not (lo < e < hi)


def test_subgene_event_is_promoted_to_whole_gene():
    g = _genic(100, [(20, 40, "geneA")])
    # a length-2 draw starting at 30 is entirely inside the gene -> promoted to [20,40)
    # (drive draw_target's snapping by hand via the private helpers it uses)
    s2 = g._snap(30, -1)
    e2 = g._snap((30 + 2) % 100, +1)
    assert (s2, e2) == (20, 40)


@pytest.mark.parametrize("seed", range(15))
def test_dense_genes_with_heavy_transposition(seed):
    """Regression: transposition into a densely-genic genome must snap its paste point out of
    genes without ever landing on the wrap boundary (a stale-index crash)."""
    tree = _tree(seed=seed)
    genes = [(2 + 20 * i, 20 + 20 * i, f"g{i}") for i in range(10)]  # ~90% genic, 2nt gaps
    res = simulate_nucleotide_genomes(
        tree, inversion=0.02, loss=0.01, duplication=0.008, transfer=0.008,
        transposition=0.02, origination=0.4, root_length=200, extension=0.9,
        gene_intervals=genes, pseudogenization=0.6, replacement=0.7, seed=seed)
    for a in res.blocks:                            # gene integrity still holds under stress
        if a.kind == "gene":
            assert a.gene_id is not None
    res.block_reconciliations()                     # trees build cleanly


@pytest.mark.parametrize("seed", range(30))
def test_genes_are_never_split_by_any_event(seed):
    """Full simulation: every gene block is a whole gene; no block straddles a gene boundary."""
    tree = _tree(seed=seed)
    genes = [(50, 90, "a"), (140, 200, "b"), (300, 340, "c"), (600, 660, "d")]
    res = simulate_nucleotide_genomes(
        tree, inversion=0.01, loss=0.006, duplication=0.004, transfer=0.004,
        transposition=0.003, origination=0.4, root_length=1000, extension=0.95,
        gene_intervals=genes, pseudogenization=0.5, replacement=0.6, seed=seed * 7 + 1)
    gene_iv = {(gi.source, gi.start, gi.end): gi.gene_id
               for gis in res.registry.genes.values() for gi in gis}
    for a in res.blocks:
        if a.kind == "gene":
            assert (a.source, a.start, a.end) in gene_iv
            assert a.gene_id == gene_iv[(a.source, a.start, a.end)]
        for (src, s, e) in gene_iv:                 # no block straddles a gene boundary
            if a.source == src and not (a.end <= s or a.start >= e):
                assert a.start >= s and a.end <= e
    # and per-leaf: each gene copy's ancestral positions form a contiguous run in the trace-back
    for leaf, genome in res.leaf_genomes.items():
        for seg in genome._segments:
            if seg.gene_id is None:
                continue
            span = seg.src_end - seg.src_start
            # a gene block is never internally cut: one segment covers the whole gene interval
            assert seg.gene_id in gene_iv.values()
            gi = next(gv for gv in
                      (g for gs in res.registry.genes.values() for g in gs)
                      if gv.gene_id == seg.gene_id)
            assert span == gi.length


# --------------------------------------------------------------------------- #
# Pseudogenization: retain sequence, flip state, one genealogy, lineage-specific
# --------------------------------------------------------------------------- #
def test_pseudogenization_retains_sequence_and_flips_state():
    g = _genic(100, [(20, 40, "geneA")], pseudogenization=1.0)
    before = g.size()
    # act on exactly the gene block [20,40)
    groups = g.apply(EventType.LOSS, Selection(genes=(), region=Region(0, 20, 20)),
                     np.random.default_rng(1), TargetParams())
    assert g.size() == before                       # sequence retained
    gene_segs = [s for s in g._segments if s.gene_id == "geneA"]
    assert gene_segs and all(not s.is_gene for s in gene_segs)   # demoted to pseudogene
    # logged as parent -> pseudogenized (a state change, not a terminal loss)
    assert groups and groups[0][0].role == "parent" and groups[0][1].role == "pseudogenized"


def test_pseudogenization_no_gene_in_arc_is_normal_deletion():
    g = _genic(100, [(20, 40, "geneA")], pseudogenization=1.0)
    before = g.size()
    # act on an intergene-only arc [50,70): no gene -> real deletion, length shrinks
    g.apply(EventType.LOSS, Selection(genes=(), region=Region(0, 50, 20)),
            np.random.default_rng(1), TargetParams())
    assert g.size() == before - 20


def test_pseudogenization_lineage_specific_and_in_tree():
    tree = _tree(n_tips=8, seed=11)
    genes = [(100, 200, "g1"), (400, 500, "g2"), (700, 780, "g3")]
    res = simulate_nucleotide_genomes(
        tree, inversion=0.003, loss=0.006, root_length=1000, extension=0.97,
        gene_intervals=genes, pseudogenization=0.6, seed=5)
    ps = res.pseudogenizations()
    assert ps, "expected at least one pseudogenization with these settings"
    # every flip is on a gene block, and it surfaces as a `G` node in that block's complete tree
    gene_trees = res.gene_trees()
    for block_id, gene_id, species, t, gid in ps:
        assert res._block_by_id[block_id].kind == "gene"
        complete, _extant = gene_trees[block_id]
        assert "|G" in complete
    # no new family/source was minted for a pseudogene (genealogy count == gene-block count logic):
    # every gene block's source is one of the original gene sources (or an originated gene)
    assert {a.gene_id for a in res.gene_blocks()} >= {p[1] for p in ps}


# --------------------------------------------------------------------------- #
# Homologous replacement transfer
# --------------------------------------------------------------------------- #
def _donor_recipient(genes, root_length=70, replacement=1.0):
    donor = _genic(root_length, genes, replacement=replacement)
    recipient, _ = donor.clone_reminting()
    return donor, recipient


def test_replacement_transfer_lands_at_homolog_and_logs_losses():
    donor, recipient = _donor_recipient([(10, 20, "L"), (30, 40, "M"), (50, 60, "R")])
    # extract the middle gene M [30,40); flanks resolve to genes L and R
    ts = donor.extract_segment(Selection(genes=(), region=Region(0, 30, 10)),
                               np.random.default_rng(2))
    assert ts.replacement and ts.left_flank == ("1", "L", 1) and ts.right_flank == ("1", "R", 1)
    at = recipient.choose_insertion_point(ts, np.random.default_rng(3))
    assert isinstance(at, tuple) and at[0] == "homolog"
    recipient.insert_segment(ts, at, np.random.default_rng(3))
    removed = recipient.pop_replaced_segments()
    # the recipient locus between L and R (intergene, gene M, intergene) was replaced by the copy
    assert any(s.gene_id == "M" for s in removed)
    assert {s.gene_id for s in recipient._segments if s.gene_id} == {"L", "M", "R"}


def test_replacement_requires_matching_flank_orientation():
    """A homologous locus must carry the flank genes in the SAME orientation. If a flank has been
    inverted in the recipient, synteny is broken and the transfer falls back to additive."""
    donor, recipient = _donor_recipient([(10, 20, "L"), (30, 40, "M"), (50, 60, "R")])
    for seg in recipient._segments:                  # invert the recipient's left flank L (+1 -> -1)
        if seg.gene_id == "L":
            seg.strand = -1
    ts = donor.extract_segment(Selection(genes=(), region=Region(0, 30, 10)),
                               np.random.default_rng(2))
    assert ts.left_flank == ("1", "L", 1)            # donor's flank is forward; recipient's is now -1
    at = recipient.choose_insertion_point(ts, np.random.default_rng(3))
    assert isinstance(at, int)                        # orientation mismatch -> additive, not a homolog
    recipient.insert_segment(ts, at, np.random.default_rng(3))
    assert recipient.pop_replaced_segments() == []    # nothing replaced


def test_replacement_falls_back_to_additive_without_homolog():
    donor, _ = _donor_recipient([(10, 20, "L"), (30, 40, "M"), (50, 60, "R")])
    # a recipient that lacks the flank genes entirely -> no homolog -> additive insertion
    other = _genic(40, [(5, 15, "X")], replacement=1.0)
    ts = donor.extract_segment(Selection(genes=(), region=Region(0, 30, 10)),
                               np.random.default_rng(2))
    at = other.choose_insertion_point(ts, np.random.default_rng(3))
    assert isinstance(at, int)                       # additive, not a homolog span
    before = other.size()
    other.insert_segment(ts, at, np.random.default_rng(3))
    assert other.pop_replaced_segments() == []        # nothing removed
    assert other.size() == before + 10                # copy added


def test_self_transfer_replacement_falls_back_to_additive():
    g = _genic(70, [(10, 20, "L"), (30, 40, "M"), (50, 60, "R")], replacement=1.0)
    ts = g.extract_segment(Selection(genes=(), region=Region(0, 30, 10)),
                           np.random.default_rng(2))
    # recipient is the donor itself: its continuation segments are present -> skip homology
    at = g.choose_insertion_point(ts, np.random.default_rng(3))
    assert isinstance(at, int)


def test_replacement_transfers_run_and_reconcile(tmp_path):
    tree = _tree(n_tips=8, seed=21)
    genes = [(100, 160, "g1"), (300, 360, "g2"), (500, 560, "g3"), (700, 760, "g4")]
    res = simulate_nucleotide_genomes(
        tree, transfer=0.01, loss=0.001, root_length=1000, extension=0.97,
        gene_intervals=genes, replacement=1.0, seed=4)
    # replacement transfers log recipient losses; the full reconciliation still builds cleanly
    summary = res.write_reconciliations(tmp_path)
    assert summary["n_blocks"] >= 1
    # every reconciliation-event species branch is a real node in the tree
    nodes = {n.name for n in tree.nodes_preorder()}
    for _block_id, rec in res.block_reconciliations().items():
        for e in rec.events:
            assert e.species in nodes


# --------------------------------------------------------------------------- #
# Gene + intergene trees recovered, classification
# --------------------------------------------------------------------------- #
def test_gene_and_intergene_trees_partition_all_blocks():
    tree = _tree(n_tips=6, seed=2)
    genes = [(100, 180, "a"), (300, 360, "b"), (500, 620, "c"), (750, 800, "d")]
    res = simulate_nucleotide_genomes(
        tree, inversion=0.003, loss=0.002, duplication=0.001, transfer=0.001,
        root_length=1000, extension=0.97, gene_intervals=genes, seed=7)
    gene_ids = {a.block_id for a in res.gene_blocks()}
    inter_ids = {a.block_id for a in res.intergene_blocks()}
    assert gene_ids.isdisjoint(inter_ids)
    assert gene_ids | inter_ids == {a.block_id for a in res.blocks}
    assert set(res.gene_trees()) == gene_ids
    assert set(res.intergene_trees()) == inter_ids
    assert res.gene_trees() and res.intergene_trees()
    # each gene block is a full user gene interval (or an originated gene)
    for a in res.gene_blocks():
        assert a.gene_id is not None


def _n_leaves(newick):
    return newick.count(",") + 1


@pytest.mark.parametrize("seed", range(12))
def test_genic_reconciliation_invariant(seed):
    """A block's extant tree has one leaf per surviving copy — gene, intergene, and pseudogene
    tips alike — even with pseudogenization (unary `G` nodes) and homologous replacement."""
    tree = _tree(n_tips=6, seed=seed)
    genes = [(80, 140, "a"), (250, 300, "b"), (450, 540, "c"), (700, 760, "d")]
    res = simulate_nucleotide_genomes(
        tree, inversion=0.004, loss=0.004, duplication=0.004, transfer=0.004,
        root_length=900, extension=0.95, gene_intervals=genes,
        pseudogenization=0.5, replacement=0.5, seed=seed * 5 + 2)
    ids, _species, M = res.profile_matrix()
    rowsum = {aid: int(M[i].sum()) for i, aid in enumerate(ids)}
    for aid, (_complete, extant) in res.block_gene_trees().items():
        if extant is not None:                     # some blocks may survive in no extant leaf
            assert _n_leaves(extant) == rowsum[aid]


def test_reproducible_with_genes():
    tree = _tree(seed=1)
    kw = dict(inversion=0.004, loss=0.003, duplication=0.002, transfer=0.002,
              root_length=800, extension=0.96, gene_intervals=[(50, 120, "g")],
              pseudogenization=0.4, replacement=0.5, seed=99)
    a = simulate_nucleotide_genomes(tree, **kw)
    b = simulate_nucleotide_genomes(tree, **kw)
    ca = {n.name: a.trace_back(n) for n in a.leaf_genomes}
    cb = {n.name: b.trace_back(n) for n in b.leaf_genomes}
    assert ca == cb
    assert len(a.event_log) == len(b.event_log)


# --------------------------------------------------------------------------- #
# Origination mints a new gene; validation; Rust guard
# --------------------------------------------------------------------------- #
def test_origination_creates_a_new_gene():
    tree = _tree(seed=8)
    res = simulate_nucleotide_genomes(
        tree, origination=1.5, root_length=500, extension=0.95,
        gene_intervals=[(50, 100, "seed")], seed=3)
    ids = {gi.gene_id for gis in res.registry.genes.values() for gi in gis}
    assert "seed" in ids and len(ids) > 1            # novel originated genes appeared
    # every originated gene is a whole gene block (never split)
    assert all(a.gene_id is not None for a in res.gene_blocks())


def test_profiles_output_rejects_genes():
    tree = _tree(seed=0)
    with pytest.raises(ValueError, match="Python engine"):
        simulate_nucleotide_genomes(tree, inversion=0.001, root_length=100,
                                    gene_intervals=[(10, 20)], output="profiles")
    with pytest.raises(ValueError, match="pseudogenization"):
        simulate_nucleotide_genomes(tree, inversion=0.001, root_length=100,
                                    pseudogenization=0.5, output="profiles")


def test_invalid_gene_intervals_raise():
    tree = _tree(seed=0)
    with pytest.raises(ValueError, match="overlap"):
        simulate_nucleotide_genomes(tree, root_length=100,
                                    gene_intervals=[(10, 40), (30, 50)])
    with pytest.raises(ValueError, match="root_length"):
        simulate_nucleotide_genomes(tree, root_length=100, gene_intervals=[(80, 120)])
