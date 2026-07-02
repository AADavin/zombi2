"""Nucleotide-level circular genome, M1 (inversions only).

The workhorse is a dead-simple O(L) **array oracle** that applies each inversion by
literal circular slice-reverse over one cell per nucleotide. It is obviously correct and
hopelessly slow; running it against the efficient segment structure on the *same*
(s, length) stream pins down every circular-geometry bug. On top of that we assert the
strong invariants inversion gives us (content conservation, bijection, involution) and
that the trace-back / atom decomposition is self-consistent.
"""

import numpy as np
import pytest

from zombi2 import BirthDeath, simulate_species_tree
from zombi2.events import EventType, Region, Selection, TargetParams
from zombi2.genome import IdManager
from zombi2.nucleotide_genome import NucleotideGenome, SegmentRegistry
from zombi2.nucleotide_sim import simulate_nucleotide_genomes


# --------------------------------------------------------------------------- #
# Independent brute-force oracle: one cell per nucleotide, fixed origin.
# --------------------------------------------------------------------------- #
class ArrayGenome:
    """One cell per nucleotide; same rotate-on-wrap origin convention as the genome."""

    def __init__(self, source, length):
        self.cells = [(source, i, 1) for i in range(length)]

    def invert(self, s, ell):
        L = len(self.cells)
        ell = max(1, min(ell, L))
        s %= L
        if s + ell > L:  # wrapping: rotate the origin to s, then reverse the prefix
            self.cells = self.cells[s:] + self.cells[:s]
            s = 0
        seg = [(src, p, -st) for (src, p, st) in self.cells[s:s + ell][::-1]]
        self.cells[s:s + ell] = seg

    def delete(self, s, ell):
        L = len(self.cells)
        ell = min(ell, L)
        s %= L
        if ell == L:
            self.cells = []
        elif s + ell > L:                       # wrapping: keep the middle [e, s)
            self.cells = self.cells[(s + ell) % L:s]
        else:                                   # non-wrapping: drop [s, s+ell)
            del self.cells[s:s + ell]

    def duplicate(self, s, ell):
        L = len(self.cells)
        ell = max(1, min(ell, L))
        s %= L
        if s + ell > L:                         # wrapping: rotate origin to s (as the genome)
            self.cells = self.cells[s:] + self.cells[:s]
            s = 0
        self.cells[s + ell:s + ell] = self.cells[s:s + ell]  # tandem copy right after

    def transpose(self, s, ell, dest):
        L = len(self.cells)
        if L <= 1:
            return
        ell = max(1, min(ell, L - 1))
        s %= L
        if s + ell > L:                         # wrapping: rotate origin to s (as the genome)
            self.cells = self.cells[s:] + self.cells[:s]
            s = 0
        block = self.cells[s:s + ell]
        del self.cells[s:s + ell]
        self.cells[dest % (len(self.cells) + 1):dest % (len(self.cells) + 1)] = block


def _fresh(length, ext=0.99):
    g = NucleotideGenome(IdManager(), root_length=length, extension=ext)
    g.originate(np.random.default_rng(0), TargetParams())  # seed the root chromosome
    return g


# --------------------------------------------------------------------------- #
# Seeding
# --------------------------------------------------------------------------- #
def test_seed_is_single_forward_segment():
    g = _fresh(1000)
    assert g.size() == 1000
    assert g.n_segments() == 1
    assert g.to_cells() == [("1", i, 1) for i in range(1000)]


# --------------------------------------------------------------------------- #
# Geometry vs the oracle (the workhorse)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("seed", range(25))
def test_inversion_matches_oracle_random(seed):
    rng = np.random.default_rng(seed)
    L = int(rng.integers(20, 200))
    g = _fresh(L, ext=0.9)
    o = ArrayGenome("1", L)
    for _ in range(80):
        s = int(rng.integers(L))
        ell = int(rng.integers(1, L + 1))
        g._apply_inversion(s, ell)
        o.invert(s, ell)
        assert g.to_cells() == o.cells       # identical content & layout
        assert g.size() == L                  # inversion preserves length


def test_wrapping_and_edge_cases_match_oracle():
    L = 12
    g = _fresh(L, ext=0.7)
    o = ArrayGenome("1", L)
    for s, ell in [(8, 6), (0, 12), (11, 3), (5, 12), (3, 9), (7, 5),
                   (0, 1), (11, 1), (10, 4), (6, 12)]:
        g._apply_inversion(s, ell)  # wrapping arcs, whole-genome, length-1
        o.invert(s, ell)
        assert g.to_cells() == o.cells


def test_whole_genome_inversion_reverses_all():
    L = 8
    g = _fresh(L)
    g._apply_inversion(0, L)
    assert g.to_cells() == [("1", p, -1) for p in range(L - 1, -1, -1)]


def test_length_one_inversion_flips_one_strand():
    g = _fresh(6)
    before = g.to_cells()
    g._apply_inversion(3, 1)
    after = g.to_cells()
    assert after[3] == (before[3][0], before[3][1], -before[3][2])
    assert after[:3] == before[:3] and after[4:] == before[4:]


def test_inversion_involution():
    g = _fresh(50, ext=0.9)
    before = g.to_cells()
    g._apply_inversion(12, 20)
    g._apply_inversion(12, 20)  # same arc twice == identity on content
    assert g.to_cells() == before


def test_split_preserves_content():
    g = _fresh(30)
    before = g.to_cells()
    g._split_at(7)
    g._split_at(19)
    g._split_at(7)  # idempotent
    assert g.n_segments() == 3
    assert g.to_cells() == before  # a breakpoint reorders nothing


# --------------------------------------------------------------------------- #
# Tree-level: content conservation, bijection, reproducibility
# --------------------------------------------------------------------------- #
def test_content_conserved_and_bijective_at_leaves():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=1)
    res = simulate_nucleotide_genomes(tree, inversion=0.02, root_length=200,
                                      extension=0.9, seed=7)
    assert any(r.event is EventType.INVERSION for r in res.event_log)
    root = {("1", i) for i in range(200)}
    for genome in res.leaf_genomes.values():
        cells = genome.to_cells()
        origins = [(src, p) for (src, p, _st) in cells]
        assert len(cells) == 200
        assert set(origins) == root            # content conserved
        assert len(origins) == len(set(origins))  # bijection: each origin once


def test_reproducible():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=2)
    a = simulate_nucleotide_genomes(tree, inversion=0.02, root_length=200, seed=9)
    b = simulate_nucleotide_genomes(tree, inversion=0.02, root_length=200, seed=9)
    assert len(a.event_log) == len(b.event_log)
    for leaf, ga in a.leaf_genomes.items():
        assert ga.to_cells() == b.leaf_genomes[leaf].to_cells()


# --------------------------------------------------------------------------- #
# Trace-back / atoms
# --------------------------------------------------------------------------- #
def test_atoms_tile_the_root():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=6, age=3.0, seed=3)
    res = simulate_nucleotide_genomes(tree, inversion=0.03, root_length=150, seed=4)
    atoms = sorted((a for a in res.atoms if a.source == "1"), key=lambda a: a.start)
    assert atoms[0].start == 0 and atoms[-1].end == 150
    for x, y in zip(atoms, atoms[1:]):
        assert x.end == y.start                # contiguous, no gaps/overlaps
    n_inv = sum(1 for r in res.event_log if r.event is EventType.INVERSION)
    assert len(atoms) <= 2 * n_inv + 1         # each inversion adds <= 2 breakpoints


def test_mosaic_reassembles_each_leaf():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=6, age=3.0, seed=5)
    res = simulate_nucleotide_genomes(tree, inversion=0.03, root_length=150,
                                      extension=0.9, seed=6)
    amap = {a.atom_id: a for a in res.atoms}
    for leaf, genome in res.leaf_genomes.items():
        mosaic = res.leaf_mosaic(leaf)
        # every atom appears exactly once per leaf (inversion loses nothing)
        assert sorted(aid for aid, _ in mosaic) == sorted(a.atom_id for a in res.atoms)
        cells = []
        for aid, strand in mosaic:
            a = amap[aid]
            if strand == 1:
                cells.extend((a.source, p, 1) for p in range(a.start, a.end))
            else:
                cells.extend((a.source, p, -1) for p in range(a.end - 1, a.start - 1, -1))
        assert cells == genome.to_cells()      # mosaic reconstructs the leaf exactly


def test_atom_histories_track_inversions():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=6, age=3.0, seed=8)
    res = simulate_nucleotide_genomes(tree, inversion=0.05, root_length=120,
                                      extension=0.9, seed=8)
    histories = res.atom_histories()
    branches = {n.name for n in tree.nodes_preorder()}
    n_inv = sum(1 for r in res.event_log if r.event is EventType.INVERSION)
    if n_inv:
        assert any(histories[a.atom_id] for a in res.atoms)  # something recorded
    for events in histories.values():
        for branch, _t in events:
            assert branch in branches           # every entry is a real branch


# --------------------------------------------------------------------------- #
# M2: deletion — content shrinks (subset, still bijective), vs the oracle
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("seed", range(25))
def test_inversion_and_deletion_match_oracle(seed):
    rng = np.random.default_rng(seed)
    L0 = int(rng.integers(40, 200))
    g = _fresh(L0, ext=0.9)
    o = ArrayGenome("1", L0)
    for _ in range(120):
        L = g.size()
        if L <= 1:
            break
        s = int(rng.integers(L))
        ell = int(rng.integers(1, L + 1))
        if rng.random() < 0.5:
            g._apply_inversion(s, ell)
            o.invert(s, ell)
        else:
            ell = min(ell, L - 1)              # keep >= 1 nt so the run continues
            g._apply_loss(s, ell)
            o.delete(s, ell)
        assert g.to_cells() == o.cells
        assert g.size() == len(o.cells)


def test_wrapping_deletion_matches_oracle():
    g = _fresh(20, ext=0.7)
    o = ArrayGenome("1", 20)
    for s, ell in [(18, 5), (0, 3), (10, 4), (6, 5)]:  # includes origin-crossing deletes
        g._apply_loss(s, ell)
        o.delete(s, ell)
        assert g.to_cells() == o.cells
        assert g.size() == len(o.cells)


def test_deletion_removes_content_but_stays_bijective():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=11)
    res = simulate_nucleotide_genomes(tree, inversion=0.01, loss=0.02, root_length=300,
                                      extension=0.95, seed=11)
    assert any(r.event is EventType.LOSS for r in res.event_log)
    root = {("1", i) for i in range(300)}
    total = 0
    for genome in res.leaf_genomes.values():
        origins = [(src, p) for (src, p, _st) in genome.to_cells()]
        assert set(origins) <= root                 # only ancestral material, no novelty
        assert len(origins) == len(set(origins))    # no duplication -> still bijective
        total += len(origins)
    # something was actually deleted somewhere (expected leaf shorter than the root)
    assert any(g.size() < 300 for g in res.leaf_genomes.values())


def test_atoms_cover_exactly_the_surviving_positions():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=6, age=3.0, seed=12)
    res = simulate_nucleotide_genomes(tree, inversion=0.02, loss=0.03, root_length=200,
                                      extension=0.95, seed=12)
    atoms = sorted((a for a in res.atoms if a.source == "1"), key=lambda a: a.start)
    for x, y in zip(atoms, atoms[1:]):
        assert x.end <= y.start                       # disjoint (gaps allowed)
    atom_positions = {p for a in atoms for p in range(a.start, a.end)}
    surviving = set()
    for genome in res.leaf_genomes.values():
        surviving.update(p for (_src, p, _st) in genome.to_cells())
    assert atom_positions == surviving                # atoms == union of survivors


def test_profile_matrix_matches_leaf_coverage():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=7, age=3.0, seed=13)
    res = simulate_nucleotide_genomes(tree, inversion=0.02, loss=0.03, root_length=200,
                                      extension=0.95, seed=13)
    atom_ids, species, matrix = res.profile_matrix()
    assert matrix.shape == (len(res.atoms), len(res.leaf_genomes))
    assert set(matrix.flatten()) <= {0, 1}            # loss only -> presence/absence
    # not every atom is universal, and none is everywhere-absent (that isn't an atom)
    assert matrix.min() == 0 and matrix.max() == 1
    assert matrix.sum(axis=1).min() >= 1              # each atom present in >= 1 leaf


def test_loss_reproducible():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=14)
    a = simulate_nucleotide_genomes(tree, inversion=0.02, loss=0.03, root_length=200, seed=15)
    b = simulate_nucleotide_genomes(tree, inversion=0.02, loss=0.03, root_length=200, seed=15)
    assert len(a.event_log) == len(b.event_log)
    for leaf, ga in a.leaf_genomes.items():
        assert ga.to_cells() == b.leaf_genomes[leaf].to_cells()


# --------------------------------------------------------------------------- #
# M3: duplication — paralogs (value-identical copies), content grows, vs oracle
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("seed", range(25))
def test_all_events_match_oracle(seed):
    rng = np.random.default_rng(seed)
    L0 = int(rng.integers(30, 120))
    g = _fresh(L0, ext=0.9)
    o = ArrayGenome("1", L0)
    for _ in range(150):
        L = g.size()
        if L <= 1:
            break
        s = int(rng.integers(L))
        ell = int(rng.integers(1, L + 1))
        r = rng.random()
        if r < 0.34 and L < 400:              # cap growth so the run stays bounded
            g._apply_duplication(s, ell)
            o.duplicate(s, ell)
        elif r < 0.67:
            g._apply_inversion(s, ell)
            o.invert(s, ell)
        else:
            g._apply_loss(s, min(ell, L - 1))
            o.delete(s, min(ell, L - 1))
        assert g.to_cells() == o.cells
        assert g.size() == len(o.cells)


def test_wrapping_duplication_matches_oracle():
    g = _fresh(20, ext=0.7)
    o = ArrayGenome("1", 20)
    for s, ell in [(18, 5), (3, 4), (15, 8), (0, 6)]:  # includes origin-crossing copies
        g._apply_duplication(s, ell)
        o.duplicate(s, ell)
        assert g.to_cells() == o.cells
        assert g.size() == len(o.cells)


def test_tandem_duplication_layout():
    g = _fresh(10)
    g._apply_duplication(2, 3)                 # copy of [2,3,4] lands right after
    cells = g.to_cells()
    assert [p for _s, p, _st in cells] == [0, 1, 2, 3, 4, 2, 3, 4, 5, 6, 7, 8, 9]
    assert g.size() == 13


def test_duplication_creates_paralogs_and_coalescences():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=6, age=2.5, seed=21)
    res = simulate_nucleotide_genomes(tree, inversion=0.004, duplication=0.006, loss=0.002,
                                      root_length=300, extension=0.95, seed=21)
    dups = [r for r in res.event_log if r.event is EventType.DUPLICATION]
    assert dups
    for r in dups:                             # each duplication is a bifurcation
        assert len(r.genes) == 3
        assert [op.role for op in r.genes] == ["parent", "left", "right"]
    _ids, _species, M = res.profile_matrix()
    assert M.max() >= 2                        # a surviving paralog -> copy number > 1


def test_duplication_grows_content_as_multiset():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=6, age=2.5, seed=22)
    res = simulate_nucleotide_genomes(tree, inversion=0.004, duplication=0.006, loss=0.001,
                                      root_length=300, extension=0.95, seed=22)
    root = set(range(300))
    grew = False
    for genome in res.leaf_genomes.values():
        cells = genome.to_cells()
        assert all(src == "1" and p in root for src, p, _st in cells)  # only ancestral nt
        if len(cells) > 300:
            grew = True                        # duplication lifted length above the root
    assert grew


def test_all_events_reproducible():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=7, age=2.5, seed=23)
    kw = dict(inversion=0.004, duplication=0.005, loss=0.004, root_length=250, extension=0.95)
    a = simulate_nucleotide_genomes(tree, **kw, seed=24)
    b = simulate_nucleotide_genomes(tree, **kw, seed=24)
    assert len(a.event_log) == len(b.event_log)
    for leaf, ga in a.leaf_genomes.items():
        assert ga.to_cells() == b.leaf_genomes[leaf].to_cells()


# --------------------------------------------------------------------------- #
# Per-atom gene trees — reconstruct the "gene" of each segment (steps 6-7)
# --------------------------------------------------------------------------- #
def _n_leaves(newick):
    return newick.count(",") + 1


def test_atom_gene_trees_match_profile_counts():
    """The reconciliation invariant: an atom's extant tree has one leaf per surviving copy."""
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=6, age=2.5, seed=31)
    res = simulate_nucleotide_genomes(tree, inversion=0.004, duplication=0.005, loss=0.004,
                                      root_length=300, extension=0.95, seed=31)
    ids, _species, M = res.profile_matrix()
    rowsum = {aid: int(M[i].sum()) for i, aid in enumerate(ids)}
    trees = res.atom_gene_trees()
    assert set(trees) == set(ids)
    for aid, (_complete, extant) in trees.items():
        assert extant is not None                 # atoms always survive in >= 1 leaf
        assert _n_leaves(extant) == rowsum[aid]    # leaves == total copies across species


def test_inversion_only_atom_trees_span_every_species():
    tree = simulate_species_tree(BirthDeath(1.0, 0.0), n_tips=6, age=2.0, seed=32)  # Yule: all extant
    res = simulate_nucleotide_genomes(tree, inversion=0.03, root_length=150, extension=0.9, seed=32)
    n_species = len(res.leaf_genomes)
    trees = res.atom_gene_trees()
    assert trees                                   # some atoms exist
    for _aid, (_complete, extant) in trees.items():
        assert _n_leaves(extant) == n_species      # present exactly once per species


def test_atom_gene_trees_record_losses_in_complete_tree():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=6, age=3.0, seed=33)
    res = simulate_nucleotide_genomes(tree, inversion=0.005, loss=0.02, root_length=300,
                                      extension=0.95, seed=33)
    assert any(r.event is EventType.LOSS for r in res.event_log)
    completes = [c for c, _e in res.atom_gene_trees().values() if c]
    assert any("LOSS" in c for c in completes)     # lost lineages appear in the complete tree


def test_atom_gene_trees_are_reproducible():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=6, age=2.5, seed=34)
    kw = dict(inversion=0.004, duplication=0.005, loss=0.004, root_length=250, extension=0.95)
    a = simulate_nucleotide_genomes(tree, **kw, seed=35).atom_gene_trees()
    b = simulate_nucleotide_genomes(tree, **kw, seed=35).atom_gene_trees()
    assert a == b


# --------------------------------------------------------------------------- #
# Transfer (HGT) — a copy travels to another lineage; discordant gene trees
# --------------------------------------------------------------------------- #
def test_transfer_extract_and_insert_mechanics():
    rng = np.random.default_rng(0)
    ids = IdManager()
    reg = SegmentRegistry()
    donor = NucleotideGenome(ids, root_length=20, extension=0.9, registry=reg)
    donor.originate(rng, TargetParams())                 # source "1"
    recipient = NucleotideGenome(ids, root_length=10, extension=0.9, registry=reg)
    recipient.originate(rng, TargetParams())             # source "2"

    seg = donor.extract_segment(Selection(genes=(), region=Region(0, 5, 4)), rng)
    assert donor.size() == 20                            # donor keeps a continuation (net 0)
    before = recipient.size()
    at = recipient.choose_insertion_point(seg, rng)
    recipient.insert_segment(seg, at, rng)
    assert recipient.size() == before + 4                # additive gain
    srcs = {src for src, _p, _st in recipient.to_cells()}
    assert srcs == {"1", "2"}                            # acquired donor material (source 1)


def test_transfer_fires_and_reconciliation_invariant_holds():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=6, age=2.5, seed=41)
    res = simulate_nucleotide_genomes(tree, inversion=0.003, transfer=0.006, loss=0.003,
                                      root_length=300, extension=0.95, seed=41)
    transfers = [r for r in res.event_log if r.event is EventType.TRANSFER]
    assert transfers
    for r in transfers:                                  # donor lineage forks; copy crosses over
        assert [op.role for op in r.genes] == ["parent", "donor_copy", "transfer_copy"]
        assert r.donor and r.recipient and r.donor != r.recipient
    ids, _sp, M = res.profile_matrix()
    rowsum = {aid: int(M[i].sum()) for i, aid in enumerate(ids)}
    for aid, (_c, extant) in res.atom_gene_trees().items():
        assert extant is not None and _n_leaves(extant) == rowsum[aid]  # xenologs counted


def test_transfer_is_additive_growth():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=6, age=2.5, seed=42)
    res = simulate_nucleotide_genomes(tree, transfer=0.008, loss=0.001, root_length=300,
                                      extension=0.95, seed=42)
    assert any(g.size() > 300 for g in res.leaf_genomes.values())  # acquired extra material


def test_transfer_reproducible():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=6, age=2.5, seed=43)
    kw = dict(inversion=0.003, transfer=0.006, loss=0.003, root_length=250, extension=0.95)
    a = simulate_nucleotide_genomes(tree, **kw, seed=44)
    b = simulate_nucleotide_genomes(tree, **kw, seed=44)
    assert len(a.event_log) == len(b.event_log)
    assert a.atom_gene_trees() == b.atom_gene_trees()
    for leaf, ga in a.leaf_genomes.items():
        assert ga.to_cells() == b.leaf_genomes[leaf].to_cells()


# --------------------------------------------------------------------------- #
# Transposition — cut-and-paste, content-preserving; vs the oracle
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("seed", range(20))
def test_transposition_matches_oracle(seed):
    rng = np.random.default_rng(seed)
    L = int(rng.integers(20, 120))
    g = _fresh(L, ext=0.9)
    o = ArrayGenome("1", L)
    for _ in range(80):
        s = int(rng.integers(L))
        ell = int(rng.integers(1, L))
        dest = int(rng.integers(L))
        g._apply_transposition(s, ell, dest)
        o.transpose(s, ell, dest)
        assert g.to_cells() == o.cells
        assert g.size() == L                       # content- and length-preserving


def test_transposition_conserves_content_at_leaves():
    tree = simulate_species_tree(BirthDeath(1.0, 0.0), n_tips=6, age=2.0, seed=51)
    res = simulate_nucleotide_genomes(tree, transposition=0.03, root_length=150,
                                      extension=0.9, seed=51)
    assert any(r.event is EventType.TRANSPOSITION for r in res.event_log)
    root = {("1", i) for i in range(150)}
    for genome in res.leaf_genomes.values():
        origins = {(src, p) for src, p, _st in genome.to_cells()}
        assert origins == root                     # pure permutation: nothing gained/lost


# --------------------------------------------------------------------------- #
# Origination — novel sequence under a fresh source namespace
# --------------------------------------------------------------------------- #
def test_origination_creates_new_sources_with_their_own_trees():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=6, age=2.5, seed=61)
    res = simulate_nucleotide_genomes(tree, origination=0.8, loss=0.001, root_length=200,
                                      extension=0.95, seed=61)
    assert sum(1 for r in res.event_log if r.event is EventType.ORIGINATION) > 1  # + the seed
    sources = {a.source for a in res.atoms}
    assert len(sources) > 1                         # novel sources beyond the root chromosome
    # every source's atoms reconstruct correctly (reconciliation invariant across sources)
    ids, _sp, M = res.profile_matrix()
    rowsum = {aid: int(M[i].sum()) for i, aid in enumerate(ids)}
    for aid, (_c, extant) in res.atom_gene_trees().items():
        assert extant is not None and _n_leaves(extant) == rowsum[aid]


def test_novel_gene_is_absent_from_lineages_that_predate_it():
    tree = simulate_species_tree(BirthDeath(1.0, 0.0), n_tips=6, age=2.0, seed=62)
    res = simulate_nucleotide_genomes(tree, origination=0.6, root_length=150, extension=0.95, seed=62)
    ids, _species, M = res.profile_matrix()
    # a source born mid-tree cannot be present in every species (unlike the root chromosome)
    by_source = {}
    for i, aid in enumerate(ids):
        a = next(x for x in res.atoms if x.atom_id == aid)
        by_source.setdefault(a.source, []).append(M[i])
    non_root = [s for s in by_source if s != "1"]
    assert non_root                                 # some novel sources exist
    assert any(row.min() == 0 for s in non_root for row in by_source[s])  # patchy presence
