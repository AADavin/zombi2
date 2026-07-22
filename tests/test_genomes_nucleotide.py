"""Tests for nucleotide genomes — the coordinate representation, inversion, the multi-chromosome,
identity-bearing karyotype wired along the species tree, and the copy-lineage genealogy log.

A chromosome is an ordered list of blocks (runs of one ancestry — not merged during the run: option
B, maximality is a recovered property). The strongest checks are the oracle (apply the same
inversions to a plain per-nucleotide array and require ``trace_back`` to agree), the strong invariant
(an ancestry-neutral event conserves ancestry, so every node carries the whole root sequence,
permuted), and — with loss / duplication — that the copy-lineage log accounts for every node's copy
numbers exactly and every copy traces back to a seed origination.
"""

import collections

import numpy as np
import pytest

from zombi2.genomes.events import node_from_label, node_label
from zombi2.species import simulate_species_tree
from zombi2.genomes import simulate_genomes_nucleotide
from zombi2.genomes.nucleotide import (
    Block,
    Chromosome,
    Duplication,
    Inversion,
    Loss,
    NucleotideGenome,
    Origination,
    Speciation,
    Transfer,
    Translocation,
    Transposition,
    _do_duplication,
    _do_loss,
    _do_translocation,
    _split_block,
)


def _chrom(length, topology="circular", cid=0, source=0):
    return Chromosome(cid, topology, [Block(source, 0, length, 1)])


def _ancestry(chrom):
    return sorted((src, pos) for (src, pos, _s) in chrom.trace_back())


# --- blocks + split -------------------------------------------------------------------------------

def test_split_block_is_strand_aware():
    assert _split_block(Block(0, 10, 20, 1), 3) == [Block(0, 10, 13, 1), Block(0, 13, 20, 1)]
    # reversed: the first o physical positions are the HIGH end of the source
    assert _split_block(Block(0, 10, 20, -1), 3) == [Block(0, 17, 20, -1), Block(0, 10, 17, -1)]


def test_chromosome_seed_and_readers():
    ch = _chrom(5)
    assert ch.length == 5
    assert ch.mosaic() == [(0, 0, 5, 1)]
    assert ch.trace_back() == [(0, 0, 1), (0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)]


def test_topology_validation():
    with pytest.raises(ValueError):
        Chromosome(0, "loop", [Block(0, 0, 5, 1)])


# --- inversion ------------------------------------------------------------------------------------

def test_inversion_hand_example():
    ch = _chrom(5)
    ch.invert(1, 3)
    assert ch.mosaic() == [(0, 0, 1, 1), (0, 1, 4, -1), (0, 4, 5, 1)]
    assert ch.trace_back() == [(0, 0, 1), (0, 3, -1), (0, 2, -1), (0, 1, -1), (0, 4, 1)]


def test_circular_inversion_wraps_the_origin():
    ch = _chrom(6, "circular")
    ch.invert(4, 4)                                          # arc 4,5,0,1 wraps the origin
    assert ch.trace_back() == [(0, 1, -1), (0, 0, -1), (0, 5, -1), (0, 4, -1), (0, 2, 1), (0, 3, 1)]
    assert ch.length == 6


def test_linear_inversion_clamps_and_never_wraps():
    ch = _chrom(6, "linear")
    ch.invert(4, 10)                                         # overruns the end -> clamp to [4, 6)
    assert ch.trace_back() == [(0, 0, 1), (0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 5, -1), (0, 4, -1)]


def test_inversion_flips_strands_in_the_arc_only():
    ch = _chrom(10)
    ch.invert(3, 4)
    assert [s for (_x, _p, s) in ch.trace_back()] == [1, 1, 1, -1, -1, -1, -1, 1, 1, 1]


# --- blocks are NOT merged during the run (option B: maximality is a recovered property) ----------

def test_inversion_leaves_breakpoints_but_conserves_trace():
    # An inversion and its inverse restore the ancestry exactly, but the seam stays split: blocks are
    # never merged during the run (that would fuse copy lineages) — maximality is recovered later.
    ch = _chrom(12)
    before = ch.trace_back()
    ch.invert(2, 6)
    ch.invert(2, 6)                                          # undone -> ancestry restored, seam remains
    assert ch.trace_back() == before                        # ancestry exactly restored
    assert len(ch.blocks) > 1                               # ...but the breakpoints are not merged away


def test_splits_preserve_copy_lineage():
    # A split is not a birth: both pieces keep the parent copy id.
    ch = Chromosome(0, "linear", [Block(0, 0, 10, 1, 7)])
    ch.invert(2, 6)
    assert len(ch.blocks) > 1
    assert {b.copy for b in ch.blocks} == {7}               # every piece is still copy 7


# --- the oracle: random inversions vs a per-nucleotide array (both topologies) --------------------

def _oracle_invert(arr, start, length, topology):
    n = len(arr)
    if n == 0:
        return
    if topology == "linear":
        start = max(0, min(start, n))
        end = min(start + max(0, length), n)
        if end <= start:
            return
        arr[start:end] = [(src, pos, -strand) for (src, pos, strand) in reversed(arr[start:end])]
        return
    ell = max(1, min(length, n))
    s = start % n
    if s + ell <= n:
        arr[s:s + ell] = [(src, pos, -strand) for (src, pos, strand) in reversed(arr[s:s + ell])]
    else:
        arr[:] = arr[s:] + arr[:s]
        arr[:ell] = [(src, pos, -strand) for (src, pos, strand) in reversed(arr[:ell])]


@pytest.mark.parametrize("topology", ["circular", "linear"])
def test_trace_back_matches_the_oracle_under_random_inversions(topology):
    rng = np.random.default_rng(0)
    for _trial in range(50):
        n = int(rng.integers(5, 40))
        ch = _chrom(n, topology)
        arr = [(0, p, 1) for p in range(n)]
        for _ in range(int(rng.integers(1, 15))):
            s, ell = int(rng.integers(0, n)), int(rng.integers(0, n))
            ch.invert(s, ell)
            _oracle_invert(arr, s, ell, topology)
            assert ch.trace_back() == arr                    # exact, every step
        assert ch.length == n


# --- the karyotype (NucleotideGenome) -------------------------------------------------------------

def test_karyotype_length_and_readers():
    g = NucleotideGenome([Chromosome(0, "circular", [Block(0, 0, 100, 1)]),
                          Chromosome(1, "linear", [Block(1, 0, 40, 1)])])
    assert g.length == 140
    assert set(g.mosaic()) == {0, 1}                         # keyed by chromosome id
    assert g.ancestry() == sorted([(0, p) for p in range(100)] + [(1, p) for p in range(40)])


# --- the tree wiring: heterogeneous seeding, identity, the strong invariant -----------------------

def _run(seed=1, inversion=2, chromosomes=1, **kw):
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=seed)
    params = dict(inversion=inversion, inversion_length=15, root_length=120, chromosomes=chromosomes,
                  seed=seed)
    params.update(kw)
    return sp, simulate_genomes_nucleotide(sp, **params)


def test_every_node_has_a_karyotype_of_conserved_length():
    sp, r = _run(seed=2, chromosomes=3, root_length=100)
    assert set(r.genomes) == set(sp.complete_tree.nodes)
    assert all(g.length == 300 for g in r.genomes.values())  # 3 x 100, inversion conserves length
    assert all(len(g.chromosomes) == 3 for g in r.genomes.values())   # no tier yet: count conserved


def test_heterogeneous_seeding_sizes_and_shapes():
    specs = [(100, "circular"), (40, "circular"), (25, "linear")]
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=4, seed=5)
    r = simulate_genomes_nucleotide(sp, inversion=0, chromosomes=specs, seed=5)
    tip = sorted(n.id for n in sp.complete_tree.extant())[0]
    chroms = r.genomes[tip].chromosomes
    assert [c.topology for c in chroms] == ["circular", "circular", "linear"]
    assert sorted(c.length for c in chroms) == [25, 40, 100]     # the three sizes, preserved
    # with no inversion, each replicon is a single block under its own source
    assert {c.blocks[0].source for c in chroms} == {0, 1, 2}


def test_every_node_carries_the_whole_root_sequence_permuted():
    # the strong invariant across chromosomes: inheritance copies ancestry, inversion changes none
    specs = [(100, "circular"), (40, "circular"), (25, "linear")]
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=3)
    r = simulate_genomes_nucleotide(sp, inversion=3, inversion_length=10, chromosomes=specs, seed=3)
    full = sorted((s, p) for s, (length, _t) in enumerate(specs) for p in range(length))
    for node_id in r.genomes:
        assert r.ancestry(node_id) == full


def test_chromosome_identity_network_is_well_formed():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=4)
    r = simulate_genomes_nucleotide(sp, inversion=1,
                                    chromosomes=[(80, "circular"), (30, "linear")], seed=4)
    ce = r.chromosome_events
    roots = [e for e in ce if e.kind == "origination"]
    spec = [e for e in ce if e.kind == "speciation"]
    assert len(roots) == 2                                    # one origination per seed replicon
    assert all(e.parents == () and len(e.children) == 1 for e in roots)
    assert spec and all(len(e.parents) == 1 and len(e.children) == 2 for e in spec)
    born = {}
    for e in ce:                                             # every chromosome id minted exactly once
        for ch in e.children:
            assert ch not in born
            born[ch] = e
    assert all(p in born for e in ce for p in e.parents)     # every speciation parent born earlier


def test_zero_inversion_leaves_each_replicon_a_single_block():
    sp, r = _run(seed=6, inversion=0, chromosomes=2, root_length=50)
    assert r.rearrangements == []
    for g in r.genomes.values():
        assert all(len(c.blocks) == 1 for c in g.chromosomes)


def test_inversions_recorded_within_their_branch_and_chromosome():
    sp, r = _run(seed=7, inversion=4, chromosomes=2, root_length=100)
    assert r.rearrangements
    node_chroms = {n: {c.id for c in r.genomes[n].chromosomes} for n in r.genomes}
    for inv in r.rearrangements:
        node = sp.complete_tree.nodes[inv.lineage]
        assert node.birth_time <= inv.time <= node.end_time
        assert inv.chromosome in node_chroms[inv.lineage]    # fired on a chromosome that node carried
        assert inv.length >= 1


def test_deterministic_given_seed():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=8)
    kw = dict(inversion=2.5, inversion_length=15, chromosomes=[(120, "circular"), (30, "linear")], seed=8)
    a = simulate_genomes_nucleotide(sp, **kw)
    b = simulate_genomes_nucleotide(sp, **kw)
    assert all(a.mosaic(n) == b.mosaic(n) for n in a.genomes)
    assert a.rearrangements == b.rearrangements
    assert a.chromosome_events == b.chromosome_events


def test_simulate_validation():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=5, seed=1)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, inversion=-1.0)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, inversion_length=0)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, chromosomes=0)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, chromosomes=[(0, "circular")])   # a bad replicon length
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, fission=-1.0)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, fusion=-1.0)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, translocation=-1.0)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, translocation_length=0)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, inversion_probability=2.0)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, loss=-1.0)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, loss_length=0)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, duplication=-1.0)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, duplication_length=0)


# --- duplication: the first birth (the event; the gene-tree recovery comes next) ------------------

def test_duplication_copies_an_arc_in_tandem():
    from collections import Counter
    g = NucleotideGenome([Chromosome(0, "circular", [Block(0, 0, 20, 1)])])
    events = []
    _do_duplication(g, 5, 1.0, 5.0, np.random.default_rng(0), events, _minter(100))
    assert g.length > 20                                    # grew
    assert isinstance(events[0], Duplication) and events[0].chromosome == 0
    counts = Counter((s, p) for (s, p, _st) in g.chromosomes[0].trace_back())
    copied = {(src, p) for (_pc, _cc, src, a, b) in events[0].copied for p in range(a, b)}
    assert all(counts[pos] == 2 for pos in copied)          # exactly the copied positions are doubled
    assert all(counts[(0, p)] == 1 for p in range(20) if (0, p) not in copied)


def test_duplication_keeps_all_ancestry_with_extra_copies():
    # duplication only removes nothing, so every root position is still present (set == full), but
    # the copied stretches now appear more than once
    specs = [(120, "circular"), (40, "linear")]
    full = {(s, p) for s, (length, _t) in enumerate(specs) for p in range(length)}
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=3)
    r = simulate_genomes_nucleotide(sp, duplication=2, duplication_length=10, inversion=1,
                                    chromosomes=specs, seed=3)
    assert any(isinstance(e, Duplication) for e in r.events)
    for node_id in r.genomes:
        assert set(r.ancestry(node_id)) == full             # nothing lost
    assert any(len(r.ancestry(n)) > len(full) for n in r.genomes)   # copies really exist


def test_duplication_with_loss_stays_a_subset():
    specs = [(120, "circular"), (40, "linear")]
    full = {(s, p) for s, (length, _t) in enumerate(specs) for p in range(length)}
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=4)
    r = simulate_genomes_nucleotide(sp, duplication=2, loss=2, inversion=1, translocation=1.5,
                                    fission=0.2, fusion=0.2, chromosomes=specs, seed=4)
    for node_id in r.genomes:
        assert set(r.ancestry(node_id)) <= full             # loss can now remove; dup adds copies


# --- loss: the first ancestry-changing event (a death) -------------------------------------------

def test_loss_deletes_an_arc_and_records_the_lost_material():
    g = NucleotideGenome([Chromosome(0, "circular", [Block(0, 0, 20, 1)])])
    events = []
    _do_loss(g, 5, 1.0, 5.0, np.random.default_rng(0), events)
    assert 1 <= g.length < 20                                # shrank, never emptied
    assert isinstance(events[0], Loss) and events[0].chromosome == 0
    survived = {(s, p) for (s, p, _st) in g.chromosomes[0].trace_back()}
    lost = {(src, p) for (_cp, src, a, b) in events[0].lost for p in range(a, b)}
    assert survived | lost == {(0, p) for p in range(20)}   # a clean partition of the original
    assert survived.isdisjoint(lost)


def test_loss_is_a_noop_on_a_length_one_chromosome():
    g = NucleotideGenome([Chromosome(0, "circular", [Block(0, 0, 1, 1)])])
    events = []
    _do_loss(g, 5, 1.0, 5.0, np.random.default_rng(0), events)
    assert events == [] and g.length == 1                   # never empties a chromosome


def test_loss_weakens_the_invariant_monotonically():
    # loss only ever removes: each node's ancestry is a subset of the root (no duplicates), and a
    # child's ancestry is a subset of its parent's — monotone non-increasing down every path
    full = {(s, p) for s, (length, _t) in enumerate(_TIER_SPECS) for p in range(length)}
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=3)
    r = simulate_genomes_nucleotide(sp, loss=2.5, loss_length=10, inversion=1,
                                    chromosomes=_TIER_SPECS, seed=3)
    assert r.events                                          # losses really happened
    for node_id in r.genomes:
        anc_list = r.ancestry(node_id)
        anc = set(anc_list)
        assert anc <= full                                  # subset of the root
        assert len(anc_list) == len(anc)                    # no duplicates
        node = sp.complete_tree.nodes[node_id]
        if node.parent is not None:
            assert anc <= set(r.ancestry(node.parent))      # child ⊆ parent (monotone)
    assert any(set(r.ancestry(n)) < full for n in r.genomes)   # something really was lost


def test_loss_composes_with_everything_and_still_only_subsets():
    full = {(s, p) for s, (length, _t) in enumerate(_TIER_SPECS) for p in range(length)}
    for seed in range(3):
        sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=seed)
        r = simulate_genomes_nucleotide(sp, loss=2, inversion=1, translocation=1.5,
                                        fission=0.3, fusion=0.3, chromosomes=_TIER_SPECS, seed=seed)
        for node_id in r.genomes:
            anc_list = r.ancestry(node_id)
            assert set(anc_list) <= full and len(anc_list) == len(set(anc_list))


# --- step 2a: the chromosome tier (fission + fusion -> the reticulating network) ------------------

_TIER_SPECS = [(200, "circular"), (80, "circular"), (40, "linear")]


def _minter(start):
    box = [start]

    def mint():
        box[0] += 1
        return box[0]
    return mint


def _tier(seed, **kw):
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=seed)
    params = dict(inversion=1, inversion_length=12, fission=0.6, fusion=0.6,
                  chromosomes=_TIER_SPECS, seed=seed)
    params.update(kw)
    return sp, simulate_genomes_nucleotide(sp, **params)


def test_fission_conserves_length_and_re_mints_both_children():
    from zombi2.genomes.nucleotide import _do_fission
    g = NucleotideGenome([Chromosome(0, "linear", [Block(0, 0, 10, 1)])])
    ce = []
    _do_fission(g, 5, 1.0, np.random.default_rng(1), ce, _minter(10))
    assert len(g.chromosomes) == 2 and all(c.topology == "linear" for c in g.chromosomes)
    assert sum(c.length for c in g.chromosomes) == 10           # length conserved
    assert ce[0].kind == "fission" and ce[0].parents == (0,) and len(set(ce[0].children)) == 2


def test_fission_circular_makes_two_rings():
    from zombi2.genomes.nucleotide import _do_fission
    g = NucleotideGenome([Chromosome(0, "circular", [Block(0, 0, 12, 1)])])
    ce = []
    _do_fission(g, 5, 1.0, np.random.default_rng(0), ce, _minter(10))
    assert len(g.chromosomes) == 2 and all(c.topology == "circular" for c in g.chromosomes)
    assert sum(c.length for c in g.chromosomes) == 12


def test_fusion_is_a_noop_on_a_mixed_topology_pair():
    from zombi2.genomes.nucleotide import _do_fusion
    g = NucleotideGenome([Chromosome(0, "circular", [Block(0, 0, 10, 1)]),
                          Chromosome(1, "linear", [Block(1, 0, 10, 1)])])
    ce = []
    _do_fusion(g, 0, 1.0, np.random.default_rng(0), ce, _minter(99))
    assert ce == [] and len(g.chromosomes) == 2                 # nothing to fuse (no same-topology partner)


def test_fusion_merges_a_same_topology_pair():
    from zombi2.genomes.nucleotide import _do_fusion
    g = NucleotideGenome([Chromosome(0, "circular", [Block(0, 0, 10, 1)]),
                          Chromosome(1, "circular", [Block(1, 0, 8, 1)])])
    ce = []
    _do_fusion(g, 3, 2.0, np.random.default_rng(0), ce, _minter(41))
    assert len(g.chromosomes) == 1 and g.chromosomes[0].topology == "circular"
    assert g.chromosomes[0].length == 18                        # blocks concatenated, length conserved
    assert ce[0].kind == "fusion" and len(ce[0].parents) == 2 and len(ce[0].children) == 1


def test_tier_fires_with_the_right_arity():
    _sp, r = _tier(seed=5)
    assert {"origination", "speciation", "fission", "fusion"} <= {e.kind for e in r.chromosome_events}
    shape = {"origination": (0, 1), "speciation": (1, 2), "fission": (1, 2), "fusion": (2, 1)}
    for e in r.chromosome_events:
        assert (len(e.parents), len(e.children)) == shape[e.kind]


def test_fusion_is_the_only_reticulation():
    _sp, r = _tier(seed=6)
    for e in r.chromosome_events:
        assert (len(e.parents) == 2) == (e.kind == "fusion")


def test_chromosome_network_well_formed_with_reticulation():
    _sp, r = _tier(seed=7)
    born = {}
    for e in r.chromosome_events:
        for ch in e.children:
            assert ch not in born                               # each chromosome id minted exactly once
            born[ch] = e
    assert all(p in born for e in r.chromosome_events for p in e.parents)


def test_strong_invariant_survives_the_tier():
    full = sorted((s, p) for s, (length, _t) in enumerate(_TIER_SPECS) for p in range(length))
    for seed in range(4):
        _sp, r = _tier(seed=seed)
        for node_id in r.genomes:
            assert r.ancestry(node_id) == full                  # ancestry-neutral: whole sequence at every node


def test_tier_changes_the_chromosome_number():
    _sp, r = _tier(seed=5)
    assert len({len(g.chromosomes) for g in r.genomes.values()}) > 1


def test_tier_is_deterministic():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=8)
    kw = dict(inversion=1, fission=0.5, fusion=0.5, chromosomes=_TIER_SPECS, seed=8)
    a = simulate_genomes_nucleotide(sp, **kw)
    b = simulate_genomes_nucleotide(sp, **kw)
    assert a.chromosome_events == b.chromosome_events
    assert all(a.mosaic(n) == b.mosaic(n) for n in a.genomes)


# --- step 2b: translocation (an arc moving between chromosomes) -----------------------------------

def _two_chroms():
    return NucleotideGenome([Chromosome(0, "circular", [Block(0, 0, 20, 1)]),
                             Chromosome(1, "circular", [Block(1, 0, 10, 1)])])


def test_translocation_moves_an_arc_between_chromosomes_conserving_ancestry():
    g = _two_chroms()
    rearr = []
    _do_translocation(g, 5, 1.0, 5.0, 0.0, np.random.default_rng(0), rearr)
    assert len(g.chromosomes) == 2                            # count unchanged
    assert sum(c.length for c in g.chromosomes) == 30         # total length conserved
    anc = sorted((s, p) for c in g.chromosomes for (s, p, _st) in c.trace_back())
    assert anc == sorted([(0, p) for p in range(20)] + [(1, p) for p in range(10)])   # ancestry-neutral
    assert isinstance(rearr[0], Translocation) and rearr[0].source != rearr[0].dest


def test_translocation_can_flip_with_inversion_probability():
    g = _two_chroms()
    rearr = []
    _do_translocation(g, 5, 1.0, 5.0, 1.0, np.random.default_rng(0), rearr)   # always flip
    assert rearr[0].flipped is True


def test_translocation_is_a_noop_with_a_single_chromosome():
    g = NucleotideGenome([Chromosome(0, "circular", [Block(0, 0, 20, 1)])])
    rearr = []
    _do_translocation(g, 5, 1.0, 5.0, 0.0, np.random.default_rng(0), rearr)
    assert rearr == [] and len(g.chromosomes) == 1


def test_translocation_conserves_the_chromosome_count():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=5)
    r = simulate_genomes_nucleotide(sp, translocation=3, translocation_length=10,
                                    chromosomes=_TIER_SPECS, seed=5)
    assert all(len(g.chromosomes) == 3 for g in r.genomes.values())   # material moves, count does not
    assert any(isinstance(x, Translocation) for x in r.rearrangements)


def test_strong_invariant_survives_translocation():
    full = sorted((s, p) for s, (length, _t) in enumerate(_TIER_SPECS) for p in range(length))
    for seed in range(3):
        sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=seed)
        r = simulate_genomes_nucleotide(sp, translocation=4, translocation_length=8, inversion=1,
                                        inversion_probability=0.5, chromosomes=_TIER_SPECS, seed=seed)
        for node_id in r.genomes:
            assert r.ancestry(node_id) == full


def test_all_four_events_compose():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=5)
    r = simulate_genomes_nucleotide(sp, inversion=1, translocation=2.5, fission=0.4, fusion=0.4,
                                    chromosomes=_TIER_SPECS, seed=5)
    assert {"Inversion", "Translocation"} <= {type(x).__name__ for x in r.rearrangements}
    assert {"fission", "fusion"} <= {e.kind for e in r.chromosome_events}
    full = sorted((s, p) for s, (length, _t) in enumerate(_TIER_SPECS) for p in range(length))
    assert all(r.ancestry(node_id) == full for node_id in r.genomes)


# --- the copy-lineage genealogy log (the input the gene-tree recovery will read) ------------------

def _log_copy_number(r, node_id):
    """The per-(source, position) copy number implied by the genealogy log alone: walk the path from
    the root to ``node_id``, adding a copy per origination / duplication row that covers the position
    and removing one per loss row, on the matching branch (rearrangements and speciation are copy-
    number neutral, so they play no part)."""
    from collections import Counter
    tree = r.complete_tree
    on_path, n = set(), node_id
    while n is not None:
        on_path.add(n)
        n = tree.nodes[n].parent
    count: Counter = Counter()
    for e in r.events:
        if e.lineage not in on_path:
            continue
        if isinstance(e, Origination):
            for p in range(e.start, e.end):
                count[(e.source, p)] += 1
        elif isinstance(e, Duplication):
            for (_pc, _cc, src, a, b) in e.copied:
                for p in range(a, b):
                    count[(src, p)] += 1
        elif isinstance(e, Loss):
            for (_cp, src, a, b) in e.lost:
                for p in range(a, b):
                    count[(src, p)] -= 1
    return {k: v for k, v in count.items() if v}


def test_copy_number_matches_the_genealogy_log():
    # The genealogy log accounts for the genome exactly: at every node the copy number of each
    # ancestral position (from trace_back) equals the number implied by the originations,
    # duplications and losses on the path from the root. End-to-end validation of the log's accounting.
    from collections import Counter
    specs = [(80, "circular"), (30, "linear")]
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=7)
    r = simulate_genomes_nucleotide(sp, duplication=2.5, loss=2.5, duplication_length=8,
                                    loss_length=8, inversion=1.5, translocation=1.5, fission=0.15,
                                    fusion=0.15, chromosomes=specs, seed=7)
    assert any(isinstance(e, Duplication) for e in r.events)   # duplications really fired
    assert any(isinstance(e, Loss) for e in r.events)          # ...and losses
    for node_id in r.genomes:
        genome = Counter((s, p) for chrom in r.genomes[node_id].chromosomes
                         for (s, p, _st) in chrom.trace_back())
        assert {k: v for k, v in genome.items() if v} == _log_copy_number(r, node_id)


def test_every_copy_lineage_is_born_and_traces_to_a_root():
    # Every copy lineage present anywhere is born by a recorded event and its parentage traces back
    # to a seed origination — no orphans, no cycles, and the unset sentinel (0) never leaks.
    specs = [(80, "circular"), (30, "linear")]
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=9)
    r = simulate_genomes_nucleotide(sp, duplication=2.5, loss=1.5, inversion=1.5,
                                    translocation=1.5, fission=0.15, fusion=0.15, chromosomes=specs,
                                    seed=9)
    roots = {e.copy for e in r.events if isinstance(e, Origination)}
    parent: dict[int, int] = {}
    for e in r.events:
        if isinstance(e, Duplication):
            for (pc, cc, *_rest) in e.copied:
                parent[cc] = pc
        elif isinstance(e, Speciation):
            for c in e.children:
                parent[c] = e.parent

    def traces_to_root(cp):
        seen = set()
        while cp not in roots:
            if cp in seen or cp not in parent:
                return False
            seen.add(cp)
            cp = parent[cp]
        return True

    present = {b.copy for g in r.genomes.values() for chrom in g.chromosomes for b in chrom.blocks}
    assert present and 0 not in present                        # copies exist; no sentinel leaked
    assert all(traces_to_root(cp) for cp in present)           # each is born and reaches a seed root


def test_speciation_re_mints_every_copy_lineage():
    # At a speciation each parent copy is re-minted into one fresh copy per daughter (recorded as a
    # Speciation), so a daughter node shares no copy id with its parent node.
    specs = [(60, "circular")]
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=6, seed=2)
    r = simulate_genomes_nucleotide(sp, inversion=1.5, duplication=1.5, chromosomes=specs, seed=2)
    assert [e for e in r.events if isinstance(e, Speciation)]   # speciations recorded
    tree = r.complete_tree
    for node_id, g in r.genomes.items():
        node = tree.nodes[node_id]
        if node.children is None:
            continue
        parent_copies = {b.copy for chrom in g.chromosomes for b in chrom.blocks}
        for c in node.children:
            child_copies = {b.copy for chrom in r.genomes[c].chromosomes for b in chrom.blocks}
            assert parent_copies.isdisjoint(child_copies)      # fully re-minted at the speciation


# --- the gene-tree recovery: root partition -> one tree per block ---------------------------------

def _tips(node):
    out, stack = [], [node]
    while stack:
        n = stack.pop()
        if n.children:
            stack.extend(n.children)
        else:
            out.append(n)
    return out


def test_recovery_ancestry_neutral_is_the_species_tree():
    # With no birth/death of copies (only ancestry-neutral events), every extant leaf carries the
    # whole root once, so every root-block's extant gene tree is the species tree: one tip per leaf.
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=1)
    r = simulate_genomes_nucleotide(sp, inversion=2.5, inversion_length=12, translocation=1.5,
                                    fission=0.1, fusion=0.1, chromosomes=[(60, "circular")], seed=1)
    leaves = {n.id for n in r.complete_tree.extant()}
    assert r.root_blocks                                       # inversions leave surviving breakpoints
    for fam, _blk in enumerate(r.root_blocks):
        ex = r.gene_trees[fam].extant
        tip_species = [t.species for t in _tips(ex) if t.kind == "extant"]
        assert sorted(tip_species) == sorted(leaves)           # exactly one copy per extant leaf


def test_recovered_extant_leaves_match_observed_copies():
    # The cross-check: for every root-block and every extant leaf, the number of extant tips the
    # recovered gene tree puts on that leaf equals the number of copies actually in the leaf genome.
    import collections
    specs = [(80, "circular"), (30, "linear")]
    for seed in range(4):
        sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=seed)
        r = simulate_genomes_nucleotide(sp, duplication=2.5, loss=2.5, duplication_length=8,
                                        loss_length=8, inversion=1.5, translocation=1.5, fission=0.1,
                                        fusion=0.1, chromosomes=specs, seed=seed)
        assert any(isinstance(e, Duplication) for e in r.events)
        assert any(isinstance(e, Loss) for e in r.events)
        leaves = [n.id for n in r.complete_tree.extant()]
        for fam, (s, a, b) in enumerate(r.root_blocks):
            ex = r.gene_trees[fam].extant
            recovered = collections.Counter(t.species for t in (_tips(ex) if ex else [])
                                            if t.kind == "extant")
            for lid in leaves:
                observed = sum(1 for chrom in r.genomes[lid].chromosomes for blk in chrom.blocks
                               if blk.source == s and blk.start <= a and b <= blk.end)
                assert observed == recovered.get(lid, 0)


def test_every_root_block_has_a_parsing_gene_tree():
    specs = [(60, "circular"), (25, "linear")]
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=8)
    r = simulate_genomes_nucleotide(sp, duplication=2.5, loss=2, inversion=1.5, chromosomes=specs,
                                    seed=8)
    assert set(r.gene_trees) == set(range(len(r.root_blocks)))    # one family per root-block, by index
    assert {s for (s, _a, _b) in r.root_blocks} <= {0, 1}        # blocks live under their seed sources
    for fam in r.gene_trees:
        nwk = r.gene_trees[fam].to_newick("complete")
        assert nwk and nwk.endswith(";") and nwk.count("(") == nwk.count(")")


def test_duplication_shows_up_as_a_branch_in_a_block_tree():
    # A run with duplication (and nothing that removes copies) must produce at least one block whose
    # extant tree has a duplication node — a leaf carrying two copies of that block.
    sp = simulate_species_tree(birth=1.0, death=0.1, n_extant=6, seed=3)
    r = simulate_genomes_nucleotide(sp, duplication=3, duplication_length=10,
                                    chromosomes=[(80, "circular")], seed=3)

    def has_kind(node, kind):
        stack = [node]
        while stack:
            n = stack.pop()
            if n.kind == kind:
                return True
            stack.extend(n.children)
        return False

    assert any(has_kind(r.gene_trees[fam].extant, "duplication")
               for fam in r.gene_trees if r.gene_trees[fam].extant is not None)


# --- transfer: the global-timeline coupling (a horizontal birth) ----------------------------------

_XFER_SPECS = [(80, "circular"), (30, "linear")]
_XFER_FULL = {(s, p) for s, (length, _t) in enumerate(_XFER_SPECS) for p in range(length)}


def test_transfer_is_additive_keeps_all_ancestry_and_grows_copies():
    # Additive transfer only adds copies (the donor keeps its own), and every lineage shares the root
    # sources, so with no loss the ancestry set stays full while copy numbers grow.
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=10, seed=2)
    r = simulate_genomes_nucleotide(sp, transfer=3.5, transfer_length=10, inversion=1,
                                    chromosomes=_XFER_SPECS, seed=2)
    assert any(isinstance(e, Transfer) for e in r.events)          # transfers really fired
    for node_id in r.genomes:
        assert set(r.ancestry(node_id)) == _XFER_FULL              # nothing lost
    assert any(len(r.ancestry(n)) > len(_XFER_FULL) for n in r.genomes)   # extra copies exist


def test_transfer_recipient_is_a_distinct_contemporaneous_lineage():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=4)
    r = simulate_genomes_nucleotide(sp, transfer=4, transfer_length=8, chromosomes=_XFER_SPECS, seed=4)
    xfers = [e for e in r.events if isinstance(e, Transfer)]
    assert xfers
    tree = r.complete_tree
    for e in xfers:
        assert e.recipient != e.lineage                            # no self-transfer by default
        donor, recip = tree.nodes[e.lineage], tree.nodes[e.recipient]
        assert donor.birth_time < e.time < donor.end_time          # fired while the donor was alive
        assert recip.birth_time < e.time < recip.end_time          # ...and the recipient too (contemporaries)


def test_recovery_cross_check_holds_with_transfer():
    # The end-to-end check, now with transfer creating horizontal edges: recovered extant tips still
    # equal the copies actually present in every extant leaf, for every root-block.
    import collections
    for seed in range(4):
        sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=seed)
        r = simulate_genomes_nucleotide(sp, transfer=2.5, transfer_length=8, duplication=1.5,
                                        loss=1.5, inversion=1.5, translocation=1, fission=0.1,
                                        fusion=0.1, chromosomes=_XFER_SPECS, seed=seed)
        assert any(isinstance(e, Transfer) for e in r.events)
        leaves = [n.id for n in r.complete_tree.extant()]
        for fam, (s, a, b) in enumerate(r.root_blocks):
            ex = r.gene_trees[fam].extant
            recovered = collections.Counter(t.species for t in (_tips(ex) if ex else [])
                                            if t.kind == "extant")
            for lid in leaves:
                observed = sum(1 for chrom in r.genomes[lid].chromosomes for blk in chrom.blocks
                               if blk.source == s and blk.start <= a and b <= blk.end)
                assert observed == recovered.get(lid, 0)


def test_transfer_is_a_horizontal_edge_in_a_block_tree():
    sp = simulate_species_tree(birth=1.0, death=0.1, n_extant=8, seed=3)
    r = simulate_genomes_nucleotide(sp, transfer=4, transfer_length=10, chromosomes=_XFER_SPECS, seed=3)

    def has_transfer(node):
        stack = [node]
        while stack:
            n = stack.pop()
            if n.kind == "transfer":
                return True
            stack.extend(n.children)
        return False

    assert any(has_transfer(r.gene_trees[fam].complete) for fam in r.gene_trees)


def test_transfer_to_distance_and_self_transfer_run():
    from zombi2.genomes.nucleotide import Distance
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=5)
    for kw in (dict(transfer_to="distance"), dict(transfer_to=Distance(decay=2.0)),
               dict(self_transfer=True)):
        r = simulate_genomes_nucleotide(sp, transfer=3, transfer_length=8, chromosomes=_XFER_SPECS,
                                        seed=5, **kw)
        assert any(isinstance(e, Transfer) for e in r.events)
        for node_id in r.genomes:                              # additive: full ancestry preserved
            assert set(r.ancestry(node_id)) == _XFER_FULL


def test_transfer_validation():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=4, seed=1)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, transfer=-0.1)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, transfer=5, transfer_length=0)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, transfer=5, transfer_to="nearest")


# --- transposition: an intra-chromosome move (ancestry-neutral) -----------------------------------

def test_transposition_is_ancestry_neutral_and_conserves_counts():
    specs = [(80, "circular"), (30, "linear")]
    full = sorted((s, p) for s, (length, _t) in enumerate(specs) for p in range(length))
    for seed in range(3):
        sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=seed)
        r = simulate_genomes_nucleotide(sp, transposition=4, transposition_length=8,
                                        inversion_probability=0.5, chromosomes=specs, seed=seed)
        assert any(isinstance(x, Transposition) for x in r.rearrangements)
        for node_id in r.genomes:
            assert r.ancestry(node_id) == full                 # conserves ancestry (permutes only)
            assert len(r.genomes[node_id].chromosomes) == 2    # ...and the chromosome count


def test_transposition_works_within_a_single_chromosome():
    # Unlike translocation (needs >=2 chromosomes), transposition moves an arc within one chromosome.
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=6, seed=2)
    r = simulate_genomes_nucleotide(sp, transposition=5, transposition_length=8, chromosomes=1,
                                    root_length=100, seed=2)
    assert any(isinstance(x, Transposition) for x in r.rearrangements)
    full = sorted((0, p) for p in range(100))
    assert all(r.ancestry(n) == full for n in r.genomes)       # single replicon, ancestry conserved


# --- origination: a de-novo birth (a fresh source, a new gene family) -----------------------------

def test_origination_adds_de_novo_sources_beyond_the_root():
    specs = [(80, "circular"), (30, "linear")]
    root_full = {(s, p) for s, (length, _t) in enumerate(specs) for p in range(length)}
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=3)
    r = simulate_genomes_nucleotide(sp, origination=0.5, origination_length=10, chromosomes=specs, seed=3)
    denovo = [e for e in r.events if isinstance(e, Origination) and e.source >= len(specs)]
    assert denovo                                              # fresh sources really arose
    assert all(e.source >= len(specs) for e in denovo)         # ...numbered past the seed sources
    # the seed material is never removed by origination; some node carries new sources too
    assert all(root_full <= set(r.ancestry(n)) for n in r.genomes)
    assert any(any(s >= len(specs) for (s, _p) in r.ancestry(n)) for n in r.genomes)


def test_origination_mints_a_gene():
    # origination lays down a GENE, not plain spacer: indivisible, its own family, its own span.
    specs = [(80, "circular"), (30, "linear")]
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=6)
    r = simulate_genomes_nucleotide(sp, origination=0.5, origination_length=10, chromosomes=specs, seed=6)
    denovo = [e for e in r.events if isinstance(e, Origination) and e.source >= len(specs)]
    assert denovo
    # every de-novo source is carried by genic blocks only, and is registered as a gene span
    denovo_sources = {e.source for e in denovo}
    for g in r.genomes.values():
        for chrom in g.chromosomes:
            for b in chrom.blocks:
                if b.source in denovo_sources:
                    assert b.is_gene                          # never plain spacer
    assert denovo_sources <= {src for (src, _a, _b) in r.gene_spans.values()}


def test_origination_family_roots_at_its_own_branch():
    # Each de-novo gene's tree is rooted on the branch the origination fired on (not the tree root).
    specs = [(80, "circular"), (30, "linear")]
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=6)
    r = simulate_genomes_nucleotide(sp, origination=0.5, origination_length=10, inversion=1,
                                    chromosomes=specs, seed=6)
    origin_branch = {e.source: e.lineage for e in r.events
                     if isinstance(e, Origination) and e.source >= len(specs)}
    assert origin_branch
    seen = 0
    for fam, gt in r.gene_trees.items():
        src, _a, _b = r.gene_spans[fam]
        if src in origin_branch:
            assert gt.complete.species == origin_branch[src]  # rooted at the origination branch
            seen += 1
    assert seen                                               # some de-novo gene survived


def test_recovery_cross_check_holds_with_origination():
    import collections
    specs = [(80, "circular"), (30, "linear")]
    for seed in range(3):
        sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=seed)
        r = simulate_genomes_nucleotide(sp, origination=0.4, origination_length=10, loss=1.5,
                                        duplication=1.5, inversion=1.5, transposition=1,
                                        chromosomes=specs, seed=seed)
        assert any(isinstance(e, Origination) and e.source >= len(specs) for e in r.events)
        leaves = [n.id for n in r.complete_tree.extant()]
        assert r.gene_trees
        for fam, gt in r.gene_trees.items():
            s, a, b = r.gene_spans[fam]
            recovered = collections.Counter(t.species for t in (_tips(gt.extant) if gt.extant else [])
                                            if t.kind == "extant")
            for lid in leaves:
                observed = sum(1 for chrom in r.genomes[lid].chromosomes for blk in chrom.blocks
                               if blk.source == s and blk.start <= a and b <= blk.end)
                assert observed == recovered.get(lid, 0)


def test_transposition_and_origination_validation():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=4, seed=1)
    for kw in (dict(transposition=-0.1), dict(transposition=5, transposition_length=0),
               dict(origination=-0.1), dict(origination=0.1, origination_length=0)):
        with pytest.raises(ValueError):
            simulate_genomes_nucleotide(sp, **kw)


# --- the chromosome tier: de-novo replicons (origination) and whole-chromosome death (loss) --------

def test_chromosome_origination_adds_a_replicon_carrying_a_gene():
    specs = [(80, "circular"), (30, "linear")]
    root_full = {(s, p) for s, (length, _t) in enumerate(specs) for p in range(length)}
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=3)
    r = simulate_genomes_nucleotide(sp, chromosome_origination=0.4, origination_length=25,
                                    chromosomes=specs, seed=3)
    orig_edges = [e for e in r.chromosome_events if e.kind == "origination"]
    assert len(orig_edges) > len(specs)                            # de-novo replicons beyond the seeds
    assert all(e.parents == () and len(e.children) == 1 for e in orig_edges)   # ...are network roots
    # a de-novo replicon is never born empty: it carries one new gene on a fresh source of its own
    # (the seeded replicons here were declared with no genes at all, so they legitimately have none)
    de_novo = {c for e in orig_edges if e.time > 0.0 for c in e.children}
    assert de_novo
    assert all(c.length > 0 and c.n_genes >= 1
               for g in r.genomes.values() for c in g.chromosomes if c.id in de_novo)
    for node_id in r.genomes:
        assert root_full <= set(r.ancestry(node_id))               # the seed material is untouched...
    assert any(any(s >= len(specs) for (s, _p) in r.ancestry(n)) for n in r.genomes)  # ...plus new
    assert any(len(g.chromosomes) > len(specs) for g in r.genomes.values())   # the count grew


def test_chromosome_loss_kills_whole_chromosomes_as_a_subset():
    specs = [(80, "circular"), (30, "linear")]
    full = {(s, p) for s, (length, _t) in enumerate(specs) for p in range(length)}
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=4)
    r = simulate_genomes_nucleotide(sp, chromosome_loss=0.25, chromosomes=specs, seed=4)
    assert any(e.kind == "loss" for e in r.chromosome_events)       # whole chromosomes died
    assert all(len(g.chromosomes) >= 1 for g in r.genomes.values())  # never the last chromosome
    tree = r.complete_tree
    for node_id, node in tree.nodes.items():
        assert set(r.ancestry(node_id)) <= full                    # material is only ever removed
        if node.parent is not None:                                # ...monotonically down every path
            assert set(r.ancestry(node_id)) <= set(r.ancestry(node.parent))
    assert any(len(r.ancestry(n)) < len(full) for n in r.genomes)  # some material really was lost


def test_chromosome_tier_network_is_well_formed_with_de_novo_and_death():
    specs = [(80, "circular"), (40, "circular")]
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=5)
    r = simulate_genomes_nucleotide(sp, chromosome_origination=0.25, chromosome_loss=0.15,
                                    fission=0.25, fusion=0.25, chromosomes=specs, seed=5)
    ev = r.chromosome_events
    minted = [cid for e in ev for cid in e.children]
    assert len(minted) == len(set(minted))                         # every chromosome id minted once
    for e in ev:
        if e.kind == "origination":
            assert e.parents == () and len(e.children) == 1        # a root (no parent)
        elif e.kind == "loss":
            assert e.children == () and len(e.parents) == 1        # a leaf (no child)
        elif e.kind == "fusion":
            assert len(e.parents) == 2 and len(e.children) == 1    # the only reticulation
        else:                                                       # speciation / fission
            assert len(e.parents) == 1 and len(e.children) == 2
    assert {"origination", "loss", "speciation"} <= {e.kind for e in ev}


def test_recovery_cross_check_holds_with_chromosome_tier():
    import collections
    specs = [(80, "circular"), (30, "linear")]
    saw_chromosome_loss = False
    for seed in range(4):
        sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=seed)
        r = simulate_genomes_nucleotide(sp, chromosome_origination=0.25, chromosome_loss=0.12,
                                        loss=1.5, duplication=1.5, transfer=1.5, inversion=1.5,
                                        fission=0.1, fusion=0.1, chromosomes=specs, seed=seed)
        saw_chromosome_loss = saw_chromosome_loss or any(e.kind == "loss" for e in r.chromosome_events)
        leaves = [n.id for n in r.complete_tree.extant()]
        for fam, gt in r.gene_trees.items():               # de-novo replicons carry genes: genic mode
            s, a, b = r.gene_spans[fam]
            recovered = collections.Counter(t.species for t in (_tips(gt.extant) if gt.extant else [])
                                            if t.kind == "extant")
            for lid in leaves:
                observed = sum(1 for chrom in r.genomes[lid].chromosomes for blk in chrom.blocks
                               if blk.source == s and blk.start <= a and b <= blk.end)
                assert observed == recovered.get(lid, 0)
    assert saw_chromosome_loss                                     # chromosome deaths really occurred


def test_chromosome_tier_validation():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=4, seed=1)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, chromosome_origination=-0.1)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, chromosome_loss=-0.1)


# --- the genic layer: declared genes are indivisible ----------------------------------------------

_GENIC = dict(genes=4, gene_length=30, chromosomes=[(300, "circular")])


def _gene_spans(r):
    """``{family: {(source, start, end), …}}`` over every block of every node. A family with more than
    one span means a gene was cut somewhere — the thing that must never happen."""
    import collections
    spans = collections.defaultdict(set)
    for g in r.genomes.values():
        for chrom in g.chromosomes:
            for b in chrom.blocks:
                if b.gene:
                    spans[b.gene].add((b.source, b.start, b.end))
    return spans


def test_split_at_refuses_to_cut_a_gene():
    from zombi2.genomes.nucleotide import _CutsGene
    ch = Chromosome(0, "linear", [Block(0, 0, 10, 1, 1), Block(0, 10, 20, 1, 1, 7)])  # intergene + gene
    ch._split_at(5)                                          # inside the intergene: fine
    assert len(ch.blocks) == 3
    with pytest.raises(_CutsGene):
        ch._split_at(15)                                     # strictly inside the gene: refused
    assert len(ch.blocks) == 3                               # ...and nothing was mutated
    with pytest.raises(_CutsGene):
        ch._check_cut(15)                                    # the pure test agrees
    assert len(ch.blocks) == 3                               # _check_cut mutates nothing
    ch._split_at(10)                                         # exactly the gene's edge: allowed
    assert len(ch.blocks) == 3


def test_seeding_lays_down_the_alternating_chain():
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=4, seed=1)
    r = simulate_genomes_nucleotide(sp, genes=6, gene_length=40, chromosomes=1, root_length=600, seed=1)
    root = r.genomes[r.complete_tree.root].chromosomes[0]
    kinds = [("gene" if b.is_gene else "intergene") for b in root.blocks]
    assert kinds == ["intergene", "gene"] * 6                # I G I G … the declared chain
    assert [b.length for b in root.blocks if b.is_gene] == [40] * 6
    assert sum(b.length for b in root.blocks) == 600         # the replicon length is preserved
    assert len({b.gene for b in root.blocks if b.is_gene}) == 6      # six distinct families


def test_no_genes_declared_is_the_uniform_model():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=2)
    r = simulate_genomes_nucleotide(sp, inversion=2.5, loss=1.5, chromosomes=1, root_length=200, seed=2)
    assert all(not b.is_gene for g in r.genomes.values() for c in g.chromosomes for b in c.blocks)


def test_genes_are_never_split_under_every_event():
    for seed in range(2):
        sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=seed)
        r = simulate_genomes_nucleotide(
            sp, inversion=2.5, transposition=1.5, translocation=1.5, loss=2, duplication=2,
            transfer=1.5, fission=0.1, fusion=0.1, inversion_probability=0.4, inversion_length=20,
            loss_length=20, duplication_length=20, transfer_length=20, seed=seed, **_GENIC)
        spans = _gene_spans(r)
        assert spans                                         # genes really are present
        for fam, seen in spans.items():
            assert len(seen) == 1, f"gene family {fam} was cut into {seen}"
            (_src, a, b), = seen
            assert b - a == 30                               # every copy is the whole gene


def test_declared_genes_come_back_as_intact_root_blocks():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=3)
    r = simulate_genomes_nucleotide(sp, inversion=2.5, transposition=1.5, loss=2, duplication=2,
                                    inversion_length=30, loss_length=30, seed=3, **_GENIC)
    root_genes = {(b.source, b.start, b.end) for chrom in r.genomes[r.complete_tree.root].chromosomes
                  for b in chrom.blocks if b.is_gene}
    assert root_genes
    assert root_genes <= set(r.root_blocks)                  # each gene is one whole root-block


def test_gene_spans_records_every_declared_gene():
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=4, seed=1)
    r = simulate_genomes_nucleotide(sp, genes=5, gene_length=40, chromosomes=1, root_length=500, seed=1)
    assert len(r.gene_spans) == 5                            # one span per declared gene
    assert all(e - a == 40 for (_s, a, e) in r.gene_spans.values())
    root = r.genomes[r.complete_tree.root].chromosomes[0]
    assert {b.gene: (b.source, b.start, b.end) for b in root.blocks if b.is_gene} == r.gene_spans


def test_gene_trees_are_one_per_gene_not_per_root_block():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=3)
    r = simulate_genomes_nucleotide(sp, inversion=2.5, transposition=1.5, loss=2.0, duplication=2.0,
                                    inversion_length=30, loss_length=30, duplication_length=30,
                                    seed=3, **_GENIC)
    assert r.gene_trees                                      # genes declared -> genic families only
    assert set(r.gene_trees) <= set(r.gene_spans)            # keyed by GENE FAMILY id
    # the intergenic root-blocks are recovered as blocks but are not built into trees
    assert len(r.root_blocks) > len(r.gene_trees)
    for fam in r.gene_trees:
        assert r.gene_spans[fam] in set(r.root_blocks)       # each gene is a whole root-block


def test_recovery_cross_check_holds_with_genes():
    # the leaves-==-observed-copies check, now per declared gene
    import collections
    for seed in range(2):
        sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=5, seed=seed)
        r = simulate_genomes_nucleotide(
            sp, inversion=1.0, transposition=0.5, loss=1.0, duplication=1.0, transfer=0.5,
            inversion_length=20, loss_length=20, duplication_length=20, seed=seed, **_GENIC)
        leaves = [n.id for n in r.complete_tree.extant()]
        assert r.gene_trees
        for fam, gt in r.gene_trees.items():
            s, a, b = r.gene_spans[fam]
            recovered = collections.Counter(t.species for t in (_tips(gt.extant) if gt.extant else [])
                                            if t.kind == "extant")
            for lid in leaves:
                observed = sum(1 for chrom in r.genomes[lid].chromosomes for blk in chrom.blocks
                               if blk.source == s and blk.start <= a and b <= blk.end)
                assert observed == recovered.get(lid, 0)


def test_gene_declaration_validation():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=4, seed=1)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, genes=-1)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, genes=3, gene_length=0)
    with pytest.raises(ValueError):                          # no room left for intergenes
        simulate_genomes_nucleotide(sp, genes=10, gene_length=100, chromosomes=1, root_length=500)


# --- declaring the seed genome from a GFF ---------------------------------------------------------

_GFF_TEXT = ("##gff-version 3\n"
             "##sequence-region chrom1 1 3000\n"
             "##sequence-region plasmid 1 800\n"
             "chrom1\tZOMBI2\tgene\t201\t500\t.\t+\t.\tID=dnaA\n"
             "chrom1\tZOMBI2\tgene\t900\t1400\t.\t-\t.\tID=recA\n"
             "chrom1\tZOMBI2\tCDS\t900\t1400\t.\t-\t0\tParent=recA\n"
             "plasmid\tZOMBI2\tgene\t51\t250\t.\t+\t.\tID=toxin\n")


def _gff(tmp_path):
    path = tmp_path / "seed.gff"
    path.write_text(_GFF_TEXT)
    return path


def test_seeding_from_a_gff(tmp_path):
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=4, seed=1)
    r = simulate_genomes_nucleotide(sp, gff=_gff(tmp_path), seed=1)
    # the GFF supplies the replicons too, in sorted seqid order: chrom1 then plasmid
    chroms = r.genomes[r.complete_tree.root].chromosomes
    assert [c.length for c in chroms] == [3000, 800]
    assert r.gene_names == {"dnaA": 1, "recA": 2, "toxin": 3}
    assert r.gene_spans == {1: (0, 200, 500), 2: (0, 899, 1400), 3: (1, 50, 250)}
    # the GFF's strand is the gene's CODING strand (annotation), recorded separately...
    assert r.gene_strands == {1: 1, 2: -1, 3: 1}
    # ...while every seed block is +1: Block.strand is orientation relative to the ancestral source,
    # and at the root nothing has been inverted yet.
    assert all(b.strand == 1 for c in chroms for b in c.blocks)
    # everything between the genes is intergene, and the chain covers the replicon exactly
    for c in chroms:
        assert sum(b.length for b in c.blocks) == c.length


def test_gff_genes_are_never_split(tmp_path):
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=2)
    r = simulate_genomes_nucleotide(sp, gff=_gff(tmp_path), inversion=2.0, transposition=1.0,
                                    translocation=1.0, loss=1.5, duplication=1.5, transfer=1.0,
                                    inversion_length=100, loss_length=100, duplication_length=100,
                                    transfer_length=100, inversion_probability=0.4, seed=2)
    spans = _gene_spans(r)
    assert spans
    for fam, seen in spans.items():
        assert len(seen) == 1, f"gene family {fam} was cut into {seen}"
        assert seen == {r.gene_spans[fam]}               # ...and it is exactly the declared span


def test_gff_genes_get_one_tree_each_lookupable_by_name(tmp_path):
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=3)
    r = simulate_genomes_nucleotide(sp, gff=_gff(tmp_path), inversion=2.0, loss=1.5, duplication=1.5,
                                    inversion_length=100, loss_length=100, duplication_length=100, seed=3)
    assert set(r.gene_trees) <= set(r.gene_spans)
    tree = r.gene_trees[r.gene_names["dnaA"]]            # look a named gene up by name
    assert tree.to_newick("extant").endswith(";")


def test_gff_and_genes_are_mutually_exclusive(tmp_path):
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=4, seed=1)
    with pytest.raises(ValueError, match="either gff= or genes="):
        simulate_genomes_nucleotide(sp, gff=_gff(tmp_path), genes=3, seed=1)


# --- drawing arcs directly from the legal set (no guess-and-redraw) -------------------------------

def _dense_chromosome(n_genes=40, gene_len=94, gap=6):
    """A gene-dense replicon: 94% genic, like a real bacterial genome."""
    blocks, at, fam = [], 0, 0
    for _ in range(n_genes):
        blocks.append(Block(0, at, at + gap, 1, 1)); at += gap
        fam += 1
        blocks.append(Block(0, at, at + gene_len, 1, 1, fam)); at += gene_len
    return Chromosome(0, "circular", blocks), at


def _inside_a_gene(chrom, c):
    pos = 0
    for b in chrom.blocks:
        if pos < c < pos + b.length:
            return b.is_gene
        pos += b.length
    return False


def test_arc_extent_always_lands_on_a_legal_breakpoint():
    chrom, total = _dense_chromosome()
    rng = np.random.default_rng(0)
    starts = [p for p in range(total) if not _inside_a_gene(chrom, p)]
    drawn = 0
    for _ in range(400):
        start = starts[int(rng.integers(len(starts)))]
        d = chrom._pick_arc_extent(start, 300, rng)
        if d is None:
            continue
        drawn += 1
        assert d >= 1
        assert not _inside_a_gene(chrom, (start + d) % total)   # the far end is always legal
    assert drawn > 300                                          # and it essentially always finds one


def test_legal_cut_never_lands_inside_a_gene():
    chrom, _total = _dense_chromosome()
    rng = np.random.default_rng(1)
    for _ in range(300):
        c = chrom._pick_legal_cut(rng)
        assert c is not None and not _inside_a_gene(chrom, c)


def test_extent_tracks_the_mean_asked_for():
    # the realised extent is conditioned by the gene structure (so shorter than asked), but a larger
    # mean must still give larger arcs
    chrom, total = _dense_chromosome(n_genes=60)
    rng = np.random.default_rng(2)
    starts = [p for p in range(total) if not _inside_a_gene(chrom, p)]

    def mean_extent(mean):
        got = [chrom._pick_arc_extent(starts[int(rng.integers(len(starts)))], mean, rng)
               for _ in range(500)]
        got = [d for d in got if d is not None]
        return sum(got) / len(got)

    assert mean_extent(100) < mean_extent(1000) < mean_extent(5000)


def test_event_rate_is_honoured_on_a_gene_dense_genome():
    # THE regression this replaced: with guess-and-redraw, a large extent on a dense genome exhausted
    # its retries and silently dropped up to ~43% of events. The rate must now come out as asked.
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=5, seed=1)
    lineage_time = sum(n.end_time - n.birth_time for n in sp.complete_tree.nodes.values())
    for extent in (200, 20_000):                     # short arcs and arcs spanning many genes
        counts = []
        for seed in range(12):
            r = simulate_genomes_nucleotide(sp, genes=100, gene_length=94, chromosomes=1,
                                            root_length=10_000, loss=3.0, loss_length=extent, seed=seed)
            counts.append(sum(1 for e in r.events if isinstance(e, Loss)))
        expected = 3.0 * lineage_time
        assert 0.7 * expected < sum(counts) / len(counts) < 1.4 * expected, (extent, counts)


def test_with_no_events_the_genome_is_exactly_what_was_declared(tmp_path):
    """The round trip: run with every rate at zero and each leaf must BE the input genome.

    In particular the coordinate space must be the identity — position i traces back to source
    position i, forward — including across genes declared on the minus strand. `Block.strand` records
    inversion relative to the ancestral source, and at the root nothing has been inverted; the GFF's
    strand is the gene's coding strand, which is annotation and lives in `gene_strands`."""
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=4, seed=1)
    r = simulate_genomes_nucleotide(sp, gff=_gff(tmp_path), seed=1)          # no events at all
    for lid in (n.id for n in r.complete_tree.extant()):
        chrom1, plasmid = r.genomes[lid].chromosomes
        assert (chrom1.length, plasmid.length) == (3000, 800)
        assert chrom1.trace_back() == [(0, i, 1) for i in range(3000)]       # the identity map
        assert plasmid.trace_back() == [(1, i, 1) for i in range(800)]
        assert all(b.strand == 1 for c in (chrom1, plasmid) for b in c.blocks)
    # the declared genes are exactly where the GFF put them, with their coding strands preserved
    assert r.gene_spans == {1: (0, 200, 500), 2: (0, 899, 1400), 3: (1, 50, 250)}
    assert r.gene_strands == {1: 1, 2: -1, 3: 1}


def test_the_even_layout_also_round_trips():
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=4, seed=2)
    r = simulate_genomes_nucleotide(sp, genes=5, gene_length=40, chromosomes=1, root_length=500, seed=2)
    for lid in (n.id for n in r.complete_tree.extant()):
        chrom, = r.genomes[lid].chromosomes
        assert chrom.trace_back() == [(0, i, 1) for i in range(500)]
    assert set(r.gene_strands.values()) == {1}


def test_inversions_along_the_tree_match_an_independent_replay():
    """The strongest end-to-end check: replay the recorded inversions on a plain per-nucleotide array,
    node by node down the tree, and demand the engine's genome match exactly at EVERY node.

    The existing oracle tests one chromosome in isolation; this one validates the whole pipeline —
    genic seeding, the tree wiring, inheritance at speciation, and the block algebra — against a naive
    simulator that knows nothing about blocks."""
    import collections
    length, n_genes, gene_len = 2_000, 12, 120
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=5, seed=1)
    r = simulate_genomes_nucleotide(sp, genes=n_genes, gene_length=gene_len, chromosomes=1,
                                    root_length=length, inversion=25.0, inversion_length=400, seed=4)
    inversions = [e for e in r.rearrangements if isinstance(e, Inversion)]
    assert len(inversions) > 15                              # the run really did rearrange

    def invert(arr, start, ell):
        n = len(arr)
        ell = max(1, min(ell, n))
        s = start % n
        if s + ell <= n:
            arr[s:s + ell] = [(src, p, -st) for (src, p, st) in reversed(arr[s:s + ell])]
        else:
            arr[:] = arr[s:] + arr[:s]
            arr[:ell] = [(src, p, -st) for (src, p, st) in reversed(arr[:ell])]

    tree = r.complete_tree
    by_node = collections.defaultdict(list)
    for e in inversions:
        by_node[e.lineage].append(e)
    stack = [(tree.root, [(0, i, 1) for i in range(length)])]
    while stack:
        nid, arr = stack.pop()
        arr = list(arr)
        for e in sorted(by_node[nid], key=lambda x: x.time):
            invert(arr, e.start, e.length)
        assert r.genomes[nid].chromosomes[0].trace_back() == arr, f"node {nid} diverged from the replay"
        if tree.nodes[nid].children:
            for c in tree.nodes[nid].children:
                stack.append((c, arr))

    full = sorted((0, i) for i in range(length))
    for leaf in tree.extant():                               # nothing gained, nothing lost
        assert r.ancestry(leaf.id) == full
        seen = {b.gene for b in r.genomes[leaf.id].chromosomes[0].blocks if b.is_gene}
        assert len(seen) == n_genes                          # every gene still there, still whole


# --- two chromosomes evolving with the chromosome tier, genes and all --------------------------

_TWO = dict(chromosomes=[(2000, "circular"), (1200, "circular")], genes=8, gene_length=120)
_TWO_FULL = sorted((s, p) for s, (length, _t) in enumerate(_TWO["chromosomes"]) for p in range(length))


def test_two_chromosomes_ancestry_neutral_tier_conserves_everything():
    """Two replicons, genes on both, evolving with the chromosome tier *and* arcs moving between
    chromosomes. Every event here is ancestry-neutral, so the karyotype may be reshaped arbitrarily —
    split, merged, material shuffled across replicons — but not one base may be gained or lost, and no
    gene may be cut."""
    import collections
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=5, seed=1)
    r = simulate_genomes_nucleotide(
        sp, inversion=6.0, inversion_length=300, translocation=6.0, translocation_length=300,
        transposition=3.0, transposition_length=300, fission=3.0, fusion=3.0, seed=5, **_TWO)

    kinds = collections.Counter(e.kind for e in r.chromosome_events)
    assert kinds["fission"] > 5 and kinds["fusion"] > 5       # the tier really ran
    assert any(isinstance(x, Translocation) for x in r.rearrangements)

    for node_id in r.genomes:                                  # nothing gained, nothing lost, anywhere
        assert r.ancestry(node_id) == _TWO_FULL

    karyotypes = {len(r.genomes[leaf.id].chromosomes) for leaf in r.complete_tree.extant()}
    assert len(karyotypes) > 1 and karyotypes != {2}           # the karyotype diverged across leaves

    # translocation/fusion mix the two replicons: a chromosome ends up carrying both sources
    assert any(len({b.source for b in c.blocks}) == 2
               for leaf in r.complete_tree.extant()
               for c in r.genomes[leaf.id].chromosomes)

    spans = _gene_spans(r)
    assert len(spans) == 16                                    # 8 genes on each replicon
    for fam, seen in spans.items():
        assert seen == {r.gene_spans[fam]}                     # every gene whole, where it was declared
    for leaf in r.complete_tree.extant():                      # ...and every gene still present
        present = {b.gene for c in r.genomes[leaf.id].chromosomes for b in c.blocks if b.is_gene}
        assert present == set(r.gene_spans)


def test_two_chromosomes_with_the_whole_event_set():
    """The same karyotype with births, deaths and whole-chromosome death piled on: the accounting must
    still close and the recovered gene trees must still match the genomes."""
    import collections
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=5, seed=1)
    r = simulate_genomes_nucleotide(
        sp, inversion=4.0, inversion_length=300, translocation=4.0, translocation_length=300,
        transposition=2.0, transposition_length=300, fission=2.0, fusion=2.0,
        chromosome_origination=1.0, chromosome_loss=0.25,
        loss=3.0, loss_length=250, duplication=3.0, duplication_length=250,
        transfer=2.0, transfer_length=250, seed=1, **_TWO)

    kinds = collections.Counter(e.kind for e in r.chromosome_events)
    assert {"fission", "fusion", "origination", "loss", "speciation"} <= set(kinds)
    assert all(len(g.chromosomes) >= 1 for g in r.genomes.values())    # never the last chromosome
    # no chromosome is ever left without a gene, so none is ever empty
    assert all(c.n_genes >= 1 for g in r.genomes.values() for c in g.chromosomes)

    minted = [cid for e in r.chromosome_events for cid in e.children]  # the network stays well-formed
    assert len(minted) == len(set(minted))
    for e in r.chromosome_events:
        if e.kind == "origination":
            assert e.parents == () and len(e.children) == 1
        elif e.kind == "loss":
            assert e.children == () and len(e.parents) == 1
        elif e.kind == "fusion":
            assert len(e.parents) == 2 and len(e.children) == 1
        else:
            assert len(e.parents) == 1 and len(e.children) == 2

    for fam, seen in _gene_spans(r).items():                   # genes survive the tier intact
        assert seen == {r.gene_spans[fam]}

    leaves = [n.id for n in r.complete_tree.extant()]          # the cross-check, per gene
    assert len(r.gene_trees) > 10                              # most genes survive: a real check
    for fam, gt in r.gene_trees.items():
        s, a, b = r.gene_spans[fam]
        recovered = collections.Counter(t.species for t in (_tips(gt.extant) if gt.extant else [])
                                        if t.kind == "extant")
        for lid in leaves:
            observed = sum(1 for chrom in r.genomes[lid].chromosomes for blk in chrom.blocks
                           if blk.source == s and blk.start <= a and b <= blk.end)
            assert observed == recovered.get(lid, 0)


# --- extinction in the species tree ---------------------------------------------------------------

def test_extinct_lineages_evolve_donate_and_are_pruned():
    """With extinction, whole lineages die. Three things must hold at once: the dead still evolve and
    are still recorded; what lived only on them never reaches the observable data; but what they
    *donated* to a surviving lineage does — the classic transfer-from-the-dead."""
    import collections
    sp = simulate_species_tree(birth=1.4, death=0.7, n_extant=5, seed=4)
    tree = sp.complete_tree
    extant = {n.id for n in tree.extant()}

    def doomed(nid):                                          # no extant descendant anywhere below
        node = tree.nodes[nid]
        return (nid not in extant) if node.children is None else all(doomed(c) for c in node.children)

    dead = {nid for nid in tree.nodes if doomed(nid)}
    assert collections.Counter(n.fate for n in tree.nodes.values())["extinct"] >= 1
    assert dead and extant and not (dead & extant)

    r = simulate_genomes_nucleotide(sp, chromosomes=[(1500, "circular")], genes=10, gene_length=100,
                                    inversion=3.0, inversion_length=200, loss=3.0, loss_length=200,
                                    duplication=3.0, duplication_length=200,
                                    transfer=6.0, transfer_length=200, seed=4)

    # the dead still evolved, and we kept their genomes
    assert set(r.genomes) == set(tree.nodes)
    assert all(r.genomes[d].length > 0 for d in dead)

    # a transfer OUT of a doomed lineage into one that survives: donated by a ghost
    from_dead = [e for e in r.events if isinstance(e, Transfer)
                 and e.lineage in dead and e.recipient not in dead]
    assert from_dead

    # what lived only on doomed lineages is not observable — but it is still history, so it still
    # gets a tree (every node votes on the partition), one with no extant tip to show for it
    in_extant = {b.gene for lid in extant for c in r.genomes[lid].chromosomes
                 for b in c.blocks if b.is_gene}
    in_dead = {b.gene for d in dead for c in r.genomes[d].chromosomes for b in c.blocks if b.is_gene}
    only_dead = in_dead - in_extant
    assert only_dead                                          # some gene really did die with them
    assert only_dead <= set(r.gene_trees)
    for fam in only_dead:
        assert r.gene_trees[fam].extant is None               # ...and none of them is observable
        assert r.gene_trees[fam].to_newick("complete").endswith(";")

    # the complete tree keeps the dead; the extant tree prunes to survivors only
    assert any(t.kind == "extinct" for gt in r.gene_trees.values() for t in _tips(gt.complete))
    for gt in r.gene_trees.values():
        if gt.extant is not None:
            assert all(t.kind == "extant" for t in _tips(gt.extant))

    # and the cross-check still closes, against the EXTANT leaves
    for fam, gt in r.gene_trees.items():
        s, a, b = r.gene_spans[fam]
        recovered = collections.Counter(t.species for t in (_tips(gt.extant) if gt.extant else [])
                                        if t.kind == "extant")
        for lid in extant:
            observed = sum(1 for chrom in r.genomes[lid].chromosomes for blk in chrom.blocks
                           if blk.source == s and blk.start <= a and b <= blk.end)
            assert observed == recovered.get(lid, 0)


# --- pathological genomes -------------------------------------------------------------------------

_HOT = dict(inversion=8.0, inversion_length=200, translocation=4.0, translocation_length=200,
            transposition=4.0, transposition_length=200, loss=6.0, loss_length=200,
            duplication=6.0, duplication_length=200, transfer=4.0, transfer_length=200,
            fission=2.0, fusion=2.0, chromosome_origination=1.0, chromosome_loss=0.5,
            origination=1.0, origination_length=50)


def _invariants_hold(r):
    """Genes never split, and the recovered trees still match the genomes."""
    import collections
    spans = _gene_spans(r)
    assert all(v == {r.gene_spans[f]} for f, v in spans.items())
    extant = [n.id for n in r.complete_tree.extant()]
    for fam, gt in r.gene_trees.items():
        s, a, b = r.gene_spans[fam]
        recovered = collections.Counter(t.species for t in (_tips(gt.extant) if gt.extant else [])
                                        if t.kind == "extant")
        for lid in extant:
            observed = sum(1 for chrom in r.genomes[lid].chromosomes for blk in chrom.blocks
                           if blk.source == s and blk.start <= a and b <= blk.end)
            assert observed == recovered.get(lid, 0)


@pytest.mark.parametrize("label, kwargs", [
    ("a 1 bp replicon",            dict(chromosomes=1, root_length=1)),
    ("a 2 bp replicon",            dict(chromosomes=1, root_length=2)),
    ("a 10 bp replicon, gene of 5", dict(chromosomes=1, root_length=10, genes=1, gene_length=5)),
    ("genes of a single base",     dict(chromosomes=1, root_length=100, genes=20, gene_length=1)),
    ("a linear replicon",          dict(chromosomes=1, root_length=200, genes=3, gene_length=40,
                                        topology="linear")),
])
def test_degenerate_genomes_survive_every_event(label, kwargs):
    """Absurdly small genomes, run flat out with every event switched on. Nothing may crash, no gene
    may be cut, and the recovery must still agree with the genomes."""
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=4, seed=1)
    r = simulate_genomes_nucleotide(sp, seed=1, **kwargs, **_HOT)
    assert all(g.length >= 1 for g in r.genomes.values())      # a genome is never wiped out entirely
    _invariants_hold(r)


def test_extents_far_larger_than_the_genome():
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=4, seed=1)
    r = simulate_genomes_nucleotide(sp, chromosomes=1, root_length=500, genes=4, gene_length=50,
                                    inversion=8.0, inversion_length=10**6,
                                    loss=6.0, loss_length=10**6,
                                    duplication=6.0, duplication_length=10**6, seed=1)
    assert all(g.length >= 1 for g in r.genomes.values())
    _invariants_hold(r)


def test_a_single_leaf_tree():
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=1, seed=1)
    r = simulate_genomes_nucleotide(sp, chromosomes=1, root_length=500, genes=4, gene_length=50,
                                    seed=1, **_HOT)
    assert len(list(r.complete_tree.extant())) == 1
    _invariants_hold(r)


def test_a_fully_genic_genome_rearranges_at_its_gene_boundaries(tmp_path):
    """A genome with **no intergenic base at all** still evolves. Every breakpoint falls on a boundary
    between two genes, so genes are inverted, moved, duplicated and lost **whole** and none is split.

    This test used to assert the opposite — that such a genome was frozen — because the sampling was
    stricter than the rule it implements. The rule is that a cut may not fall *strictly inside* a gene;
    a gene's own edge has always been legal, and the code that *lands* an arc always knew it. The code
    that *chose* one said "in the spacer" instead, so with no spacer every event was a silent no-op,
    whatever the rates. Both now read the one legal cut set."""
    path = tmp_path / "allgenic.gff"
    path.write_text("##gff-version 3\n##sequence-region c 1 600\n" +
                    "".join(f"c\t.\tgene\t{1 + 100 * i}\t{100 * (i + 1)}\t.\t+\t.\tID=g{i}\n"
                            for i in range(6)))
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=4, seed=1)
    r = simulate_genomes_nucleotide(sp, gff=path, **_HOT, seed=1)
    assert r.rearrangements
    assert [e for e in r.events if isinstance(e, (Loss, Duplication, Transfer))]
    # every breakpoint is a gene boundary, so every block is exactly one declared gene and nothing
    # was ever cut. (Genome lengths are not multiples of 100: origination mints de-novo genes of its
    # own geometric length, and those are whole genes too.)
    spans = set(r.gene_spans.values())
    for g in r.genomes.values():
        for c in g.chromosomes:
            for b in c.blocks:
                assert b.is_gene and (b.source, b.start, b.end) in spans
    _invariants_hold(r)


def test_a_genome_packed_with_genes_needs_no_gap_between_them():
    """`genes=` used to demand room for intergenes, so ten 100 bp genes in 1000 bp was refused. There
    is no such requirement: genes may sit flush, and the run then breaks at the joins between them."""
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=6, seed=1)
    r = simulate_genomes_nucleotide(sp, root_length=1000, genes=10, gene_length=100, inversion=8.0,
                                    inversion_length=200, duplication=8.0, loss=8.0, seed=1)
    assert not [b for c in r.initial_genome.chromosomes for b in c.blocks if not b.is_gene]
    assert r.rearrangements and [e for e in r.events if isinstance(e, Duplication)]
    assert all(g.length % 100 == 0 for g in r.genomes.values())
    with pytest.raises(ValueError, match="do not fit"):        # 11 x 100 genuinely does not fit
        simulate_genomes_nucleotide(sp, root_length=1000, genes=11, gene_length=100, seed=1)


def test_legal_cuts_include_gene_edges_and_empty_replicons():
    """A breakpoint is illegal only *strictly inside* a gene, so a gene's own edge is a legal cut, and
    an empty replicon has position 0. Sampling used to be stricter than the rule, which left de-novo
    plasmids unable to ever receive material and froze fully-genic chromosomes completely."""
    rng = np.random.default_rng(0)
    assert Chromosome(0, "circular", [])._pick_legal_cut(rng) == 0          # an empty replicon
    genic = Chromosome(1, "circular", [Block(0, 0, 100, 1, 1, 1), Block(0, 100, 200, 1, 1, 2)])
    assert {genic._pick_legal_cut(rng) for _ in range(200)} == {0, 100}     # only the gene edges
    mixed = Chromosome(2, "circular", [Block(0, 0, 10, 1, 1), Block(0, 10, 60, 1, 1, 1)])
    seen = {mixed._pick_legal_cut(rng) for _ in range(2000)}
    assert seen == set(range(11))              # the spacer's interior plus the gene's left edge only


def test_a_de_novo_replicon_is_born_with_a_gene_and_can_grow():
    """A plasmid from `chromosome_origination` is born carrying a gene, and material can still arrive
    on it afterwards."""
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=5, seed=3)
    r = simulate_genomes_nucleotide(sp, chromosomes=[(800, "circular")], genes=4, gene_length=80,
                                    chromosome_origination=3.0, translocation=6.0,
                                    translocation_length=150, transfer=4.0, transfer_length=150,
                                    origination=2.0, origination_length=60, inversion=2.0,
                                    inversion_length=150, seed=3)
    seeds = {c for e in r.chromosome_events if e.kind == "origination" and e.time == 0.0
             for c in e.children}
    de_novo = {c for e in r.chromosome_events if e.kind == "origination" and e.time > 0.0
               for c in e.children}
    assert de_novo
    born = {c.id for g in r.genomes.values() for c in g.chromosomes if c.id in de_novo}
    assert born and not (born & seeds)
    assert all(c.n_genes >= 1 for g in r.genomes.values() for c in g.chromosomes)
    # and they are not frozen at their birth size: some grew past their single gene
    assert any(c.length > 0 and c.n_genes >= 1 and len(c.blocks) > 1
               for g in r.genomes.values() for c in g.chromosomes if c.id in de_novo)
    _invariants_hold(r)


def test_no_chromosome_is_ever_left_without_a_gene():
    """A chromosome never exists without a gene. A replicon is born with one, and any event that would
    strip a chromosome of its last — a loss, a translocation carrying it away, a fission splitting off
    a geneless half — simply does not happen. Run flat out, with every event that could violate it."""
    for seed in range(4):
        sp = simulate_species_tree(birth=1.2, death=0.4, n_extant=6, seed=seed)
        r = simulate_genomes_nucleotide(
            sp, chromosomes=[(800, "circular"), (600, "circular")], genes=4, gene_length=100,
            loss=12.0, loss_length=400,                  # brutal loss: tries hard to empty things
            translocation=8.0, translocation_length=400, # tries to carry the last gene away
            fission=6.0, fusion=3.0,                     # tries to split off a geneless half
            transposition=4.0, transposition_length=200, inversion=4.0, inversion_length=200,
            duplication=3.0, duplication_length=200, transfer=3.0, transfer_length=200,
            chromosome_origination=2.0, chromosome_loss=1.0, origination_length=80, seed=seed)
        for node_id, g in r.genomes.items():
            assert g.chromosomes                          # a genome always keeps a chromosome
            assert g.length > 0                           # ...and a lineage never ends up with no DNA
            for c in g.chromosomes:
                assert c.n_genes >= 1, f"node {node_id} chromosome {c.id} lost its last gene"
                assert c.length > 0                       # so a chromosome is never empty either
        _invariants_hold(r)


# --- writing the outputs --------------------------------------------------------------------------

_ALL_OUTPUTS = ("events", "blocks", "genes", "chromosome_events")


def _read(path):
    lines = path.read_text().splitlines()
    cols = lines[0].split("\t")
    rows = [dict(zip(cols, row.split("\t"))) for row in lines[1:] if row]
    for row in lines[1:]:
        assert len(row.split("\t")) == len(cols), f"{path.name} is ragged: {row!r}"
    return cols, rows


def _written(tmp_path, **kw):
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=6, seed=3)
    r = simulate_genomes_nucleotide(sp, root_length=400, inversion=1.5, inversion_length=40,
                                    duplication=2.0, loss=1.5, transfer=1.5, origination=1.0,
                                    seed=4, **{"transposition": 1.0, **kw})
    r.write(tmp_path, outputs=_ALL_OUTPUTS)
    return r


def test_write_emits_the_selected_outputs(tmp_path):
    _written(tmp_path)
    assert {p.name for p in tmp_path.iterdir()} == {
        "genome_events.tsv", "blocks.tsv", "genes.tsv", "chromosome_events.tsv"}


def test_write_defaults_to_every_table(tmp_path):
    # as at the ordered resolution, a default run writes the whole history — blocks.tsv included,
    # big though it is: a replayable history missing one of its tables is not one
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=4, seed=1)
    simulate_genomes_nucleotide(sp, root_length=200, inversion=1.0, seed=1).write(tmp_path)
    written = {p.name for p in tmp_path.iterdir()}
    assert {"genome_events.tsv", "genes.tsv", "blocks.tsv", "initial_genome.tsv",
            "chromosome_events.tsv"} <= written
    assert any(p.name.startswith("gene_tree_fam") for p in tmp_path.iterdir())


def test_written_blocks_tile_every_chromosome_of_every_node(tmp_path):
    # the layout file must be the genome, not a summary of it: each chromosome's rows have to run
    # end to end from 0 with no gap and no overlap, and reproduce the in-memory blocks exactly
    r = _written(tmp_path)
    _, rows = _read(tmp_path / "blocks.tsv")

    by_chromosome = {}
    for row in rows:
        by_chromosome.setdefault((node_from_label(row["lineage"]), int(row["chromosome"])), []).append(row)

    seen = set()
    for (node, chrom_id), chrom_rows in by_chromosome.items():
        seen.add(node)
        chrom = next(c for c in r.genomes[node].chromosomes if c.id == chrom_id)
        at = 0
        for row, block in zip(chrom_rows, chrom.blocks, strict=True):
            assert int(row["position"]) == at, f"node {node} chromosome {chrom_id} does not tile"
            assert (int(row["source"]), int(row["start"]), int(row["end"]), int(row["strand"])) == \
                   (block.source, block.start, block.end, block.strand)
            assert (int(row["copy"]), int(row["gene"])) == (block.copy, block.gene)
            at += block.end - block.start
        assert at == chrom.length
    assert seen == set(r.genomes), "ancestors as well as tips must be written"


def test_written_events_account_for_every_recorded_event(tmp_path):
    # an event can touch several blocks at once, so it writes several rows — but no event may go
    # unwritten, and no row may be invented
    r = _written(tmp_path)
    _, rows = _read(tmp_path / "genome_events.tsv")

    expected = collections.Counter()
    for e in r.events:
        kind = type(e).__name__.lower()
        n = {"loss": len(getattr(e, "lost", ())), "duplication": len(getattr(e, "copied", ())),
             "transfer": len(getattr(e, "transferred", ())),
             "speciation": len(getattr(e, "children", ()))}.get(kind, 1)
        expected[kind] += n
    for r_ in r.rearrangements:                      # they share the table; count them too
        expected[type(r_).__name__.lower()] += 1
    assert collections.Counter(row["kind"] for row in rows) == expected
    assert {"origination", "loss", "duplication", "transfer", "speciation"} <= set(expected), \
        f"the fixture should exercise every event kind, got {sorted(expected)}"


def test_written_genes_match_the_declared_spans(tmp_path):
    gff_dir = tmp_path / "in"
    gff_dir.mkdir()
    r = _written(tmp_path, gff=_gff(gff_dir))
    _, rows = _read(tmp_path / "genes.tsv")
    written = {int(row["family"]): (int(row["source"]), int(row["start"]), int(row["end"]))
               for row in rows}
    assert written == r.gene_spans
    named = {row["name"]: int(row["family"]) for row in rows if row["name"]}
    assert named == r.gene_names
    assert {int(row["family"]): int(row["strand"]) for row in rows} == r.gene_strands


def test_genes_file_is_header_only_when_none_were_declared(tmp_path):
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=3, seed=1)
    simulate_genomes_nucleotide(sp, root_length=200, seed=1).write(tmp_path, outputs=("genes",))
    assert (tmp_path / "genes.tsv").read_text() == "family\tname\tsource\tstart\tend\tstrand\n"


def test_written_rearrangements_share_the_event_table(tmp_path):
    """They used to be their own file. They end no gene lineage — which is why they were apart — but
    they are events on the same branches at the same clock, and a reader replaying one has to
    interleave them anyway. So: same table, their own kinds, and no ancestral columns."""
    r = _written(tmp_path, chromosomes=3, translocation=1.5, transposition=3.0)
    assert not (tmp_path / "rearrangements.tsv").exists()
    _, rows = _read(tmp_path / "genome_events.tsv")
    kinds = {"inversion", "transposition", "translocation"}
    rear = [row for row in rows if row["kind"] in kinds]
    assert len(rear) == len(r.rearrangements) and {row["kind"] for row in rear} == kinds
    for row, rec in zip(rear, sorted(r.rearrangements, key=lambda x: x.time), strict=True):
        assert float(row["time"]) == rec.time and node_from_label(row["lineage"]) == rec.lineage
        assert int(row["length"]) == rec.length and int(row["position"]) == rec.start
        # physical, not ancestral: a rearrangement moves material without changing where it came from
        assert row["source"] == row["start"] == row["end"] == row["copy"] == ""
    assert [float(row["time"]) for row in rows] == sorted(float(row["time"]) for row in rows)


def test_both_resolutions_write_the_same_chromosome_network_format(tmp_path):
    # one writer, one home: the network file must not drift between the two engines
    from zombi2.genomes import simulate_genomes_ordered

    _written(tmp_path / "nt")
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=6, seed=3)
    simulate_genomes_ordered(sp, initial_families=4, chromosomes=2, fission=0.2, seed=4).write(
        tmp_path / "ord", outputs=("chromosome_events",))

    nt_cols, _ = _read(tmp_path / "nt" / "chromosome_events.tsv")
    ord_cols, _ = _read(tmp_path / "ord" / "chromosome_events.tsv")
    assert nt_cols == ord_cols == ["time", "kind", "lineage", "parents", "children"]


def test_block_trees_cover_the_whole_genome_and_agree_with_the_gene_trees(tmp_path):
    # a block never splits, so its genealogy is in the log exactly as a gene's is: the same recovery
    # pointed at every block reconstructs the spacer too, which is what makes a whole ancestral
    # genome recoverable rather than a few declared loci.
    import re
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=6, seed=1)
    gff = tmp_path / "mini.gff"
    gff.write_text("##gff-version 3\n##sequence-region c 1 3000\n"
                   "c\tt\tgene\t1\t300\t.\t+\t.\tID=a\n"
                   "c\tt\tgene\t601\t900\t.\t+\t.\tID=b\n"
                   "c\tt\tgene\t1201\t1500\t.\t-\t.\tID=c\n")
    g = simulate_genomes_nucleotide(sp, gff=gff, duplication=0.2, loss=0.2, inversion=0.3, seed=1)

    assert len(g.block_trees) == len(g.root_blocks)          # every block, not just the genes
    assert len(g.gene_trees) < len(g.block_trees)            # the spacer is the difference

    # a declared gene recovers the same genealogy either way. The g<id> labels differ — segment ids
    # are handed out as the recovery walks its targets — so compare the shape, which is the claim.
    shape = lambda nwk: re.sub(r"g\d+", "g", nwk)            # noqa: E731
    span_to_gene = {span: fam for fam, span in g.gene_spans.items()}
    checked = 0
    for i, interval in enumerate(g.root_blocks):
        if interval in span_to_gene:
            a = g.block_trees[i].to_newick("complete")
            b = g.gene_trees[span_to_gene[interval]].to_newick("complete")
            assert shape(a) == shape(b)
            checked += 1
    assert checked == len(g.gene_trees)


def _expand(result, node_id):
    """The assembly spelled back out as ``{chromosome: [(source, position, strand), …]}`` — one entry
    per nucleotide, the same shape ``trace_back`` returns. A ``-1`` piece is read down its source."""
    blocks = result.root_blocks
    out = {}
    for cid, pieces in result.assembly(node_id).items():
        seq = []
        for (i, _gene, strand) in pieces:
            src, a, z = blocks[i]
            span = range(a, z)
            seq.extend((src, p, strand) for p in (span if strand == 1 else reversed(span)))
        out[cid] = seq
    return out


def test_assembly_tiles_every_node_exactly_as_its_trace_back():
    # assembly() says which recovered block each piece of a genome is; expanding that back to one entry
    # per nucleotide has to give the per-nucleotide ancestry the genome already records. Nucleotide for
    # nucleotide, at EVERY node — this is what a reconstructed genome rests on.
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=4)
    g = simulate_genomes_nucleotide(sp, inversion=3.0, inversion_length=80, loss=0.4, loss_length=40,
                                    duplication=0.4, duplication_length=40, transfer=0.6,
                                    transfer_length=60, root_length=600, genes=3, gene_length=90,
                                    seed=4)
    for node_id in sorted(g.genomes):                    # no node is skipped and none refuses
        assert _expand(g, node_id) == g.trace_back(node_id)
    assert len(g.genomes) > len(g.complete_tree.extant()), "the tree has no ancestors to check"


def test_the_partition_is_at_least_as_fine_as_every_nodes_blocks():
    # the invariant the whole reconstruction rests on: because every node votes on where the cuts go,
    # each of its blocks is a whole number of root blocks — never a fragment of one. That is what makes
    # a piece always a whole block, and what a partition cut from the survivors alone cannot give.
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=4)
    g = simulate_genomes_nucleotide(sp, inversion=3.0, inversion_length=80, loss=1.5, loss_length=60,
                                    duplication=0.6, duplication_length=40, transfer=0.6,
                                    transfer_length=60, root_length=600, genes=3, gene_length=90,
                                    seed=4)
    cuts = collections.defaultdict(set)
    for src, a, b in g.root_blocks:
        cuts[src].update((a, b))
    spans = 0
    for node_id, genome in g.genomes.items():
        for chrom in genome.chromosomes:
            for blk in chrom.blocks:
                assert blk.start in cuts[blk.source] and blk.end in cuts[blk.source]
                spans += sum(1 for (s, a, b) in g.root_blocks
                             if s == blk.source and blk.start <= a and b <= blk.end) > 1
    assert spans, "no block spanned several root blocks — the multi-piece path went untested"


# --- GFF / BED export ------------------------------------------------------------------------------

def _export_run(tmp_path, **kw):
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=4, seed=5)
    params = dict(root_length=1200, genes=4, gene_length=150, inversion=3.0, inversion_length=200,
                  duplication=1.5, loss=0.5, seed=5)
    params.update(kw)
    g = simulate_genomes_nucleotide(sp, **params)
    g.write(tmp_path)
    return g


def test_gff_and_bed_are_written_for_every_genome_and_name_the_fasta_records(tmp_path):
    # the whole point of the seqid: a genome and its annotation have to join without renaming
    from zombi2.sequences import simulate_sequences
    from zombi2.sequences.substitution_models import jc69

    g = _export_run(tmp_path)
    simulate_sequences(g, model=jc69(), substitution=0.05, seed=5).write(tmp_path)
    labels = [node_label(i) for i in g.genomes] + ["initial"]
    assert {p.stem for p in tmp_path.glob("genome_*.gff")} == {f"genome_{lab}" for lab in labels}
    assert {p.stem for p in tmp_path.glob("genome_*.bed")} == {f"genome_{lab}" for lab in labels}
    for lab in labels:
        fasta = {ln[1:] for ln in (tmp_path / f"genome_{lab}.fasta").read_text().splitlines()
                 if ln.startswith(">")}
        bed = {ln.split("\t")[0] for ln in (tmp_path / f"genome_{lab}.bed").read_text().splitlines()}
        gff = {ln.split("\t")[0] for ln in (tmp_path / f"genome_{lab}.gff").read_text().splitlines()
               if not ln.startswith("#")}
        assert bed == fasta and gff <= fasta


def test_bed_tiles_the_genome_and_names_each_block_by_its_ancestry(tmp_path):
    g = _export_run(tmp_path)
    for node_id, genome in g.genomes.items():
        rows = [ln.split("\t") for ln in
                (tmp_path / f"genome_{node_label(node_id)}.bed").read_text().splitlines()]
        assert len(rows) == sum(len(c.blocks) for c in genome.chromosomes)
        at = collections.defaultdict(int)
        for (chrom, start, end, name, _score, strand) in rows:
            assert int(start) == at[chrom]                 # 0-based half-open, tiling from 0
            at[chrom] = int(end)
            assert strand in ("+", "-")
        for chrom in genome.chromosomes:                   # ...to exactly the chromosome's length
            assert at[f"{node_label(node_id)}_chr{chrom.id}"] == chrom.length
        expected = {f"{b.source}:{b.start}-{b.end}" for c in genome.chromosomes for b in c.blocks}
        assert {r[3] for r in rows} == expected


def test_gff_gives_every_gene_unique_id_right_coordinates_and_the_strand_it_now_reads_on(tmp_path):
    g = _export_run(tmp_path)
    flipped = 0
    for node_id, genome in g.genomes.items():
        rows = [ln.split("\t") for ln in
                (tmp_path / f"genome_{node_label(node_id)}.gff").read_text().splitlines()
                if not ln.startswith("#")]
        genes = [(c, at, b) for c in genome.chromosomes
                 for at, b in [(sum(x.length for x in c.blocks[:i]), c.blocks[i])
                               for i in range(len(c.blocks))] if b.is_gene]
        assert len(rows) == len(genes)
        assert len({r[8].split(";")[0] for r in rows}) == len(rows)          # IDs unique in the file
        for (chrom, at, b), row in zip(genes, rows):
            assert (int(row[3]), int(row[4])) == (at + 1, at + b.length)     # GFF is 1-based inclusive
            attrs = dict(kv.split("=", 1) for kv in row[8].split(";"))
            assert int(attrs["family"]) == b.gene and int(attrs["copy"]) == b.copy
            assert attrs["source"] == f"{b.source}:{b.start}-{b.end}"
            # the strand it reads on HERE: coding strand, flipped if the block has been inverted
            expected = "+" if g.gene_strands.get(b.gene, 1) * b.strand == 1 else "-"
            assert row[6] == expected
            flipped += b.strand == -1
    assert flipped, "no gene sits on an inverted block — the flip went untested"


def test_a_gff_we_wrote_reads_back_through_our_own_gff_reader(tmp_path):
    # the strongest thing to check about a format: our reader accepts what our writer produced, and
    # puts the genes back where they were
    from zombi2.genomes.gff import read_gff

    g = _export_run(tmp_path, inversion=0.0, duplication=0.0, loss=0.0)   # nothing moved yet
    lengths, genes = read_gff(tmp_path / "genome_initial.gff")
    chrom = g.initial_genome.chromosomes[0]
    assert lengths == {f"initial_chr{chrom.id}": chrom.length}
    declared = sorted((a, b) for (_src, a, b) in g.gene_spans.values())
    assert sorted((gene.start, gene.end) for gene in genes) == declared
