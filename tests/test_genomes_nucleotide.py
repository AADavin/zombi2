"""Tests for nucleotide genomes — the coordinate representation, inversion, the multi-chromosome,
identity-bearing karyotype wired along the species tree, and the copy-lineage genealogy log.

A chromosome is an ordered list of blocks (runs of one ancestry — not merged during the run: option
B, maximality is a recovered property). The strongest checks are the oracle (apply the same
inversions to a plain per-nucleotide array and require ``trace_back`` to agree), the strong invariant
(an ancestry-neutral event conserves ancestry, so every node carries the whole root sequence,
permuted), and — with loss / duplication — that the copy-lineage log accounts for every node's copy
numbers exactly and every copy traces back to a seed origination.
"""

import numpy as np
import pytest

from zombi2.species import simulate_species_tree
from zombi2.genomes import simulate_genomes_nucleotide
from zombi2.genomes.nucleotide import (
    Block,
    Chromosome,
    Duplication,
    Loss,
    NucleotideGenome,
    Origination,
    Speciation,
    Transfer,
    Translocation,
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

def _run(seed=1, inversion=0.04, chromosomes=1, **kw):
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
    r = simulate_genomes_nucleotide(sp, inversion=0.0, chromosomes=specs, seed=5)
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
    r = simulate_genomes_nucleotide(sp, inversion=0.06, inversion_length=10, chromosomes=specs, seed=3)
    full = sorted((s, p) for s, (length, _t) in enumerate(specs) for p in range(length))
    for node_id in r.genomes:
        assert r.ancestry(node_id) == full


def test_chromosome_identity_network_is_well_formed():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=4)
    r = simulate_genomes_nucleotide(sp, inversion=0.02,
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
    sp, r = _run(seed=6, inversion=0.0, chromosomes=2, root_length=50)
    assert r.rearrangements == []
    for g in r.genomes.values():
        assert all(len(c.blocks) == 1 for c in g.chromosomes)


def test_inversions_recorded_within_their_branch_and_chromosome():
    sp, r = _run(seed=7, inversion=0.08, chromosomes=2, root_length=100)
    assert r.rearrangements
    node_chroms = {n: {c.id for c in r.genomes[n].chromosomes} for n in r.genomes}
    for inv in r.rearrangements:
        node = sp.complete_tree.nodes[inv.lineage]
        assert node.birth_time <= inv.time <= node.end_time
        assert inv.chromosome in node_chroms[inv.lineage]    # fired on a chromosome that node carried
        assert inv.length >= 1


def test_deterministic_given_seed():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=8)
    kw = dict(inversion=0.05, inversion_length=15, chromosomes=[(120, "circular"), (30, "linear")], seed=8)
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
    r = simulate_genomes_nucleotide(sp, duplication=0.04, duplication_length=10, inversion=0.02,
                                    chromosomes=specs, seed=3)
    assert any(isinstance(e, Duplication) for e in r.events)
    for node_id in r.genomes:
        assert set(r.ancestry(node_id)) == full             # nothing lost
    assert any(len(r.ancestry(n)) > len(full) for n in r.genomes)   # copies really exist


def test_duplication_with_loss_stays_a_subset():
    specs = [(120, "circular"), (40, "linear")]
    full = {(s, p) for s, (length, _t) in enumerate(specs) for p in range(length)}
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=4)
    r = simulate_genomes_nucleotide(sp, duplication=0.04, loss=0.04, inversion=0.02, translocation=0.03,
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
    r = simulate_genomes_nucleotide(sp, loss=0.05, loss_length=10, inversion=0.02,
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
        r = simulate_genomes_nucleotide(sp, loss=0.04, inversion=0.02, translocation=0.03,
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
    params = dict(inversion=0.02, inversion_length=12, fission=0.6, fusion=0.6,
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
    kw = dict(inversion=0.02, fission=0.5, fusion=0.5, chromosomes=_TIER_SPECS, seed=8)
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
    r = simulate_genomes_nucleotide(sp, translocation=0.06, translocation_length=10,
                                    chromosomes=_TIER_SPECS, seed=5)
    assert all(len(g.chromosomes) == 3 for g in r.genomes.values())   # material moves, count does not
    assert any(isinstance(x, Translocation) for x in r.rearrangements)


def test_strong_invariant_survives_translocation():
    full = sorted((s, p) for s, (length, _t) in enumerate(_TIER_SPECS) for p in range(length))
    for seed in range(3):
        sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=seed)
        r = simulate_genomes_nucleotide(sp, translocation=0.08, translocation_length=8, inversion=0.02,
                                        inversion_probability=0.5, chromosomes=_TIER_SPECS, seed=seed)
        for node_id in r.genomes:
            assert r.ancestry(node_id) == full


def test_all_four_events_compose():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=5)
    r = simulate_genomes_nucleotide(sp, inversion=0.02, translocation=0.05, fission=0.4, fusion=0.4,
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
    r = simulate_genomes_nucleotide(sp, duplication=0.05, loss=0.05, duplication_length=8,
                                    loss_length=8, inversion=0.03, translocation=0.03, fission=0.15,
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
    r = simulate_genomes_nucleotide(sp, duplication=0.05, loss=0.03, inversion=0.03,
                                    translocation=0.03, fission=0.15, fusion=0.15, chromosomes=specs,
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
    r = simulate_genomes_nucleotide(sp, inversion=0.03, duplication=0.03, chromosomes=specs, seed=2)
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
    r = simulate_genomes_nucleotide(sp, inversion=0.05, inversion_length=12, translocation=0.03,
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
        r = simulate_genomes_nucleotide(sp, duplication=0.05, loss=0.05, duplication_length=8,
                                        loss_length=8, inversion=0.03, translocation=0.03, fission=0.1,
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
    r = simulate_genomes_nucleotide(sp, duplication=0.05, loss=0.04, inversion=0.03, chromosomes=specs,
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
    r = simulate_genomes_nucleotide(sp, duplication=0.06, duplication_length=10,
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
    r = simulate_genomes_nucleotide(sp, transfer=0.07, transfer_length=10, inversion=0.02,
                                    chromosomes=_XFER_SPECS, seed=2)
    assert any(isinstance(e, Transfer) for e in r.events)          # transfers really fired
    for node_id in r.genomes:
        assert set(r.ancestry(node_id)) == _XFER_FULL              # nothing lost
    assert any(len(r.ancestry(n)) > len(_XFER_FULL) for n in r.genomes)   # extra copies exist


def test_transfer_recipient_is_a_distinct_contemporaneous_lineage():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=4)
    r = simulate_genomes_nucleotide(sp, transfer=0.08, transfer_length=8, chromosomes=_XFER_SPECS, seed=4)
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
        r = simulate_genomes_nucleotide(sp, transfer=0.05, transfer_length=8, duplication=0.03,
                                        loss=0.03, inversion=0.03, translocation=0.02, fission=0.1,
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
    r = simulate_genomes_nucleotide(sp, transfer=0.08, transfer_length=10, chromosomes=_XFER_SPECS, seed=3)

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
        r = simulate_genomes_nucleotide(sp, transfer=0.06, transfer_length=8, chromosomes=_XFER_SPECS,
                                        seed=5, **kw)
        assert any(isinstance(e, Transfer) for e in r.events)
        for node_id in r.genomes:                              # additive: full ancestry preserved
            assert set(r.ancestry(node_id)) == _XFER_FULL


def test_transfer_validation():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=4, seed=1)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, transfer=-0.1)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, transfer=0.1, transfer_length=0)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, transfer=0.1, transfer_to="nearest")
