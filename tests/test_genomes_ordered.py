"""Tests for ordered genomes — genes with a position and a strand, on chromosomes (zombi2.genomes).

Slice 1: the chromosome container + D/T/L/O made position-aware + inversions + chromosome identity.
The gene genealogy is the shared, position-blind event log, so the strong gene-tree/profile invariant
must survive unchanged; on top of it we check the two new things — the inversion operator and the
chromosome genealogy.
"""

import inspect

import pytest

from zombi2.rates import scope
from zombi2.rates.distributions import Geometric
from zombi2.rates.modifiers import OnTime
from zombi2.species import Node, Tree, simulate_species_tree
from zombi2.genomes import (
    Chromosome,
    Gene,
    Inversion,
    Transposition,
    Translocation,
    simulate_genomes_ordered,
    simulate_genomes_unordered,
)
from zombi2.genomes.ordered import _duplicate, _fission, _fusion, _invert, _transpose, _translocate


def _tier(seed, death=0.5, n_extant=15, **kw):
    sp = simulate_species_tree(birth=1.0, death=death, n_extant=n_extant, seed=seed)
    params = dict(duplication=0.4, transfer=0.3, loss=0.2, origination=0.6, inversion=0.3,
                  chromosomes=3, fission=0.15, fusion=0.15, chromosome_origination=0.08,
                  chromosome_loss=0.12, initial_families=14, seed=seed)
    params.update(kw)
    return sp, simulate_genomes_ordered(sp, **params)


def _run(seed=1, death=0.4, n_extant=15, **kw):
    sp = simulate_species_tree(birth=1.0, death=death, n_extant=n_extant, seed=seed)
    params = dict(duplication=0.3, transfer=0.2, loss=0.25, origination=0.6, inversion=0.4,
                  chromosomes=4, initial_families=16, seed=seed)
    params.update(kw)
    return sp, simulate_genomes_ordered(sp, **params)


def _extant_leaves(node):
    if node is None:
        return 0
    if node.is_leaf:
        return 1 if node.kind == "extant" else 0
    return sum(_extant_leaves(c) for c in node.children)


# --- the genome is chromosomes of oriented genes -------------------------------------------------

def test_genome_is_a_tuple_of_chromosomes_of_oriented_genes():
    _, r = _run(seed=2)
    for chroms in r.genomes.values():
        assert isinstance(chroms, tuple)
        for ch in chroms:
            assert isinstance(ch, Chromosome) and ch.topology in ("circular", "linear")
            assert all(isinstance(g, Gene) and g.strand in (1, -1) for g in ch.genes)


def test_seeded_chromosome_count_and_topology():
    nodes = {0: Node(0, None, 0.0, 1.0, None, "extant")}  # a lone leaf: its genome is the seed
    r = simulate_genomes_ordered(Tree(nodes, 0), chromosomes=5, topology="linear",
                                 initial_families=0, seed=1)
    assert len(r.genomes[0]) == 5
    assert all(ch.topology == "linear" for ch in r.genomes[0])


def test_initial_families_dealt_round_robin_across_chromosomes():
    nodes = {0: Node(0, None, 0.0, 1.0, None, "extant")}  # no events -> genome is exactly the seed
    r = simulate_genomes_ordered(Tree(nodes, 0), chromosomes=3, initial_families=7, seed=1)
    # 7 genes over 3 chromosomes, round-robin: 3, 2, 2
    assert [len(ch.genes) for ch in r.genomes[0]] == [3, 2, 2]


def test_shared_params_are_a_subset_of_the_ordered_signature():
    # the layering contract (genome-api.md): unordered ⊂ ordered, only the extra args differ
    shared = set(inspect.signature(simulate_genomes_unordered).parameters) - {"tree"}
    ordered = set(inspect.signature(simulate_genomes_ordered).parameters) - {"tree"}
    assert shared <= ordered                                 # unordered ⊂ ordered: nothing dropped
    assert ordered - shared == {                            # ordered's own additions:
        "chromosomes", "topology", "inversion", "transposition", "translocation",
        "fission", "fusion", "chromosome_origination", "chromosome_loss", "inversion_probability",
        "duplication_extension", "loss_extension", "transfer_extension",
        "inversion_extension", "transposition_extension", "translocation_extension"}


# --- the shared gene genealogy still holds -------------------------------------------------------

def test_extant_gene_tree_leaves_equal_the_extant_copy_total():
    # the strongest invariant, inherited from the unordered core: surviving gene-tree leaves == copies
    sp, r = _run(seed=5, death=0.5)
    extant_sp = {n.id for n in sp.complete_tree.extant()}
    for fam, tree in r.gene_trees.items():
        copies = sum(r.profiles.counts.get((fam, s), 0) for s in extant_sp)
        assert _extant_leaves(tree.extant) == copies


def test_per_node_gene_ids_are_unique():
    _, r = _run(seed=6)
    for chroms in r.genomes.values():
        ids = [g.id for ch in chroms for g in ch.genes]
        assert len(ids) == len(set(ids))


def test_family_counts_and_gene_order_agree():
    _, r = _run(seed=7)
    for node_id in r.genomes:
        order = r.gene_order(node_id)
        assert len(order) == sum(r.family_counts(node_id).values())
        # gene_order lists (chromosome, position, strand, family, gid); positions run 0..len-1 per chrom
        for ch in r.genomes[node_id]:
            rows = [row for row in order if row[0] == ch.id]
            assert [row[1] for row in rows] == list(range(len(ch.genes)))


# --- inversions ----------------------------------------------------------------------------------

def test_invert_reverses_the_span_and_flips_each_strand():
    ch = Chromosome(0, "linear", [Gene(0, 0, 1), Gene(1, 1, 1), Gene(2, 2, -1), Gene(3, 3, 1)])
    node = Node(7, None, 0.0, 1.0, None, "extant")
    rearr = []
    _invert(ch, 1, 2, node, 3.0, rearr)
    assert [g.id for g in ch.genes] == [0, 2, 1, 3]        # the span reversed, ids preserved
    assert [g.strand for g in ch.genes] == [1, 1, -1, 1]   # the two inverted genes flipped strand
    assert rearr == [Inversion(3.0, 7, 0, 1, 2)] and rearr[0].length == 2


def test_no_inversions_when_the_rate_is_zero():
    _, r = _run(seed=3, inversion=0.0)
    assert r.rearrangements == []


def test_inversions_never_remint_gene_ids():
    # a single branch with only inversions: the genes at the tip are exactly the seeded ids
    nodes = {0: Node(0, None, 0.0, 1.0, None, "extant")}
    r = simulate_genomes_ordered(Tree(nodes, 0), inversion=6.0, chromosomes=1,
                                 initial_families=8, seed=2)
    assert r.rearrangements                                     # inversions really fired
    assert {g.id for ch in r.genomes[0] for g in ch.genes} == set(range(8))


def test_recorded_inversions_are_well_formed():
    _, r = _run(seed=4)
    assert r.rearrangements
    for inv in r.rearrangements:
        assert 0 <= inv.start <= inv.end
        assert inv.length == inv.end - inv.start + 1


# --- the chromosome genealogy --------------------------------------------------------------------

def test_seed_chromosomes_are_origination_roots():
    _, r = _run(seed=8, chromosomes=4)
    roots = [e for e in r.chromosome_events if e.kind == "origination"]
    assert len(roots) == 4
    assert all(e.parents == () and len(e.children) == 1 for e in roots)


def test_speciation_chromosome_edges_are_one_parent_to_two_children():
    _, r = _run(seed=9)
    spec = [e for e in r.chromosome_events if e.kind == "speciation"]
    assert spec
    assert all(len(e.parents) == 1 and len(e.children) == 2 for e in spec)


def test_chromosome_genealogy_is_a_connected_forest():
    # every chromosome that is a speciation parent was itself born by an earlier event (no orphans),
    # and no chromosome id is ever produced twice (identity is re-minted, never reused)
    _, r = _run(seed=10)
    born = {}
    for e in r.chromosome_events:
        for ch in e.children:
            assert ch not in born                              # each id minted exactly once
            born[ch] = e
    for e in r.chromosome_events:
        for p in e.parents:
            assert p in born                                   # every parent was born earlier


def test_chromosome_count_is_conserved_through_speciation():
    # slice 1 has no fission/fusion, so a daughter has exactly its parent's chromosome count
    _, r = _run(seed=11, chromosomes=4)
    for node_id, chroms in r.genomes.items():
        node = r.complete_tree.nodes[node_id]
        if node.children is None:
            assert len(chroms) == 4                            # inherited unchanged down every branch


# --- inherited mechanics (transfer) --------------------------------------------------------------

def test_transfer_events_appear_and_cross_species_branches():
    _, r = _run(seed=7, self_transfer=False)
    xfer_rows = [e for e in r.events if e.kind == "transfer"]
    assert xfer_rows                                           # transfers really fired
    # the recipient row carries a recipient different from the donor lineage
    recip = [e for e in xfer_rows if e.recipient is not None]
    assert recip and all(e.lineage == e.recipient for e in recip)


def test_no_transfer_events_when_transfer_is_zero():
    _, r = _run(seed=8, transfer=0.0)
    assert not any(e.kind == "transfer" for e in r.events)


def test_replacement_run_stays_consistent():
    # replacement overwrites a homologous copy; the strong invariant must still hold
    sp, r = _run(seed=2, replacement=True)
    extant_sp = {n.id for n in sp.complete_tree.extant()}
    for fam, tree in r.gene_trees.items():
        assert _extant_leaves(tree.extant) == sum(r.profiles.counts.get((fam, s), 0) for s in extant_sp)


# --- determinism, validation, writing ------------------------------------------------------------

def test_deterministic_given_seed():
    sp, r = _run(seed=3)
    r2 = simulate_genomes_ordered(sp, duplication=0.3, transfer=0.2, loss=0.25, origination=0.6,
                                  inversion=0.4, chromosomes=4, initial_families=16, seed=3)
    assert all(r.gene_order(x) == r2.gene_order(x) for x in r.genomes)
    assert r.rearrangements == r2.rearrangements
    assert r.chromosome_events == r2.chromosome_events


def test_ontime_skyline_modifier_is_accepted():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=1)
    r = simulate_genomes_ordered(sp, duplication=0.3 * OnTime({0: 1.0, 1.0: 0.2}),
                                 inversion=0.2 * OnTime({0: 0.5, 1.0: 2.0}),
                                 chromosomes=2, initial_families=6, seed=1)
    assert r.genomes                                           # ran without complaint


def test_scope_override_is_rejected_this_slice():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=1)
    with pytest.raises(ValueError, match="scope"):
        simulate_genomes_ordered(sp, duplication=scope.Global(0.3), chromosomes=2, seed=1)
    with pytest.raises(ValueError, match="scope"):
        simulate_genomes_ordered(sp, inversion=scope.PerLineage(0.3), chromosomes=2, seed=1)


def test_topology_validation():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=1)
    with pytest.raises(ValueError, match="topology"):
        simulate_genomes_ordered(sp, chromosomes=3, topology=["circular", "linear"], seed=1)  # wrong length
    with pytest.raises(ValueError, match="topology"):
        simulate_genomes_ordered(sp, chromosomes=2, topology="loop", seed=1)                  # bad label
    with pytest.raises(ValueError, match="chromosomes"):
        simulate_genomes_ordered(sp, chromosomes=0, seed=1)


def test_mixed_topology_per_chromosome():
    nodes = {0: Node(0, None, 0.0, 1.0, None, "extant")}
    r = simulate_genomes_ordered(Tree(nodes, 0), chromosomes=2, topology=["circular", "linear"],
                                 initial_families=4, seed=1)
    assert [ch.topology for ch in r.genomes[0]] == ["circular", "linear"]


def test_write_emits_the_selected_outputs(tmp_path):
    _, r = _run(seed=5)
    r.write(tmp_path, outputs=("events", "profiles", "gene_order", "rearrangements",
                               "chromosome_events"))
    for name in ("genome_events.tsv", "profiles.tsv", "gene_order.tsv", "rearrangements.tsv",
                 "chromosome_events.tsv"):
        assert (tmp_path / name).exists()
    head = (tmp_path / "gene_order.tsv").read_text().splitlines()[0]
    assert head.split("\t") == ["species", "chromosome", "position", "strand", "family", "gene"]


def _written_gene_order(path):
    """``gene_order.tsv`` -> ``{node: [(chromosome, position, strand, family, gene), ...]}``."""
    lines = (path / "gene_order.tsv").read_text().splitlines()[1:]
    written = {}
    for row in lines:
        s, *rest = (int(c) for c in row.split("\t"))
        written.setdefault(s, []).append(tuple(rest))
    return written


def test_gene_order_is_written_for_every_node_not_only_the_tips(tmp_path):
    # ancestral layouts are what make rearrangements.tsv replayable: an inversion's (start, length)
    # only means something against the genome its branch started from — its parent's rows.
    _, r = _run(seed=5)
    r.write(tmp_path, outputs=("gene_order",))
    written = _written_gene_order(tmp_path)

    internal = {n.id for n in r.complete_tree.nodes.values() if n.children is not None}
    assert internal, "the fixture tree should have internal nodes to write"
    # every node with genes is present — root and internal branches included, not just the tips
    assert set(written) == {s for s in r.genomes if r.gene_order(s)}
    assert internal & set(written)
    # and each node's written rows are that node's actual layout
    for s, rows in written.items():
        assert rows == r.gene_order(s)


def test_empty_run_has_chromosomes_but_no_genes():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=1)
    r = simulate_genomes_ordered(sp, chromosomes=3, seed=1)   # no families, no events
    assert r.events == [] and r.rearrangements == []
    assert all(len(chroms) == 3 and sum(len(ch.genes) for ch in chroms) == 0
               for chroms in r.genomes.values())
    assert r.gene_trees == {}


# --- slice 2: the chromosome tier (the reticulating network) -------------------------------------

def test_tier_events_fire_with_the_right_network_arity():
    _, r = _tier(seed=1)
    seen = {e.kind for e in r.chromosome_events}
    assert {"origination", "speciation", "fission", "fusion", "loss"} <= seen  # all five kinds occur
    shape = {"origination": (0, 1), "speciation": (1, 2), "fission": (1, 2),
             "fusion": (2, 1), "loss": (1, 0)}
    for e in r.chromosome_events:
        assert (len(e.parents), len(e.children)) == shape[e.kind]


def test_fusion_is_the_only_reticulation():
    # only a fusion has in-degree 2; every other event has exactly one parent (or none)
    _, r = _tier(seed=2)
    for e in r.chromosome_events:
        assert (len(e.parents) == 2) == (e.kind == "fusion")


def test_chromosome_network_is_a_connected_dag():
    # every chromosome id is minted exactly once (identity re-minted, never reused), and every parent
    # of every event was itself born by an earlier event — the fusion parents included
    _, r = _tier(seed=3)
    born = {}
    for e in r.chromosome_events:
        for ch in e.children:
            assert ch not in born
            born[ch] = e
    assert all(p in born for e in r.chromosome_events for p in e.parents)


def test_strong_invariant_survives_the_tier():
    # chromosome loss ends its genes as gene losses; they must become gene-tree death, not extant
    # leaves — so the surviving-leaves == profile-copies invariant must still hold under heavy tier
    for seed in range(5):
        sp, r = _tier(seed=seed)
        extant = {n.id for n in sp.complete_tree.extant()}
        for fam, tree in r.gene_trees.items():
            assert _extant_leaves(tree.extant) == sum(r.profiles.counts.get((fam, s), 0) for s in extant)


def test_a_genome_never_loses_its_last_chromosome():
    _, r = _tier(seed=4, chromosome_loss=1.0, chromosome_origination=0.0)  # push loss hard
    assert all(len(chroms) >= 1 for chroms in r.genomes.values())


def test_the_tier_changes_chromosome_number():
    _, r = _tier(seed=5)
    counts = {len(chroms) for chroms in r.genomes.values()}
    assert len(counts) > 1                                    # not the conserved single value of slice 1


def test_tier_rates_zero_is_byte_identical_to_a_no_tier_call():
    # the added firing branches must not perturb the RNG stream when the tier is off
    sp = simulate_species_tree(birth=1.0, death=0.4, n_extant=15, seed=7)
    base = dict(duplication=0.3, transfer=0.15, loss=0.25, origination=0.5, inversion=0.3,
                chromosomes=3, initial_families=12, seed=7)
    a = simulate_genomes_ordered(sp, **base)
    b = simulate_genomes_ordered(sp, **base, fission=0.0, fusion=0.0, chromosome_origination=0.0,
                                 chromosome_loss=0.0)
    assert all(a.gene_order(x) == b.gene_order(x) for x in a.genomes)
    assert a.chromosome_events == b.chromosome_events and a.rearrangements == b.rearrangements


def test_de_novo_replicon_is_an_empty_origination_root():
    # only chromosome origination on a lone branch: each de-novo replicon is a rootless-parent, empty
    nodes = {0: Node(0, None, 0.0, 5.0, None, "extant")}
    r = simulate_genomes_ordered(Tree(nodes, 0), chromosome_origination=1.0, chromosomes=1,
                                 initial_families=3, seed=1)
    de_novo = [e for e in r.chromosome_events if e.kind == "origination" and e.parents == ()][1:]
    assert de_novo                                           # some replicons appeared past the seed
    # the seed carried the 3 genes; the de-novo replicons are empty
    assert sum(len(ch.genes) for ch in r.genomes[0]) == 3
    assert len(r.genomes[0]) == 1 + len(de_novo)


def _minter(start):
    box = [start]

    def mint():
        box[0] += 1
        return box[0]
    return mint


def test_fission_partitions_genes_in_order_preserving_ids():
    import numpy as np
    genome = [Chromosome(0, "linear", [Gene(0, 0, 1), Gene(1, 1, 1), Gene(2, 2, -1), Gene(3, 3, 1)])]
    node = Node(5, None, 0.0, 1.0, None, "extant")
    ce = []
    dc, dg = _fission(genome, 0, node, 2.0, ce, _minter(99), np.random.default_rng(0))
    assert dc == 1 and dg == 0
    assert len(genome) == 2
    assert [g.id for ch in genome for g in ch.genes] == [0, 1, 2, 3]   # order preserved across split
    assert ce[0].kind == "fission" and ce[0].parents == (0,) and len(set(ce[0].children)) == 2


def test_fusion_concatenates_two_chromosomes_into_one():
    import numpy as np
    genome = [Chromosome(10, "circular", [Gene(0, 0, 1), Gene(1, 1, 1)]),
              Chromosome(11, "circular", [Gene(2, 2, 1)])]
    node = Node(5, None, 0.0, 1.0, None, "extant")
    ce = []
    dc, dg = _fusion(genome, 0, node, 3.0, ce, _minter(20), np.random.default_rng(0))
    assert dc == -1 and dg == 0
    assert len(genome) == 1
    assert [g.id for g in genome[0].genes] == [0, 1, 2]              # a.genes + b.genes
    assert ce[0].kind == "fusion" and ce[0].parents == (10, 11) and len(ce[0].children) == 1


def test_tier_rate_scope_override_is_rejected():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=1)
    with pytest.raises(ValueError, match="scope"):
        simulate_genomes_ordered(sp, fission=scope.Global(0.1), chromosomes=2, seed=1)


# --- slice 3: segmental events (the extension) + transposition + translocation --------------------

def test_geometric_mean_one_is_always_a_single_gene():
    import numpy as np
    rng = np.random.default_rng(0)
    assert all(Geometric(mean=1).sample(rng) == 1.0 for _ in range(100))
    with pytest.raises(ValueError):
        Geometric(mean=0.5)


def test_duplicate_copies_a_block_in_tandem():
    ch = Chromosome(0, "linear", [Gene(0, 0, 1), Gene(1, 1, 1), Gene(2, 2, 1)])
    events, positions, counter = [], [], [10]

    def ng(fam, strand):
        counter[0] += 1
        return Gene(counter[0], fam, strand)
    added = _duplicate(ch, 0, 2, Node(3, None, 0.0, 1.0, None, "extant"), 1.0, events, positions, ng)
    assert added == 2 and len(ch.genes) == 5
    assert [g.family for g in ch.genes] == [0, 1, 0, 1, 2]   # conts in place, then the tandem copy block
    assert len(events) == 4 and all(e.kind == "duplication" for e in events)
    # one position row for the whole event, naming where the copy block landed
    assert len(positions) == 1
    p = positions[0]
    assert (p.kind, p.chromosome, p.start, p.length, p.dest_position) == ("duplication", 0, 0, 2, 2)


def test_transpose_relocates_a_segment_within_the_chromosome_preserving_ids():
    import numpy as np
    ch = Chromosome(0, "linear", [Gene(i, i, 1) for i in range(5)])
    rearr = []
    _transpose(ch, 0, 2, Node(3, None, 0.0, 1.0, None, "extant"), 1.0, rearr,
               np.random.default_rng(0), 0.0)
    assert sorted(g.id for g in ch.genes) == [0, 1, 2, 3, 4]  # same genes, reordered — nothing lost
    assert len(ch.genes) == 5
    assert isinstance(rearr[0], Transposition) and rearr[0].length == 2 and rearr[0].flipped is False


def test_transpose_flips_the_segment_when_inversion_probability_is_one():
    import numpy as np
    ch = Chromosome(0, "linear", [Gene(0, 0, 1), Gene(1, 1, 1), Gene(2, 2, 1)])
    rearr = []
    _transpose(ch, 0, 2, Node(3, None, 0.0, 1.0, None, "extant"), 1.0, rearr,
               np.random.default_rng(0), 1.0)
    strands = {g.id: g.strand for g in ch.genes}
    assert rearr[0].flipped is True
    assert strands[0] == -1 and strands[1] == -1 and strands[2] == 1  # only the moved block flipped


def test_translocate_moves_a_segment_to_a_different_chromosome():
    import numpy as np
    genome = [Chromosome(0, "linear", [Gene(0, 0, 1), Gene(1, 1, 1), Gene(2, 2, 1)]),
              Chromosome(1, "linear", [Gene(3, 3, 1)])]
    rearr = []
    _translocate(genome, 0, 0, 2, Node(5, None, 0.0, 1.0, None, "extant"), 1.0, rearr,
                 np.random.default_rng(0), 0.0)
    assert sorted(g.id for ch in genome for g in ch.genes) == [0, 1, 2, 3]   # nothing gained/lost
    assert len(genome[0].genes) == 1 and {g.id for g in genome[1].genes} == {0, 1, 3}
    assert isinstance(rearr[0], Translocation) and rearr[0].source == 0 and rearr[0].dest == 1


def test_translocate_is_a_noop_with_a_single_chromosome():
    import numpy as np
    genome = [Chromosome(0, "linear", [Gene(0, 0, 1), Gene(1, 1, 1)])]
    rearr = []
    _translocate(genome, 0, 0, 1, Node(5, None, 0.0, 1.0, None, "extant"), 1.0, rearr,
                 np.random.default_rng(0), 0.0)
    assert rearr == [] and len(genome[0].genes) == 2         # nowhere to move to


def test_all_three_rearrangements_fire_and_are_typed():
    _, r = _run(seed=3, inversion=0.3, transposition=0.3, translocation=0.3,
                transposition_extension=Geometric(mean=3), translocation_extension=Geometric(mean=2))
    kinds = {type(x).__name__ for x in r.rearrangements}
    assert {"Inversion", "Transposition", "Translocation"} <= kinds


def test_default_extension_is_a_single_gene_and_scales_up():
    # default Geometric(mean=1): every inversion spans exactly one gene
    _, small = _run(seed=4, inversion=2.0, transposition=0.0, translocation=0.0,
                    duplication=0.0, transfer=0.0, loss=0.0, origination=0.0)
    assert small.rearrangements and all(x.length == 1 for x in small.rearrangements)
    # dial the extension up: longer blocks appear
    _, big = _run(seed=4, inversion=2.0, transposition=0.0, translocation=0.0,
                  duplication=0.0, transfer=0.0, loss=0.0, origination=0.0,
                  inversion_extension=Geometric(mean=6))
    assert max(x.length for x in big.rearrangements) > 1


def test_inversion_probability_governs_flips():
    _, always = _run(seed=2, transposition=1.0, translocation=1.0, inversion=0.0, duplication=0.3,
                     transposition_extension=Geometric(mean=3), translocation_extension=Geometric(mean=3),
                     inversion_probability=1.0)
    moves = [x for x in always.rearrangements if isinstance(x, (Transposition, Translocation))]
    assert moves and all(x.flipped for x in moves)
    _, never = _run(seed=2, transposition=1.0, translocation=1.0, inversion=0.0, duplication=0.3,
                    inversion_probability=0.0)
    moves0 = [x for x in never.rearrangements if isinstance(x, (Transposition, Translocation))]
    assert moves0 and not any(x.flipped for x in moves0)


def test_strong_invariant_holds_under_segmental_everything():
    # segmental duplication/loss/transfer record correctly, so surviving leaves == profile copies.
    # loss >= duplication keeps genomes bounded (segmental dup would otherwise blow up fast)
    for seed in range(3):
        sp, r = _run(seed=seed, n_extant=8, duplication=0.3, loss=0.4, transfer=0.25, inversion=0.2,
                     transposition=0.2, translocation=0.2, inversion_probability=0.5,
                     duplication_extension=Geometric(mean=3), loss_extension=Geometric(mean=3),
                     transfer_extension=Geometric(mean=2))
        extant = {n.id for n in sp.complete_tree.extant()}
        for fam, tree in r.gene_trees.items():
            assert _extant_leaves(tree.extant) == sum(r.profiles.counts.get((fam, s), 0) for s in extant)


def test_rearrangements_tsv_has_one_table_for_all_kinds(tmp_path):
    _, r = _run(seed=3, inversion=0.3, transposition=0.3, translocation=0.3)
    r.write(tmp_path, outputs=("rearrangements",))
    lines = (tmp_path / "rearrangements.tsv").read_text().splitlines()
    assert lines[0].split("\t") == ["time", "kind", "lineage", "chromosome", "start", "length",
                                    "dest_chromosome", "dest_position", "flipped"]
    kinds = {row.split("\t")[1] for row in lines[1:]}
    assert kinds <= {"inversion", "transposition", "translocation"}
