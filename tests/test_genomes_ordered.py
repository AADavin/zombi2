"""Tests for ordered genomes — genes with a position and a strand, on chromosomes (zombi2.genomes).

Slice 1: the chromosome container + D/T/L/O made position-aware + inversions + chromosome identity.
The gene genealogy is the shared, position-blind event log, so the strong gene-tree/profile invariant
must survive unchanged; on top of it we check the two new things — the inversion operator and the
chromosome genealogy.
"""

import inspect

import pytest

from zombi2.genomes.events import node_from_label
from zombi2.rates import scope
from zombi2.rates.distributions import Fixed, Geometric
from zombi2.rates.modifiers import OnTime
from zombi2.species import simulate_species_tree
from zombi2.tree import Node, Tree
from zombi2.genomes import (
    Chromosome,
    Gene,
    Inversion,
    Transposition,
    Translocation,
    simulate_genomes_ordered,
    simulate_genomes_unordered,
)
from zombi2.genomes.ordered import (
    _duplicate,
    _extent,
    _fission,
    _fusion,
    _invert,
    _lose_at,
    _transpose,
    _translocate,
)


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
    # the layering contract: unordered ⊂ ordered, only the extra args differ. `parallel` is the one
    # documented exception — it is an unordered-only engine, because per-family parallelism needs the
    # families to be independent, and the ordered resolution couples them by position (an inversion or
    # translocation spans several families), so it can never have a per-family engine.
    shared = set(inspect.signature(simulate_genomes_unordered).parameters) - {"tree", "parallel"}
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
    assert rearr == [Inversion(3.0, 7, 0, 1, 2)]           # start 1, length 2
    assert (rearr[0].start, rearr[0].length) == (1, 2)


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
    # every rearrangement names its run the same way: a start position and a length in genes
    _, r = _run(seed=4)
    assert r.rearrangements
    for inv in r.rearrangements:
        assert inv.start >= 0 and inv.length >= 1


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
    r.write(tmp_path, outputs=("events", "profiles", "gene_order",
                               "chromosome_events"))
    for name in ("genome_events.tsv", "profiles.tsv", "gene_order.tsv",
                 "chromosome_events.tsv"):
        assert (tmp_path / name).exists()
    head = (tmp_path / "gene_order.tsv").read_text().splitlines()[0]
    assert head.split("\t") == ["lineage", "chromosome", "position", "strand", "family", "copy"]


def _written_gene_order(path):
    """``gene_order.tsv`` -> ``{node: [(chromosome, position, strand, family, gene), ...]}``."""
    lines = (path / "gene_order.tsv").read_text().splitlines()[1:]
    written = {}
    for row in lines:
        s, *rest = row.split("\t")
        s, rest = node_from_label(s), [int(c) for c in rest]
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
    r = simulate_genomes_ordered(sp, chromosomes=3, initial_families=0, seed=1)   # no families, no events
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


def test_one_table_carries_the_genealogy_its_places_and_the_rearrangements(tmp_path):
    """Three files became one. The genealogy is unchanged, each event carries the arc it acted on —
    once, on its first row, because the arc is the event's and not each copy's — and the
    ancestry-neutral rearrangements are interleaved by time."""
    _, r = _run(seed=3, inversion=0.3, transposition=0.3, translocation=0.3)
    r.write(tmp_path, outputs=("events",))
    lines = (tmp_path / "genome_events.tsv").read_text().splitlines()
    assert lines[0].split("\t") == ["time", "kind", "lineage", "family", "copy", "parent",
                                    "recipient", "donor", "dest_lineage", "chromosome", "position",
                                    "length", "dest_chromosome", "dest_position", "flipped"]
    rows = [ln.split("\t") for ln in lines[1:]]
    assert [float(x[0]) for x in rows] == sorted(float(x[0]) for x in rows)   # in the order it happened

    genealogy = {"origination", "duplication", "loss", "transfer", "speciation"}
    assert len([x for x in rows if x[1] in genealogy]) == len(r.events)
    assert len([x for x in rows if x[1] not in genealogy]) == len(r.rearrangements)
    assert {x[1] for x in rows} - genealogy <= {"inversion", "transposition", "translocation"}
    # one arc per positioned event, not one per copy it touched
    assert len([x for x in rows if x[1] in genealogy and x[10]]) == len(r.event_positions)
    assert not [x for x in rows if x[1] == "speciation" and x[10]]     # a speciation moves nothing
    # a rearrangement ends no gene lineage, so it names no family, copy or parent
    assert all(x[3] == x[4] == x[5] == "" for x in rows if x[1] not in genealogy)


# --- topology: a circular chromosome has no ends, so a run wraps past position 0 ------------------

def _lone_branch(total_time):
    return Tree({0: Node(0, None, 0.0, total_time, None, "extant")}, 0)


def _inversion_coverage(topology, n=8, mean=4.0, total_time=3000.0, seed=1):
    """Inversions only, on one chromosome of ``n`` genes: how often each position ends up inside a
    run. An inversion never creates or destroys a gene, so the chromosome keeps its ``n`` positions
    all run long and the tally is comparable across them. Returns ``(result, coverage per position)``."""
    r = simulate_genomes_ordered(_lone_branch(total_time), inversion=1.0, chromosomes=1,
                                 topology=topology, initial_families=n,
                                 inversion_extension=Geometric(mean=mean), seed=seed)
    cov = [0] * n
    for x in r.rearrangements:
        for k in range(x.length):
            cov[(x.start + k) % n] += 1
    return r, cov


def test_extent_wraps_on_a_circle_and_stops_at_the_end_of_a_line():
    import numpy as np
    circ = Chromosome(0, "circular", [Gene(i, i, 1) for i in range(8)])
    lin = Chromosome(1, "linear", [Gene(i, i, 1) for i in range(8)])
    rng = np.random.default_rng(0)
    big = Geometric(mean=20)
    # from position 6 a linear run reaches the last gene and stops: at most 2 genes
    assert all(_extent(rng, big, lin, 6) <= 2 for _ in range(50))
    # a circular one carries on past position 0
    assert max(_extent(rng, big, circ, 6) for _ in range(50)) > 2


def test_a_run_never_exceeds_the_whole_chromosome():
    # m >= n: a run cannot wrap onto itself, so it is clamped to the whole chromosome
    import numpy as np
    rng = np.random.default_rng(0)
    circ = Chromosome(0, "circular", [Gene(i, i, 1) for i in range(5)])
    lin = Chromosome(1, "linear", [Gene(i, i, 1) for i in range(5)])
    huge = Fixed(1000)
    assert [_extent(rng, huge, circ, s) for s in range(5)] == [5, 5, 5, 5, 5]
    assert [_extent(rng, huge, lin, s) for s in range(5)] == [5, 4, 3, 2, 1]


def test_a_wrapped_inversion_reverses_the_run_across_the_origin():
    ch = Chromosome(0, "circular", [Gene(i, i, 1) for i in range(4)])
    _invert(ch, 3, 3, Node(7, None, 0.0, 1.0, None, "extant"), 2.0, rearr := [])
    # the run is positions 3, 0, 1 — genes 3, 0, 1 — reversed to 1, 0, 3 with strands flipped
    assert [g.id for g in ch.genes] == [1, 0, 3, 2]
    assert [g.strand for g in ch.genes] == [-1, -1, -1, 1]
    # recorded in the frame the chromosome had before the event, so start + length exceeds 4 genes
    assert rearr == [Inversion(2.0, 7, 0, 3, 3)]


def test_a_whole_chromosome_inversion_reverses_the_ring():
    ch = Chromosome(0, "circular", [Gene(i, i, 1) for i in range(4)])
    _invert(ch, 2, 4, Node(7, None, 0.0, 1.0, None, "extant"), 1.0, rearr := [])
    # every gene is in the run: the whole ring reverses — the same molecule read the other way
    assert [g.id for g in ch.genes] == [1, 0, 3, 2]
    assert all(g.strand == -1 for g in ch.genes)
    assert rearr == [Inversion(1.0, 7, 0, 2, 4)]


def test_a_wrapped_duplication_keeps_the_block_together():
    ch = Chromosome(0, "circular", [Gene(i, i, 1) for i in range(4)])
    events, positions, counter = [], [], [10]

    def ng(fam, strand):
        counter[0] += 1
        return Gene(counter[0], fam, strand)
    added = _duplicate(ch, 3, 2, Node(3, None, 0.0, 1.0, None, "extant"), 1.0, events, positions, ng)
    # the run is families 3 then 0, across the origin; its tandem copy lands right behind it
    assert added == 2 and [g.family for g in ch.genes] == [3, 0, 3, 0, 1, 2]
    assert len(events) == 4 and all(e.kind == "duplication" for e in events)
    # the position is recorded in the re-anchored frame, where the run starts at 0
    assert [(p.start, p.length) for p in positions] == [(0, 2)]


def test_a_wrapped_loss_removes_the_genes_on_both_sides_of_the_origin():
    ch = Chromosome(0, "circular", [Gene(i, i, 1) for i in range(5)])
    events, positions = [], []
    removed = _lose_at(ch, 4, 3, Node(3, None, 0.0, 1.0, None, "extant"), 1.0, events, positions)
    assert removed == 3 and [g.id for g in ch.genes] == [2, 3]   # genes 4, 0 and 1 went
    assert sorted(e.copy for e in events) == [0, 1, 4]
    assert [(p.start, p.length) for p in positions] == [(0, 3)]


def test_a_whole_chromosome_loss_empties_it_but_leaves_the_chromosome():
    # a run covering every gene is legal. The chromosome survives as an empty replicon, exactly as a
    # de-novo one starts out; only chromosome_loss takes a chromosome out of the karyotype.
    r = simulate_genomes_ordered(_lone_branch(5.0), loss=2.0, chromosomes=2, initial_families=6,
                                 loss_extension=Fixed(50), seed=1)
    assert len(r.genomes[0]) == 2                                # both chromosomes still there
    assert sum(len(ch.genes) for ch in r.genomes[0]) == 0        # and both empty
    assert sorted(e.copy for e in r.events if e.kind == "loss") == list(range(6))


def test_a_linear_chromosome_still_clamps_at_its_end():
    r = simulate_genomes_ordered(_lone_branch(200.0), inversion=1.0, chromosomes=1,
                                 topology="linear", initial_families=8,
                                 inversion_extension=Geometric(mean=6), seed=1)
    assert r.rearrangements
    assert all(x.start + x.length <= 8 for x in r.rearrangements)


def test_a_circular_chromosome_really_wraps():
    r = simulate_genomes_ordered(_lone_branch(200.0), inversion=1.0, chromosomes=1,
                                 topology="circular", initial_families=8,
                                 inversion_extension=Geometric(mean=6), seed=1)
    assert any(x.start + x.length > 8 for x in r.rearrangements)  # runs cross position 0
    assert all(x.length <= 8 for x in r.rearrangements)           # never more than the whole ring


def test_segmental_events_cover_a_circle_evenly():
    # Translation invariance: a ring has no special position, so every gene must be covered at the
    # same rate. Clamping runs at the end of the gene array broke this — the first gene was covered
    # only when a run started exactly on it, a factor of the mean extension less often than an
    # interior gene. That asymmetry is the bug; it must be gone.
    _, circ = _inversion_coverage("circular")
    mean = sum(circ) / len(circ)
    assert max(abs(c - mean) for c in circ) < 0.10 * mean
    # a linear chromosome keeps the edge behaviour, which is real for a replicon with ends
    _, lin = _inversion_coverage("linear")
    assert lin[0] < 0.5 * lin[-1]


def test_realised_extension_on_a_circle_matches_the_nominal_one():
    # with no end to truncate them, runs on a circle realise E[min(M, n)] — everything the extension
    # distribution asks for, short only of what the chromosome cannot hold
    q = 1 - 1 / 4.0                                   # M ~ Geometric(mean=4); E[min(M, 8)] = sum q^k
    expected = sum(q ** k for k in range(8))
    r, _ = _inversion_coverage("circular")
    realised = sum(x.length for x in r.rearrangements) / len(r.rearrangements)
    assert abs(realised - expected) < 0.05 * expected
    # clamping at the array end instead loses a good fraction of it
    r2, _ = _inversion_coverage("linear")
    short = sum(x.length for x in r2.rearrangements) / len(r2.rearrangements)
    assert short < 0.85 * expected


def test_the_strong_invariant_survives_wrapped_runs():
    # runs longer than the chromosome, on circles, so most events cross the origin: the gene
    # genealogy must still account for every surviving copy
    exts = {f"{e}_extension": Geometric(mean=6) for e in
            ("duplication", "loss", "transfer", "inversion", "transposition", "translocation")}
    for seed in range(3):
        sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=seed)
        r = simulate_genomes_ordered(sp, duplication=0.25, loss=0.35, transfer=0.2, inversion=0.3,
                                     transposition=0.2, translocation=0.2, chromosomes=3,
                                     initial_families=9, inversion_probability=0.5, seed=seed,
                                     **exts)
        extant = {n.id for n in sp.complete_tree.extant()}
        for fam, tree in r.gene_trees.items():
            assert _extant_leaves(tree.extant) == sum(r.profiles.counts.get((fam, s), 0)
                                                      for s in extant)


def test_wrapped_runs_stay_deterministic_given_a_seed():
    exts = {f"{e}_extension": Geometric(mean=5) for e in
            ("duplication", "loss", "transfer", "inversion", "transposition", "translocation")}
    kw = dict(duplication=0.3, loss=0.35, transfer=0.2, inversion=0.4, transposition=0.3,
              translocation=0.3, chromosomes=3, initial_families=9, topology="circular",
              inversion_probability=0.5, seed=17, **exts)
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=17)
    a = simulate_genomes_ordered(sp, **kw)
    b = simulate_genomes_ordered(sp, **kw)
    assert a.rearrangements and a.rearrangements == b.rearrangements
    assert all(a.gene_order(x) == b.gene_order(x) for x in a.genomes)
    assert a.events == b.events


# --- scope: a rearrangement is counted per gene, not per chromosome -------------------------------

def _inversions(result):
    return sum(1 for r in result.rearrangements if isinstance(r, Inversion))


def test_a_rearrangement_starts_at_a_gene_not_at_a_chromosome():
    # Inversion and transposition act on a run of genes, so they are scoped per copy: the drawn
    # gene IS the run's start. Picking a chromosome first and a position inside it would make a gene
    # on a small replicon far likelier to be a breakpoint than a gene on a large one.
    from zombi2.genomes import ordered

    calls = {"gene": 0, "chromosome": 0}
    real_gene, real_chrom = ordered._pick_gene, ordered._pick_chromosome

    def spy_gene(*a, **k):
        calls["gene"] += 1
        return real_gene(*a, **k)

    def spy_chrom(*a, **k):
        calls["chromosome"] += 1
        return real_chrom(*a, **k)

    ordered._pick_gene, ordered._pick_chromosome = spy_gene, spy_chrom
    try:
        sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=6, seed=2)
        r = simulate_genomes_ordered(sp, inversion=0.4, transposition=0.3, chromosomes=3,
                                     initial_families=12, seed=2)
    finally:
        ordered._pick_gene, ordered._pick_chromosome = real_gene, real_chrom
    assert r.rearrangements, "the run produced no rearrangements to check"
    assert calls["gene"] > 0
    assert calls["chromosome"] == 0, "a rearrangement must not draw a chromosome first"


def test_rearrangement_count_ignores_how_the_genes_are_split_into_chromosomes():
    # The same genes carved into 1, 2 or 4 chromosomes is the same amount of DNA, so it must give the
    # same number of inversions. Under per-chromosome scope this tripled from left to right, which is
    # also why a fission used to double a genome's inversion rate without creating a single gene.
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=6, seed=1)
    counts = [_inversions(simulate_genomes_ordered(sp, inversion=0.05, initial_families=40,
                                                   chromosomes=c, seed=3))
              for c in (1, 2, 4)]
    assert counts[0] > 0
    assert len(set(counts)) == 1, f"chromosome number changed the inversion count: {counts}"


def test_rearrangement_count_scales_with_gene_count():
    # Twice the DNA, twice the chances to start a run. Averaged over seeds, doubling the genome
    # doubles the inversions.
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=6, seed=1)

    def mean_at(fams):
        runs = [_inversions(simulate_genomes_ordered(sp, inversion=0.05, initial_families=fams,
                                                     chromosomes=1, seed=s))
                for s in range(40)]
        return sum(runs) / len(runs)

    small, large = mean_at(20), mean_at(40)
    assert small > 0
    assert 1.7 < large / small < 2.3, f"expected ~2x, got {large / small:.2f} ({small} -> {large})"


def test_the_initial_genome_is_the_layout_the_run_started_with(tmp_path):
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=5, seed=1)
    g = simulate_genomes_ordered(sp, duplication=0.5, loss=0.5, inversion=2.0,
                                 initial_families=8, seed=1)
    assert [len(c.genes) for c in g.initial_genome] == [8]
    assert [gene.strand for c in g.initial_genome for gene in c.genes] == [1] * 8
    root = sp.complete_tree.root
    assert g.gene_order(root) != [(c.id, p, gn.strand, gn.family, gn.id)
                                  for c in g.initial_genome for p, gn in enumerate(c.genes)], \
        "the stem was quiet — pick another seed"
    g.write(tmp_path)
    rows = (tmp_path / "initial_genome.tsv").read_text().splitlines()
    assert rows[0] == "chromosome\tposition\tstrand\tfamily\tcopy" and len(rows) == 9
