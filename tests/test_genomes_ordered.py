"""Tests for ordered genomes — genes with a position and a strand, on chromosomes (zombi2.genomes).

Slice 1: the chromosome container + D/T/L/O made position-aware + inversions + chromosome identity.
The gene genealogy is the shared, position-blind event log, so the strong gene-tree/profile invariant
must survive unchanged; on top of it we check the two new things — the inversion operator and the
chromosome genealogy.
"""

import hashlib
import inspect

import numpy as np
import pytest

from zombi2.rates import modifiers as mod

from zombi2.genomes.events import gene_from_label, node_from_label
from zombi2.rates import scope
from zombi2.rates.distributions import Fixed, Geometric
from zombi2.rates.modifiers import ByFamily, ByLineage, FromParent, OnTime, OnTotalDiversity
from zombi2.species import simulate_species_tree
from zombi2.tree import Node, Tree
from zombi2.genomes import (
    Between,
    Chromosome,
    Clades,
    Distance,
    Gene,
    Inversion,
    Transposition,
    Translocation,
    simulate_genomes_ordered,
    simulate_genomes_family,
)
from zombi2.genomes.ordered import (
    _do_transfer,
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
    # the layering contract: family ⊂ ordered, only the extra args differ. `parallel` / `stream_to` /
    # `outputs` are the documented exceptions — the per-family engine (and its streaming form) is
    # family-only, because per-family parallelism needs the families to be independent, and the
    # ordered resolution couples them by position (an inversion or translocation spans several families).
    shared = (set(inspect.signature(simulate_genomes_family).parameters)
              - {"tree", "parallel", "stream_to", "outputs"})
    ordered = set(inspect.signature(simulate_genomes_ordered).parameters) - {"tree"}
    assert shared <= ordered                                 # family ⊂ ordered: nothing dropped
    assert ordered - shared == {                            # ordered's own additions:
        "chromosomes", "topology", "inversion", "transposition", "translocation",
        "fission", "fusion", "chromosome_origination", "chromosome_loss", "inversion_probability",
        "duplication_extent", "loss_extent", "transfer_extent",
        "inversion_extent", "transposition_extent", "translocation_extent"}


# --- the shared gene genealogy still holds -------------------------------------------------------

def test_extant_gene_tree_leaves_equal_the_extant_copy_total():
    # the strongest invariant, inherited from the family core: surviving gene-tree leaves == copies
    sp, r = _run(seed=5, death=0.5)
    extant_sp = {n.id for n in sp.complete_tree.extant_leaves()}
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

def test_initial_chromosomes_are_roots_of_their_own_kind():
    # a replicon the run STARTS with is not something it did, so it is `initial`, not `origination` —
    # the same word the gene log and the block log use at t=0. Both are roots (no parent).
    _, r = _run(seed=8, chromosomes=4)
    roots = [e for e in r.chromosome_events if e.kind == "initial"]
    assert len(roots) == 4
    assert all(e.parents == () and len(e.children) == 1 for e in roots)
    assert not [e for e in r.chromosome_events if e.kind == "origination"]  # none arose de novo


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
    xfer_rows = [e for e in r.edges if e.kind == "transfer"]
    assert xfer_rows                                           # transfers really fired
    # the recipient row carries a recipient different from the donor lineage
    recip = [e for e in xfer_rows if e.recipient is not None]
    assert recip and all(e.lineage == e.recipient for e in recip)


def test_no_transfer_events_when_transfer_is_zero():
    _, r = _run(seed=8, transfer=0.0)
    assert not any(e.kind == "transfer" for e in r.edges)


def test_replacement_run_stays_consistent():
    # replacement overwrites a homologous copy; the strong invariant must still hold
    sp, r = _run(seed=2, replacement=True)
    extant_sp = {n.id for n in sp.complete_tree.extant_leaves()}
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


#: The digest of one seeded run exercising **every** event class, captured BEFORE ``DrivenBy`` was
#: wired into this engine. Every driven-path addition is behind ``if any_driven`` / a ``w`` argument
#: and draws nothing from the rng, so an undriven run must hash the same: the draw order of the plain
#: path is what a hundred seeded tests, the gallery and every analysis depend on.
_UNDRIVEN_ORDERED_DIGEST = "2c2b782b7bd55a2197dbb153aabfd5e34ccc39a7f557bc03490bd3184f509c06"


def _ordered_digest(r) -> str:
    """Everything the run produced, hashed: the gene genealogy and where each event happened, the
    rearrangements, the chromosome network, every node's layout, and the initial genome."""
    import hashlib
    key = repr([
        [(round(e.time, 12), e.kind, e.lineage, e.family, e.copy, e.parent, e.recipient)
         for e in r.edges],
        [(round(p.time, 12), p.kind, p.lineage, p.chromosome, p.start, p.length, p.family,
          p.donor, p.recipient, p.dest_position) for p in r.event_positions],
        [(type(x).__name__, round(x.time, 12), x.lineage, tuple(sorted(vars(x).items())))
         for x in r.rearrangements],
        [(round(c.time, 12), c.kind, c.lineage, c.parents, c.children) for c in r.chromosome_events],
        {k: r.gene_order(k) for k in sorted(r.genomes)},
        [(c.id, c.topology, [(g.id, g.family, g.strand) for g in c.genes]) for c in r.initial_genome],
    ])
    return hashlib.sha256(key.encode()).hexdigest()


def test_undriven_ordered_run_is_byte_identical():
    tree = simulate_species_tree(birth=1.2, death=0.2, total_time=2.0, seed=17).complete_tree
    r = simulate_genomes_ordered(
        tree, duplication=0.2, transfer=0.3, loss=0.15, origination=0.25,
        inversion=0.4, transposition=0.2, translocation=0.2,
        chromosomes=3, fission=0.1, fusion=0.1,
        chromosome_origination=0.05, chromosome_loss=0.05,
        duplication_extent=3, loss_extent=2, transfer_extent=2, inversion_extent=4,
        transposition_extent=2, translocation_extent=2, inversion_probability=0.3,
        initial_families=30, seed=23)
    assert _ordered_digest(r) == _UNDRIVEN_ORDERED_DIGEST, (
        "an undriven ordered run changed: the rng draw order of the plain path must not move")


def test_ontime_skyline_modifier_is_accepted():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=1)
    r = simulate_genomes_ordered(sp, duplication=0.3 * OnTime({0: 1.0, 1.0: 0.2}),
                                 inversion=0.2 * OnTime({0: 0.5, 1.0: 2.0}),
                                 chromosomes=2, initial_families=6, seed=1)
    assert r.genomes                                           # ran without complaint


@pytest.mark.parametrize("modifier", [ByLineage(spread=0.5), FromParent(spread=0.5),
                                      OnTotalDiversity(cap=100)])
def test_unsupported_modifier_is_rejected(modifier):
    """The gate had no test at all. A modifier the engine cannot read must raise — one that returns
    1.0 because nothing looks at it is a run quietly not the model asked for (SPEC §5) — and the
    message must name what this engine *does* take, so the reader knows where to go next."""
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=1)
    with pytest.raises(ValueError, match="ordered genome engine does not support"):
        simulate_genomes_ordered(sp, inversion=0.3 * modifier, initial_families=6, seed=1)


@pytest.mark.parametrize("modifier", [ByLineage(spread=0.5), FromParent(spread=0.5)])
def test_unsupported_modifier_on_an_extent_is_rejected(modifier):
    """An extent takes the same modifiers a rate does (SPEC §6), so the same refusal applies — and
    names the two it takes rather than pointing at another resolution."""
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=1)
    with pytest.raises(ValueError, match="does not support on an extent"):
        simulate_genomes_ordered(sp, inversion=0.3, inversion_extent=3 * modifier,
                                 initial_families=6, seed=1)


def test_the_extent_declaration_is_the_rate_declaration_minus_byfamily():
    """The one difference between the two lists is a modelling fact, not an accident: ``ByFamily``
    attaches to the contents, and an extent is drawn before the run's genes are known."""
    from zombi2.genomes.ordered import WIRED_EXTENT_MODIFIERS, WIRED_MODIFIERS
    assert set(WIRED_MODIFIERS) - set(WIRED_EXTENT_MODIFIERS) == {ByFamily}


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
    head = (tmp_path / "gene_order.tsv").read_text(encoding="utf-8").splitlines()[0]
    assert head.split("\t") == ["lineage", "chromosome", "topology", "position", "strand",
                                "family", "copy"]


def _written_gene_order(path):
    """``gene_order.tsv`` -> ``{node: [(chromosome, position, strand, family, gene), ...]}``.

    The ``topology`` column is dropped here on purpose: this helper exists to compare the file
    against ``result.gene_order()``, which is the arrangement and carries no topology. What the
    column says is checked on its own, in ``test_gene_order_records_each_chromosome_topology``."""
    lines = (path / "gene_order.tsv").read_text(encoding="utf-8").splitlines()[1:]
    written = {}
    for row in lines:
        s, chrom, _topology, *rest = row.split("\t")
        *ints, copy = rest                       # position, strand, family are ints; copy is g<n>
        written.setdefault(node_from_label(s), []).append(
            (int(chrom), *(int(c) for c in ints), gene_from_label(copy)))
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
    assert {"initial", "origination", "speciation", "fission", "fusion", "loss"} <= seen  # all occur
    shape = {"initial": (0, 1), "origination": (0, 1), "speciation": (1, 2), "fission": (1, 2),
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
        extant = {n.id for n in sp.complete_tree.extant_leaves()}
        for fam, tree in r.gene_trees.items():
            assert _extant_leaves(tree.extant) == sum(r.profiles.counts.get((fam, s), 0) for s in extant)


def test_a_genome_never_loses_its_last_chromosome():
    _, r = _tier(seed=4, chromosome_loss=1.0, chromosome_origination=0.0)  # push loss hard
    assert all(len(chroms) >= 1 for chroms in r.genomes.values())


def test_a_genome_never_loses_its_last_genes_to_a_chromosome_loss():
    """"Not the genome's last chromosome" is not enough to keep a genome alive.

    A lineage can be carrying **empty** replicons — `chromosome_origination` mints one empty, and a
    translocation can empty one — so a genome of one gene-bearing chromosome beside an empty plasmid
    used to lose everything the moment `chromosome_loss` picked the one with the genes on it. The
    parameters below are the ones that reproduced it: 3 of 6 extant genomes came out with no genes."""
    for seed in range(40):
        sp = simulate_species_tree(birth=1.0, n_extant=6, seed=seed)
        r = simulate_genomes_ordered(sp.complete_tree, loss=3.0, translocation=3.0,
                                     chromosome_loss=1.0, chromosome_origination=0.5,
                                     initial_families=4, chromosomes=2, seed=seed)
        for node in sp.complete_tree.extant_leaves():
            assert any(c.genes for c in r.genomes[node.id]), (
                f"seed {seed}, lineage {node.id}: every chromosome came out empty")


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
    de_novo = [e for e in r.chromosome_events if e.kind == "origination"]
    assert de_novo                                    # some replicons appeared past the initial one
    # the initial replicon carried the 3 genes; the de-novo ones are empty
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


def test_a_circular_chromosome_never_fuses_with_a_linear_one():
    # A ring and a molecule with two ends cannot become one molecule. This used to draw the partner
    # uniformly over the whole karyotype and hand the child `a.topology`, so the two silently fused
    # into one chromosome whose shape was whichever end the chromosome pick landed on first.
    import numpy as np
    genome = [Chromosome(10, "circular", [Gene(0, 0, 1), Gene(1, 1, 1)]),
              Chromosome(11, "linear", [Gene(2, 2, 1)])]
    node = Node(5, None, 0.0, 1.0, None, "extant")
    ce = []
    for ci in (0, 1):                                  # neither one has a legal partner
        assert _fusion(genome, ci, node, 3.0, ce, _minter(20), np.random.default_rng(0)) == (0, 0)
    assert [c.id for c in genome] == [10, 11]          # the karyotype is exactly as it was
    assert [c.topology for c in genome] == ["circular", "linear"]
    assert ce == []                                    # a declined event logs nothing


def test_fusion_picks_a_partner_of_its_own_topology():
    # With one legal partner among two candidates the choice is forced, so this holds at every seed:
    # the circular chromosome fuses with the other circular one and the linear one is left alone.
    import numpy as np
    for seed in range(20):
        genome = [Chromosome(0, "circular", [Gene(0, 0, 1)]),
                  Chromosome(1, "linear", [Gene(1, 1, 1)]),
                  Chromosome(2, "circular", [Gene(2, 2, 1)])]
        node = Node(5, None, 0.0, 1.0, None, "extant")
        ce = []
        assert _fusion(genome, 0, node, 3.0, ce, _minter(90),
                       np.random.default_rng(seed)) == (-1, 0)
        assert ce[0].parents == (0, 2)                          # never (0, 1)
        assert [c.topology for c in genome] == ["linear", "circular"]   # the linear one survives
        assert [g.id for g in genome[1].genes] == [0, 2]


def test_a_single_topology_genome_draws_its_partner_exactly_as_before():
    # The rule must cost nothing where it cannot bite. In a genome of one topology `partners` is
    # every other chromosome in index order, so it is the same single draw from a pool of n-1 that
    # the old `rng.integers(n - 1)` + skip-past-ci arithmetic made, and it maps to the same
    # chromosome. Pinning it here is what makes "a seeded circular run is unchanged" a test rather
    # than a claim.
    import numpy as np
    node = Node(5, None, 0.0, 1.0, None, "extant")
    for seed in range(10):
        for ci in range(6):
            genome = [Chromosome(k, "circular", [Gene(k, k, 1)]) for k in range(6)]
            ce = []
            _fusion(genome, ci, node, 3.0, ce, _minter(90), np.random.default_rng(seed))
            i = int(np.random.default_rng(seed).integers(5))    # the one draw, from a pool of n-1
            expected = i if i < ci else i + 1                   # ... mapped as the old code mapped it
            assert ce[0].parents == (ci, expected)


def test_a_mixed_topology_karyotype_keeps_one_chromosome_of_each_shape():
    # The end-to-end regression. Under a hard fusion rate a karyotype of two rings and two linear
    # molecules used to collapse to a single chromosome carrying every gene; now it can only ever
    # collapse to one of each shape, because the last ring and the last linear molecule have no
    # legal partner left.
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=4, seed=1)
    r = simulate_genomes_ordered(sp, fusion=1.0, fission=0.5, chromosomes=4,
                                 topology=["circular", "linear", "circular", "linear"],
                                 initial_families=12, seed=5)
    for n in sp.complete_tree.extant_leaves():
        assert sorted(c.topology for c in r.genomes[n.id]) == ["circular", "linear"]
    # ... and no recorded fusion edge ever joined two shapes. Ancestral chromosomes are not in the
    # result, so track each id's topology down the network from the roots the run laid down.
    topology = {}
    initial = iter(["circular", "linear", "circular", "linear"])
    for e in r.chromosome_events:
        if e.kind == "initial":
            topology[e.children[0]] = next(initial)
        elif e.kind == "origination":
            topology[e.children[0]] = "circular"            # a de-novo replicon is a plasmid
        else:
            assert len({topology[p] for p in e.parents}) == 1
            for ch in e.children:
                topology[ch] = topology[e.parents[0]]
    assert any(e.kind == "fusion" for e in r.chromosome_events)      # the rule was actually exercised


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
                transposition_extent=Geometric(mean=3), translocation_extent=Geometric(mean=2))
    kinds = {type(x).__name__ for x in r.rearrangements}
    assert {"Inversion", "Transposition", "Translocation"} <= kinds


def test_default_extent_is_a_single_gene_and_scales_up():
    # default Geometric(mean=1): every inversion spans exactly one gene
    _, small = _run(seed=4, inversion=2.0, transposition=0.0, translocation=0.0,
                    duplication=0.0, transfer=0.0, loss=0.0, origination=0.0)
    assert small.rearrangements and all(x.length == 1 for x in small.rearrangements)
    # dial the extension up: longer blocks appear
    _, big = _run(seed=4, inversion=2.0, transposition=0.0, translocation=0.0,
                  duplication=0.0, transfer=0.0, loss=0.0, origination=0.0,
                  inversion_extent=Geometric(mean=6))
    assert max(x.length for x in big.rearrangements) > 1


def test_inversion_probability_governs_flips():
    _, always = _run(seed=2, transposition=1.0, translocation=1.0, inversion=0.0, duplication=0.3,
                     transposition_extent=Geometric(mean=3), translocation_extent=Geometric(mean=3),
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
                     duplication_extent=Geometric(mean=3), loss_extent=Geometric(mean=3),
                     transfer_extent=Geometric(mean=2))
        extant = {n.id for n in sp.complete_tree.extant_leaves()}
        for fam, tree in r.gene_trees.items():
            assert _extant_leaves(tree.extant) == sum(r.profiles.counts.get((fam, s), 0) for s in extant)


def test_the_genealogy_and_the_rearrangements_are_two_tables(tmp_path):
    """One row per event, and one file per kind of statement. ``genome_events.tsv`` is the shared
    genealogy — ``parents`` and ``children`` naming copies where they sit — with the arc each event
    acted on beside it. ``rearrangement_events.tsv`` is what moved without beginning or ending a
    lineage, and so has nothing to say about parents and children at all."""
    _, r = _run(seed=3, inversion=0.3, transposition=0.3, translocation=0.3)
    r.write(tmp_path, outputs=("events",))
    lines = (tmp_path / "genome_events.tsv").read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    assert header == ["time", "kind", "family", "parents", "children",
                      "chromosome", "start", "length", "dest_chromosome", "dest_position"]
    # keyed by name, not position, so adding a column cannot silently shift what this asserts on
    at = {c: i for i, c in enumerate(header)}
    rows = [ln.split("\t") for ln in lines[1:]]
    col = lambda x, c: x[at[c]]                                              # noqa: E731
    assert ([float(col(x, "time")) for x in rows]
            == sorted(float(col(x, "time")) for x in rows))   # in the order it happened

    # a genome event is never a rearrangement now: those have a file of their own
    assert {col(x, "kind") for x in rows} <= {"origination", "duplication", "loss", "speciation",
                                              "transfer_additive", "transfer_replacing"}
    # one row per EVENT, not per gene-tree edge: the kinds that end a gene and start two write one
    assert len(rows) == len({(e.kind, e.event) for e in r.edges})
    # every event that moved genes carries where it happened; a speciation moves nothing and never has
    assert all(col(x, "chromosome") for x in rows if col(x, "kind") != "speciation")
    assert not [x for x in rows if col(x, "kind") == "speciation" and col(x, "chromosome")]
    # every participant carries its own branch, so no row needs a lineage column
    assert all(all(tok.startswith("n") or tok.startswith("e")
                   for cell in (col(x, "parents"), col(x, "children"))
                   for tok in cell.split(";") if tok) for x in rows)

    rear = (tmp_path / "rearrangement_events.tsv").read_text(encoding="utf-8").splitlines()
    assert rear[0].split("\t") == ["time", "kind", "lineage", "chromosome", "start", "length",
                                   "dest_chromosome", "dest_position", "flipped"]
    assert len(rear) - 1 == len(r.rearrangements)
    assert {ln.split("\t")[1] for ln in rear[1:]} == {"inversion", "transposition", "translocation"}


# --- topology: a circular chromosome has no ends, so a run wraps past position 0 ------------------

def _lone_branch(total_time):
    return Tree({0: Node(0, None, 0.0, total_time, None, "extant")}, 0)


def _inversion_coverage(topology, n=8, mean=4.0, total_time=3000.0, seed=1):
    """Inversions only, on one chromosome of ``n`` genes: how often each position ends up inside a
    run. An inversion never creates or destroys a gene, so the chromosome keeps its ``n`` positions
    all run long and the tally is comparable across them. Returns ``(result, coverage per position)``."""
    r = simulate_genomes_ordered(_lone_branch(total_time), inversion=1.0, chromosomes=1,
                                 topology=topology, initial_families=n,
                                 inversion_extent=Geometric(mean=mean), seed=seed)
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


def test_a_loss_never_takes_a_chromosome_below_its_last_gene():
    # A run covering every gene left does not fire — the same floor `Chromosome.delete()` enforces at
    # the nucleotide resolution, so the two agree on what a chromosome is. Emptying the karyotype is
    # chromosome_loss's job, at the chromosome tier.
    ch = Chromosome(0, "circular", [Gene(i, i, 1) for i in range(3)])
    events, positions = [], []
    node = Node(3, None, 0.0, 1.0, None, "extant")
    assert _lose_at(ch, 0, 3, node, 1.0, events, positions) == 0   # the whole chromosome: refused
    assert [g.id for g in ch.genes] == [0, 1, 2]                   # and the order is untouched
    assert events == [] and positions == []                        # a declined event logs nothing
    assert _lose_at(ch, 0, 2, node, 1.0, events, positions) == 2   # one short: fires
    assert [g.id for g in ch.genes] == [2]


def test_a_crushing_loss_rate_leaves_every_chromosome_standing():
    r = simulate_genomes_ordered(_lone_branch(5.0), loss=2.0, chromosomes=2, initial_families=6,
                                 loss_extent=Fixed(50), seed=1)
    assert len(r.genomes[0]) == 2                                # both chromosomes still there
    assert all(len(ch.genes) >= 1 for ch in r.genomes[0])        # and neither was emptied


def test_a_linear_chromosome_still_clamps_at_its_end():
    r = simulate_genomes_ordered(_lone_branch(200.0), inversion=1.0, chromosomes=1,
                                 topology="linear", initial_families=8,
                                 inversion_extent=Geometric(mean=6), seed=1)
    assert r.rearrangements
    assert all(x.start + x.length <= 8 for x in r.rearrangements)


def test_a_circular_chromosome_really_wraps():
    r = simulate_genomes_ordered(_lone_branch(200.0), inversion=1.0, chromosomes=1,
                                 topology="circular", initial_families=8,
                                 inversion_extent=Geometric(mean=6), seed=1)
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


def test_realised_extent_on_a_circle_matches_the_nominal_one():
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
    exts = {f"{e}_extent": Geometric(mean=6) for e in
            ("duplication", "loss", "transfer", "inversion", "transposition", "translocation")}
    for seed in range(3):
        sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=seed)
        r = simulate_genomes_ordered(sp, duplication=0.25, loss=0.35, transfer=0.2, inversion=0.3,
                                     transposition=0.2, translocation=0.2, chromosomes=3,
                                     initial_families=9, inversion_probability=0.5, seed=seed,
                                     **exts)
        extant = {n.id for n in sp.complete_tree.extant_leaves()}
        for fam, tree in r.gene_trees.items():
            assert _extant_leaves(tree.extant) == sum(r.profiles.counts.get((fam, s), 0)
                                                      for s in extant)


def test_wrapped_runs_stay_deterministic_given_a_seed():
    exts = {f"{e}_extent": Geometric(mean=5) for e in
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
    rows = (tmp_path / "initial_genome.tsv").read_text(encoding="utf-8").splitlines()
    assert rows[0] == "chromosome\ttopology\tposition\tstrand\tfamily\tcopy" and len(rows) == 9


# --- transfer_to: steering who receives -----------------------------------------------------------
# The recipient rule is the family core's choice slot (SPEC §5), and the whole of it works here: the
# numbers are weights normalised over the contemporaneous candidates, so they change who receives and
# never how many transfers happen. What is ordered about an ordered transfer is the *block* that
# moves, which is the extent — a separate axis (SPEC §6).

_STEERED = dict(transfer=1.0, initial_families=12, max_family_size=8, seed=11)


def _arrivals(result):
    """Every fired transfer's arrival row — the recipient's half of the horizontal edge."""
    return [e for e in result.edges if e.kind == "transfer" and e.recipient is not None]


def _sha(obj):
    return hashlib.sha256(repr(obj).encode()).hexdigest()


def _event_digest(result):
    return _sha([(round(e.time, 12), e.kind, e.lineage, e.family, e.copy, e.parent, e.recipient)
                 for e in result.edges])


def _layout_digest(result):
    """Every node's layout. The genealogy alone is not enough here: rotating a circular chromosome
    moves no gene between lineages and would slip past an events-only digest, and rotation is exactly
    what the recipient-pick hoist below changes the timing of."""
    return _sha({i: result.gene_order(i) for i in sorted(result.genomes)})


# Captured from the engine BEFORE the choice slot was wired at this resolution. Wiring it added a
# branch to `_do_transfer` and moved the recipient pick ABOVE the donor's anchoring, so these are the
# guard that a run which steers nothing did not move: the pick consumes the rng and anchoring does
# not, which is what makes the two free to swap. Both rules are pinned because they take different
# branches of the pick.
_UNDRIVEN_ORDERED_DIGESTS = {
    "uniform": ("67da8d3c71f6f4a58235efb854e1cbddd5bf6660514a59cc0107375ac55ed087",
                "c865dfe27d521e67788e3d2b6f367435282fa754db89809ee864808adfc0643b"),
    "distance": ("4c1a7cd9699709e1bea8144108a44732f9e499328385058aed70d4536d9b41b2",
                 "68c7f1c4b720b34b6f2cea9d84623595aa083b58638f5b7e048e4e53834de035"),
}


@pytest.mark.parametrize("rule", ["uniform", "distance"])
def test_an_undriven_ordered_transfer_is_unchanged(rule):
    tree = simulate_species_tree(birth=1.2, death=0.2, total_time=2.5, seed=17).complete_tree
    r = simulate_genomes_ordered(tree, duplication=0.2, transfer=0.4, loss=0.15, origination=0.3,
                                 inversion=0.2, transposition=0.1, translocation=0.1, chromosomes=2,
                                 transfer_to=rule, initial_families=8, seed=23)
    assert (_event_digest(r), _layout_digest(r)) == _UNDRIVEN_ORDERED_DIGESTS[rule], (
        f"an undriven {rule} transfer changed: the rng draw order of the undriven path must not move")


def test_a_wrapping_thinned_transfer_run_is_unchanged():
    """The same guard on the configuration that exercises the two things the hoist sits between: a
    transfer extent of three genes on a one-chromosome circular genome, so runs wrap position 0 and
    are anchored, and a cap of three copies, so some of them are thinned after the anchoring."""
    tree = simulate_species_tree(birth=1.3, death=0.1, total_time=2.0, seed=4).complete_tree
    r = simulate_genomes_ordered(tree, duplication=0.4, transfer=0.6, loss=0.1, origination=0.2,
                                 transfer_extent=3.0, chromosomes=1, transfer_to="distance",
                                 max_family_size=3, initial_families=6, seed=7)
    assert _event_digest(r) == "fd7d4ab02d283f81ce5072a37eb97d8d2a0afc671ab44771aa5c93126ac46dc4"
    assert _layout_digest(r) == "c0d19a80a8788ad75d0490628095600103ebac5a28cd16c2b6db4c121914a65d"


def test_clades_steer_an_ordered_transfer_between_two_clades():
    """The case a per-recipient weight cannot express, at this resolution: transfer runs strictly
    BETWEEN A and B — never A→A, B→B, or anything touching the rest of the tree. Same two clades and
    the same kernel as the family test, so what is being shown is that the resolution does not matter
    to the rule."""
    from test_genomes_family import _clade_pairs, _two_clades

    sp = simulate_species_tree(birth=1.0, death=0.5, n_extant=20, seed=7)
    a, b, _, _, lab = _two_clades(sp)
    r = simulate_genomes_ordered(
        sp, transfer_to=Clades({"A": a, "B": b},
                               Between({("A", "B"): 1.0, ("B", "A"): 1.0}, default=0.0)),
        **_STEERED)
    pairs = _clade_pairs(r.edges, lab)
    assert pairs
    assert all(p in {("A", "B"), ("B", "A")} for p in pairs)


def test_clades_steer_an_ordered_transfer_in_one_direction_only():
    """A donates, B receives, nothing else — the kernel is a directed one."""
    from test_genomes_family import _clade_pairs, _two_clades

    sp = simulate_species_tree(birth=1.0, death=0.5, n_extant=20, seed=7)
    a, b, _, _, lab = _two_clades(sp)
    r = simulate_genomes_ordered(
        sp, transfer_to=Clades({"A": a, "B": b}, Between({("A", "B"): 1.0}, default=0.0)),
        **_STEERED)
    pairs = _clade_pairs(r.edges, lab)
    assert pairs
    assert all(p == ("A", "B") for p in pairs)


def test_steering_composes_with_replacement():
    """Who receives (the choice slot) and which homolog inside the recipient is overwritten
    (``replacement``) are two independent questions, so steering must leave replacement working."""
    from test_genomes_family import _clade_pairs, _two_clades

    sp = simulate_species_tree(birth=1.0, death=0.5, n_extant=20, seed=7)
    a, b, _, _, lab = _two_clades(sp)
    r = simulate_genomes_ordered(
        sp, replacement=True,
        transfer_to=Clades({"A": a, "B": b},
                           Between({("A", "B"): 1.0, ("B", "A"): 1.0}, default=0.0)),
        **_STEERED)
    pairs = _clade_pairs(r.edges, lab)
    assert pairs
    assert all(p in {("A", "B"), ("B", "A")} for p in pairs)
    assert any(e.replaced is not None for e in _arrivals(r)), "no homolog was displaced — retune"


def test_a_kernel_that_lets_nobody_receive_fires_no_transfer_at_all():
    """Every candidate at weight 0 means the transfer cannot happen, so the event is dropped whole:
    no donor continuation, no arrival, no gene minted. The same run under 'uniform' transfers freely,
    so what is being shown is the weighting and not a dead setup."""
    sp = simulate_species_tree(birth=1.0, death=0.5, n_extant=20, seed=7)
    tips = [i for i, n in sp.complete_tree.nodes.items() if n.children is None]
    blocked = simulate_genomes_ordered(
        sp, transfer_to=Clades({"A": tips[:1]}, Between({("A", "rest"): 0.0}, default=0.0)),
        **_STEERED)
    free = simulate_genomes_ordered(sp, transfer_to="uniform", **_STEERED)
    assert not [e for e in blocked.edges if e.kind == "transfer"]
    assert [e for e in free.edges if e.kind == "transfer"]


def test_a_dropped_transfer_leaves_the_donor_chromosome_untouched():
    """The claim the Poisson-thinning argument rests on, checked directly: a transfer that does not
    fire changes *nothing*.

    It used to be false in a way no run-level assertion would catch. ``_do_transfer`` anchored the
    donor's chromosome — rotating its gene list in place, so a run wrapping position 0 becomes a plain
    slice — *before* the recipient was chosen. A rotation moves no gene between lineages and is
    biologically nothing on a ring, but it renumbers every position the run writes out, so a transfer
    dropped afterwards left the donor visibly changed by an event that did not happen. The recipient
    pick now comes first. The run below wraps (positions 4..7 of a six-gene ring), so the old order
    would have rotated it."""
    gen = [[Chromosome(0, "circular", [Gene(i, i, 1) for i in range(6)])],
           [Chromosome(1, "circular", [Gene(10 + i, i, 1) for i in range(6)])]]
    before = [list(c.genes) for c in gen[0]]
    # both lineages painted "A", and the only pair the kernel weighs is (A, B): nobody can receive
    blocked = Clades({"A": 0, "B": 1}, Between({("A", "B"): 1.0}, default=0.0))
    events, positions = [], []
    delta = _do_transfer(np.random.default_rng(0), None, [0, 1], gen, 0, 0, 4, 4, 1.0,
                         events, positions, None, blocked, False, False, 1.0, None,
                         None, {0: "A", 1: "A"})
    assert delta == 0
    assert not events and not positions
    assert [list(c.genes) for c in gen[0]] == before, "a dropped transfer rotated the donor"


def test_the_ordered_choice_slot_refuses_what_the_family_one_refuses():
    """One validator for all three resolutions, so the words are the same wherever you meet them."""
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=4, seed=1)
    with pytest.raises(ValueError, match="transfer_to must be"):
        simulate_genomes_ordered(sp, transfer=0.1, transfer_to="closest", initial_families=2, seed=1)
    with pytest.raises(ValueError, match="on its own, not a rate"):
        simulate_genomes_ordered(sp, transfer=0.1, initial_families=2, seed=1,
                                 transfer_to=1.0 * mod.DrivenBy("f.tsv", {"a": 2.0}))
    with pytest.raises(ValueError, match="one recipient rule"):
        simulate_genomes_ordered(sp, transfer=0.1, initial_families=2, seed=1,
                                 transfer_to=(Distance(), mod.DrivenBy("f.tsv", {"a": 2.0})))
    with pytest.raises(ValueError, match="silently do nothing"):
        simulate_genomes_ordered(sp, transfer=0.1, initial_families=2, seed=1,
                                 transfer_to=Clades({"A": 0}, Between({("A", "rest"): 1.0})))


def test_gene_order_records_each_chromosome_topology(tmp_path):
    """Topology was recoverable from nothing the run wrote.

    It is load-bearing — it decides where a segmental event stops and which chromosomes may fuse —
    and it is what a rearrangement format's per-chromosome terminator depends on, so a reader handed
    the output directory alone could not export one. Checked on a MIXED karyotype, which is the case
    that made it unrecoverable: on a single-topology run you could at least assume."""
    sp = simulate_species_tree(birth=1.0, n_extant=4, seed=1)
    r = simulate_genomes_ordered(sp.complete_tree, duplication=0.2, loss=0.1, initial_families=6,
                                 chromosomes=2, topology=["circular", "linear"], inversion=0.3,
                                 seed=1)
    r.write(tmp_path, outputs=("gene_order", "initial_genome"))
    rows = (tmp_path / "gene_order.tsv").read_text(encoding="utf-8").splitlines()
    header = rows[0].split("\t")
    assert header[2] == "topology"

    # every row's topology is the one the in-memory chromosome actually has
    seen = set()
    for row in rows[1:]:
        cells = row.split("\t")
        lineage, chrom, topology = node_from_label(cells[0]), int(cells[1]), cells[2]
        actual = {c.id: c.topology for c in r.genomes[lineage]}
        assert topology == actual[chrom], f"{cells[0]} chr{chrom}: wrote {topology!r}"
        seen.add(topology)
    assert seen == {"circular", "linear"}, f"the mixed karyotype should show both, saw {seen}"

    initial = (tmp_path / "initial_genome.tsv").read_text(encoding="utf-8").splitlines()
    assert initial[0].split("\t")[1] == "topology"
    started = {c.id: c.topology for c in r.initial_genome}
    for row in initial[1:]:
        chrom, topology = row.split("\t")[:2]
        assert topology == started[int(chrom)]
