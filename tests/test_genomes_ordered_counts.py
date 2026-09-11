"""The ordered engine's running counts of genes per family (issue #433).

The engine keeps, for every living lineage, how many genes of each family it holds, and changes that
count wherever it adds or removes genes. The family-size cap and every rule that reads gene content
look a family up there instead of walking a genome. A count that drifts would change a run without
any error, so these tests turn on the engine's check, which recounts every genome before each step of
the loop and once after it. Each setup first requires the events it is there to exercise.
"""

import collections

import pytest

import zombi2.genomes.ordered as ordered
from zombi2.genomes import family, simulate_genomes_ordered
from zombi2.params import Fixed, Geometric, PerCopy, Recipients
from zombi2.species import simulate_species_tree

W = {"present": 3.0, "absent": 1.0}


@pytest.fixture(scope="module")
def tree():
    return simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=4).complete_tree


@pytest.fixture
def checked(monkeypatch):
    monkeypatch.setattr(ordered, "_CHECK_COUNTS", True)


def _run(tree, **kw):
    base = dict(duplication=0.2, transfer=0.1, loss=0.25, origination=0.3, initial_families=12, seed=3)
    return simulate_genomes_ordered(tree, **{**base, **kw})


def _kinds(result):
    return collections.Counter(e.kind for e in result.edges)


def test_duplication_loss_and_origination_over_segments(tree, checked):
    g = _run(tree, duplication_extent=Geometric(3.0), loss_extent=Geometric(3.0), max_family_size=4)
    kinds = _kinds(g)
    assert kinds["duplication"] and kinds["loss"] and kinds["origination"]


def test_transfers_that_replace_and_return_to_the_donor(tree, checked):
    g = _run(tree, transfer=0.6, transfer_extent=Geometric(2.0), replacement=True, self_transfer=True,
             max_family_size=3)
    arrivals = [e for e in g.edges if e.kind == "transfer" and e.recipient is not None]
    assert arrivals
    assert any(e.replaced is not None for e in arrivals)
    assert any(e.recipient == e.donor for e in arrivals)


def test_chromosome_events_on_linear_chromosomes(tree, checked):
    g = _run(tree, chromosomes=3, topology="linear", fission=0.2, fusion=0.2, chromosome_origination=0.1,
             chromosome_loss=0.2, translocation=0.3, transposition=0.3, inversion=0.3)
    kinds = collections.Counter(e.kind for e in g.chromosome_events)
    assert kinds["fission"] and kinds["fusion"] and kinds["origination"] and kinds["loss"]
    assert g.rearrangements


def test_placed_families_speciation_and_extinction(tree, checked):
    late = sorted(n for n in tree.nodes if n != tree.root)[6]
    g = _run(tree, families=[family("late", origin=late)])
    assert any(tree.nodes[n].fate == "extinct" for n in tree.nodes)
    assert any(e.kind == "origination" and e.lineage == late and e.family == g.family_names["late"]
               for e in g.edges)


def test_rules_that_read_gene_content(tree, checked):
    clade = sorted(n for n in tree.nodes if tree.nodes[n].parent == tree.root)[0]
    g = _run(tree, joint=True, transfer=0.4, transfer_extent=Fixed(3), loss_extent=Fixed(2),
             loss=PerCopy(0.3).scaled_by("genomes:count", lambda n: min(2.0, n / 20.0)),
             transfer_to=Recipients().weighted_by("genomes:A", W),
             families=[family("A", origin=clade),
                       family("B", duplication=PerCopy(0.4).scaled_by("genomes:A", W),
                              transfer_to=Recipients().weighted_by("genomes:module:m", lambda s: 1.0 + s)),
                       family("x", module="m"), family("y", module="m", origin=clade)])
    kinds = _kinds(g)
    assert kinds["loss"] and kinds["transfer"] and kinds["duplication"]


def test_the_check_catches_a_count_that_drifts(tree, checked, monkeypatch):
    """A duplication that forgets to count its copies is caught at the next step."""
    monkeypatch.setattr(ordered._GeneCounts, "added_all", lambda self, k, families: None)
    with pytest.raises(AssertionError, match="running counts differ"):
        _run(tree, transfer=0.0)


def test_the_cap_still_binds(tree):
    g = _run(tree, duplication=0.8, loss=0.05, duplication_extent=Geometric(2.0), max_family_size=3)
    per_genome = [collections.Counter(gene.family for chrom in genome for gene in chrom.genes)
                  for genome in g.node_genomes.values()]
    assert max(max(c.values(), default=0) for c in per_genome) == 3
