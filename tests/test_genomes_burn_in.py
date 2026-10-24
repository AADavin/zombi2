"""A burn-in before the root (issue #473): the genome process runs on one lineage before the species
tree begins, so the root holds families of many ages. The two closed forms it reaches are checked
in ``test_validation.py``; here are the mechanics and the refusals."""

import collections

import pytest

from zombi2.genomes import family, group, simulate_genomes_family, simulate_genomes_ordered
from zombi2.params import PerCopy, PerLineage, LogNormal, Recipients
from zombi2.species import simulate_species_tree
from zombi2.tree import Node, Tree

LEAF = Tree({0: Node(0, None, 0.0, 1.0, (), "extant")}, 0)


@pytest.fixture(scope="module")
def tree():
    return simulate_species_tree(birth=1.0, death=0.3, n_extant=12, seed=2).complete_tree


def test_the_burn_in_runs_before_the_root_and_the_tree_starts_from_what_it_leaves(tree):
    r = simulate_genomes_family(tree, origination=20.0, duplication=0.2, loss=1.0, transfer=0.2,
                                initial_families=0, burn_in=10.0, seed=1)
    before = [e for e in r.edges if e.time < 0]
    assert before and min(e.time for e in before) >= -10.0
    assert all(e.lineage == tree.root for e in before)
    assert all(e.kind in ("origination", "duplication", "loss") for e in before)   # no transfer
    families_at_root = {c.family for c in r.initial_genome}
    born_before = {e.family for e in before if e.kind == "origination"}
    assert families_at_root and families_at_root <= born_before      # every root family came from it
    assert len(born_before) > len(families_at_root)                  # and some were lost on the way
    # the gene tree of a root family reaches above the root
    fid = next(iter(families_at_root))
    assert any(e.time < 0 for e in r.edges if e.family == fid)


def test_initial_families_and_declared_families_add_to_the_burn_in(tree):
    r = simulate_genomes_family(tree, origination=5.0, duplication=0.1, loss=1.0, initial_families=7,
                                families=[family("toxin")], burn_in=5.0, seed=3)
    at_root = collections.Counter(e.time for e in r.edges if e.kind == "origination"
                                  and e.time <= tree.nodes[tree.root].birth_time)
    assert at_root[tree.nodes[tree.root].birth_time] == 8          # 7 anonymous + the declared one
    assert "toxin" in r.family_names


def test_a_run_without_a_burn_in_is_the_run_it_was(tree):
    a = simulate_genomes_family(tree, origination=1.0, duplication=0.2, loss=0.3, initial_families=20, seed=4)
    b = simulate_genomes_family(tree, origination=1.0, duplication=0.2, loss=0.3, initial_families=20, seed=4,
                                burn_in=None)
    assert [(e.time, e.kind, e.family, e.copy) for e in a.edges] == \
        [(e.time, e.kind, e.family, e.copy) for e in b.edges]


def test_the_ordered_burn_in_rearranges_and_duplicates_segments_before_the_root(tree):
    r = simulate_genomes_ordered(tree, origination=20.0, duplication=0.2, loss=1.0, inversion=0.3,
                                 transposition=0.1, translocation=0.1, chromosomes=2,
                                 duplication_extent=2, initial_families=0, burn_in=10.0, seed=1,
                                 fission=0.2, fusion=0.2)
    assert sum(len(c.genes) for c in r.initial_genome) > 0
    assert any(x.time < 0 for x in r.rearrangements)
    assert all(e.time >= 0 for e in r.chromosome_events if e.kind != "initial")   # not before the root
    assert len(r.initial_genome) == 2
    # a declared group still starts at the root, whole
    g = simulate_genomes_ordered(tree, origination=5.0, duplication=0.2, loss=1.0, initial_families=0,
                                 burn_in=5.0, seed=1, families=[family("A"), family("B")],
                                 groups=[group("op", ["A", "B"])])
    seq = [x.family for c in g.initial_genome for x in c.genes]
    i = seq.index(g.family_names["A"])
    assert seq[i + 1] == g.family_names["B"]


def test_a_per_family_draw_weights_the_burn_in_at_the_family_resolution():
    r = simulate_genomes_family(LEAF, origination=10.0, duplication=PerCopy(0.2).varying_among(
        "families", LogNormal(0.0, 0.5)), loss=1.0, initial_families=0, burn_in=10.0, seed=1,
        max_family_size=None)
    assert r.initial_genome
    assert set(r.family_multipliers) >= {c.family for c in r.initial_genome}


def test_the_burn_in_is_refused_where_it_has_no_meaning(tree):
    base = dict(origination=1.0, initial_families=0, seed=1)
    with pytest.raises(ValueError, match="below loss"):
        simulate_genomes_family(tree, duplication=1.0, loss=0.5, burn_in=5.0, **base)
    with pytest.raises(ValueError, match="positive time"):
        simulate_genomes_family(tree, duplication=0.2, loss=1.0, burn_in=-1.0, **base)
    with pytest.raises(ValueError, match="joint"):
        simulate_genomes_family(tree, duplication=0.2, loss=1.0, burn_in=5.0, families=[family("A")],
                                transfer=0.1, transfer_to=Recipients().weighted_by(
                                    "genomes:A", {"present": 2.0, "absent": 1.0}), joint=True, **base)
    from zombi2.traits import simulate_discrete
    habitat = simulate_discrete(tree, states=["x", "y"], switch=0.5, seed=1)
    with pytest.raises(ValueError, match="driven"):
        simulate_genomes_family(tree, duplication=0.2, burn_in=5.0,
                                loss=PerCopy(1.0).scaled_by(habitat, {"x": 1.0, "y": 0.5}), **base)
    with pytest.raises(ValueError, match="lineages"):
        simulate_genomes_family(tree, duplication=0.2, loss=PerCopy(1.0).varying_among(
            "lineages", LogNormal(0.0, 0.3)), burn_in=5.0, **base)
    with pytest.raises(ValueError, match="PerLineage"):
        simulate_genomes_family(tree, duplication=0.2, loss=PerLineage(1.0), burn_in=5.0, **base)
    with pytest.raises(ValueError, match="per-family engine"):
        simulate_genomes_family(tree, duplication=0.2, loss=1.0, burn_in=5.0, parallel=True, **base)
    with pytest.raises(ValueError, match="not built at the ordered"):
        simulate_genomes_ordered(tree, duplication=PerCopy(0.2).varying_among(
            "families", LogNormal(0.0, 0.5)), loss=1.0, burn_in=5.0, **base)
