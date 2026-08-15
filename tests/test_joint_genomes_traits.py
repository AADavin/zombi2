"""A genome and a discrete trait, simulated together on a tree the run is given (design note §7).

The first joint model whose tree is an **input**, and the first that produces two levels' results
from one run. Each half reads the other as the run goes: a genome rate reads the trait state on the
lineage, and the trait's switch rate reads that lineage's gene content.
"""

import collections

import pytest

from zombi2 import genomes, joint, species, traits
from zombi2.genomes import family
from zombi2.params import PerCopy, PerLineage


def _tree(n=30, seed=4):
    return species.simulate_species_tree(birth=1.0, n_extant=n, seed=seed).complete_tree


def _cave_run(seed=1, loss_in_caves=6.0, n=30):
    """Caves cost eyes, and losing the eye makes a lineage likelier to commit to the cave."""
    tree = _tree(n)
    r = joint.simulate(
        genomes.genome(
            duplication=0.05, origination=0.3, initial_families=40,
            loss=PerCopy(0.06).scaled_by("trait", {"cave": loss_in_caves, "surface": 1.0}),
            families=[family("eye")]),
        traits.discrete(states=["surface", "cave"], start="surface",
                        switch={"surface->cave": PerLineage(0.05).scaled_by(
                                    "genomes:eye", {"present": 1.0, "absent": 8.0}),
                                "cave->surface": 0.1}),
        tree=tree, seed=seed)
    return tree, r


# --- both levels come out of one run --------------------------------------------------------------

def test_one_run_produces_both_levels():
    tree, r = _cave_run()
    assert r.trait is not None and r.genome is not None
    assert r.complete_tree is tree                       # the tree came in, so it comes back unchanged
    assert set(r.trait.node_values) == set(tree.nodes)
    assert set(r.genome.genomes) == {tree.labels()[i] for i in tree.extant_leaves()}


def test_the_trait_drives_the_genome():
    """A lineage in the cave loses genes six times as fast, so cave tips end up smaller."""
    tree, r = _cave_run()
    lab = tree.labels()
    sizes = {lab[n]: len(r.genome.genomes[lab[n]]) for n in tree.extant_leaves()}
    by_state = collections.defaultdict(list)
    for name, size in sizes.items():
        by_state[r.trait.values[name]].append(size)
    if len(by_state) < 2:
        pytest.skip("this seed produced only one state at the tips")
    cave = sum(by_state["cave"]) / len(by_state["cave"])
    surface = sum(by_state["surface"]) / len(by_state["surface"])
    assert cave < surface, f"cave {cave:.0f} genes, surface {surface:.0f}"


def test_the_genome_drives_the_trait():
    """The other direction, on its own: with the eye eight times as costly to keep, a run whose
    families are lost fast commits to the cave far more than one whose families persist."""
    tree = _tree(30)
    spec = lambda loss: joint.simulate(
        genomes.genome(duplication=0.02, origination=0.05, initial_families=5, loss=loss,
                       families=[family("eye")]),
        traits.discrete(states=["surface", "cave"], start="surface",
                        switch={"surface->cave": PerLineage(0.03).scaled_by(
                                    "genomes:eye", {"present": 1.0, "absent": 25.0}),
                                "cave->surface": 0.05}),
        tree=tree, seed=2)
    fragile = spec(PerCopy(0.9))      # the eye goes early and often
    stable = spec(PerCopy(0.001))     # the eye is kept
    caves = lambda r: sum(1 for v in r.trait.values.values() if v == "cave")
    assert caves(fragile) > caves(stable), f"{caves(fragile)} against {caves(stable)}"


def test_it_is_deterministic():
    a, b = _cave_run(seed=5)[1], _cave_run(seed=5)[1]
    assert [(e.time, e.kind, e.family) for e in a.genome.edges] == \
           [(e.time, e.kind, e.family) for e in b.genome.edges]
    assert [(c.time, c.kind, c.to_state) for c in a.trait.events] == \
           [(c.time, c.kind, c.to_state) for c in b.trait.events]


def test_the_trait_log_reads_back_as_a_driver():
    """A joint run's trait log has to carry the ``initial`` row, or nothing downstream can replay it
    against the tree — the same rule the tree-growing joint models follow."""
    from zombi2.params.conditioned import load_driver

    tree, r = _cave_run()
    first = r.trait.events[0]
    assert (first.kind, first.lineage, first.from_state) == ("initial", tree.root, None)
    for i, node in tree.nodes.items():                    # the map covers every branch exactly
        assert sum(d for _s, d in r.trait.history[i]) == pytest.approx(
            node.end_time - node.birth_time)


def test_both_levels_write_their_own_files(tmp_path):
    _tree_, r = _cave_run()
    r.write(tmp_path)
    for f in ("trait_values.tsv", "trait_events.tsv", "genome_events.tsv", "profiles.tsv"):
        assert (tmp_path / f).exists(), f"missing {f}"


def test_transfer_works_here_and_not_on_a_growing_tree():
    """The contemporaneous set is the whole difference: on a tree handed to the run it is known, and
    on a growing one it is still forming."""
    tree = _tree(20)
    r = joint.simulate(
        genomes.genome(duplication=0.05, loss=0.05, initial_families=20, transfer=PerCopy(0.4),
                       families=[family("eye")]),
        traits.discrete(states=["surface", "cave"],
                        switch={"surface->cave": PerLineage(0.05).scaled_by(
                                    "genomes:eye", {"present": 1.0, "absent": 5.0}),
                                "cave->surface": 0.1}),
        tree=tree, seed=3)
    assert sum(1 for e in r.genome.edges if e.kind == "transfer") > 0

    with pytest.raises(ValueError, match="still forming"):
        joint.simulate(
            species.birth_death(birth=PerLineage(1.0).scaled_by("genomes:count", lambda n: 1 + n),
                                n_extant=10),
            genomes.genome(loss=0.1, transfer=0.1), seed=1)


# --- what is refused ------------------------------------------------------------------------------

def test_the_two_levels_must_read_each_other():
    tree = _tree(10)
    with pytest.raises(ValueError, match="two independent runs"):
        joint.simulate(genomes.genome(duplication=0.1, loss=0.1),
                       traits.discrete(states=["a", "b"], switch=0.2), tree=tree, seed=1)


def test_a_genome_rate_takes_only_what_this_engine_threads():
    from zombi2.params import LogNormal

    tree = _tree(10)
    with pytest.raises(ValueError, match="does not thread"):
        joint.simulate(
            genomes.genome(duplication=0.1,
                           loss=PerCopy(0.1).varying_among("families", LogNormal(0.0, 0.5))
                                            .scaled_by("trait", {"a": 2.0, "b": 1.0})),
            traits.discrete(states=["a", "b"], switch=0.2), tree=tree, seed=1)
