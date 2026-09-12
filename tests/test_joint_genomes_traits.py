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


# --- a declared family's own rates (issue #437) ---------------------------------------------------
# On a given tree the genome is a target as well as a driver, so a declared family can carry its own
# duplication, transfer and loss, and those rates can read the trait. The checks use factors of 0, so
# each is exact, and each first requires the event it is about to have happened.

def _state_at(r, tree, lineage, time):
    """The trait state on ``lineage`` at ``time``, read back from the run's own trait history."""
    at = tree.nodes[lineage].birth_time
    state = None
    for state, duration in r.trait.history[lineage]:
        if time < at + duration:
            return state
        at += duration
    return state


def test_the_run_in_issue_437_is_accepted():
    """Alyssa Henderson's run, as she wrote it: a family whose loss reads the trait, on a tree passed
    in. It used to be refused with a message about the tree being simulated."""
    tree = species.simulate_species_tree(birth=1.0, n_extant=30, seed=1).complete_tree
    r = joint.simulate(
        genomes.genome(
            initial_families=10, duplication=0.05, origination=0.0,
            loss=PerCopy(0.30),
            families=[family("A", loss=PerCopy(0.30).scaled_by(
                "trait", {"on": 6.0, "off": 1.0}))]),
        traits.discrete(states=["off", "on"], start="off",
                        switch={"off->on": 0.1, "on->off": 0.1}),
        tree=tree, seed=1)
    assert "A" in r.genome.family_names


def test_two_families_are_lost_in_opposite_states():
    """The benchmark the issue asks for: one family lost only where the trait is on, another only where
    it is off."""
    tree = _tree(40, seed=3)
    r = joint.simulate(
        genomes.genome(duplication=0.4, origination=0.0, initial_families=10, loss=PerCopy(0.2),
                       families=[family("on_only", loss=PerCopy(0.6).scaled_by(
                                     "trait", {"on": 1.0, "off": 0.0})),
                                 family("off_only", loss=PerCopy(0.6).scaled_by(
                                     "trait", {"on": 0.0, "off": 1.0}))]),
        traits.discrete(states=["off", "on"], start="off", switch={"off->on": 0.4, "on->off": 0.4}),
        tree=tree, seed=4)
    ids = r.genome.family_names
    for name, allowed in (("on_only", "on"), ("off_only", "off")):
        losses = [e for e in r.genome.edges if e.kind == "loss" and e.family == ids[name]]
        assert losses, f"{name} was never lost"
        assert {_state_at(r, tree, e.lineage, e.time) for e in losses} == {allowed}


def test_a_family_whose_own_rates_are_zero_never_duplicates_or_leaves():
    r = joint.simulate(
        genomes.genome(duplication=0.3, transfer=0.3, loss=0.05, initial_families=15,
                       families=[family("still", duplication=0.0, transfer=0.0), family("eye")]),
        traits.discrete(states=["a", "b"],
                        switch={"a->b": PerLineage(0.2).scaled_by("genomes:eye",
                                                                  {"present": 1.0, "absent": 3.0}),
                                "b->a": 0.2}),
        tree=_tree(25), seed=6)
    still = r.genome.family_names["still"]
    kinds = collections.Counter(e.kind for e in r.genome.edges)
    assert kinds["duplication"] and kinds["transfer"]
    assert not [e for e in r.genome.edges if e.kind in ("duplication", "transfer") and e.family == still]


def test_a_family_rate_can_read_another_family():
    """B is never lost where A is present, and A is never lost, so B is never lost at all."""
    r = joint.simulate(
        genomes.genome(loss=0.3, initial_families=10,
                       families=[family("A", loss=0.0),
                                 family("B", loss=PerCopy(0.8).scaled_by(
                                     "genomes:A", {"present": 0.0, "absent": 1.0}))]),
        traits.discrete(states=["a", "b"],
                        switch={"a->b": PerLineage(0.1).scaled_by("genomes:A",
                                                                  {"present": 1.0, "absent": 2.0}),
                                "b->a": 0.1}),
        tree=_tree(20), seed=2)
    b = r.genome.family_names["B"]
    assert [e for e in r.genome.edges if e.kind == "loss"]
    assert not [e for e in r.genome.edges if e.kind == "loss" and e.family == b]


def test_a_family_rate_on_a_schedule_stops_where_the_schedule_does():
    """The race has to stop at a family's own breakpoint: stepping over it would lose this family
    after 0.5 at the rate it had before."""
    r = joint.simulate(
        genomes.genome(duplication=0.5, loss=0.05, initial_families=8,
                       families=[family("early", loss=PerCopy(5.0).changing_at({0.0: 1.0, 0.5: 0.0})
                                                  .scaled_by("trait", {"a": 1.0, "b": 1.0}))]),
        traits.discrete(states=["a", "b"], switch=0.3),
        tree=_tree(25), seed=3)
    early = r.genome.family_names["early"]
    losses = [e.time for e in r.genome.edges if e.kind == "loss" and e.family == early]
    assert losses and max(losses) < 0.5


def test_family_rates_keep_a_seed_reproducible():
    run = lambda: joint.simulate(
        genomes.genome(duplication=0.2, loss=PerCopy(0.2), initial_families=10,
                       families=[family("A", loss=PerCopy(0.5).scaled_by("trait", {"on": 4.0, "off": 1.0}),
                                        duplication=0.1)]),
        traits.discrete(states=["off", "on"], switch={"off->on": 0.2, "on->off": 0.2}),
        tree=_tree(20), seed=8)
    a, b = run(), run()
    assert [(e.time, e.kind, e.family) for e in a.genome.edges] == \
           [(e.time, e.kind, e.family) for e in b.genome.edges]


def test_genome_accepts_a_family_with_rates():
    """Whether a family's rates are read is decided by the run, which knows whether the tree is given."""
    spec = genomes.genome(families=[family("A", loss=0.3)])
    assert spec.family_names == ("A",)


def test_a_family_origin_or_transfer_to_is_refused_on_a_given_tree():
    from zombi2.params import Recipients

    tree = _tree(10)
    trait = traits.discrete(states=["a", "b"], switch=0.2)
    for declared, what in ((family("A", loss=PerCopy(0.2).scaled_by("trait", {"a": 2.0, "b": 1.0}),
                                   origin=(tree.root, None)), "an origin"),
                           (family("A", loss=PerCopy(0.2).scaled_by("trait", {"a": 2.0, "b": 1.0}),
                                   transfer_to=Recipients().weighted_by("genomes:A",
                                                                        {"present": 2.0, "absent": 1.0})),
                            "a transfer_to")):
        with pytest.raises(ValueError, match=f"sets {what}, which a joint run does not read"):
            joint.simulate(genomes.genome(loss=0.1, families=[declared]), trait, tree=tree, seed=1)


def test_a_growing_tree_still_reads_only_a_family_name():
    with pytest.raises(ValueError, match="does not read: the gene content drives speciation"):
        joint.simulate(
            species.birth_death(birth=PerLineage(1.0).scaled_by("genomes:A", {"present": 2.0, "absent": 1.0}),
                                n_extant=10),
            genomes.genome(loss=0.1, families=[family("A", loss=0.3)]), seed=1)


def test_a_family_rate_reading_an_undeclared_family_is_refused():
    with pytest.raises(ValueError, match="family 'A''s loss reads 'genomes:Z'"):
        joint.simulate(
            genomes.genome(loss=0.1, families=[family("A", loss=PerCopy(0.2).scaled_by(
                "genomes:Z", {"present": 2.0, "absent": 1.0}))]),
            traits.discrete(states=["a", "b"],
                            switch={"a->b": PerLineage(0.1).scaled_by("genomes:A",
                                                                      {"present": 1.0, "absent": 2.0}),
                                    "b->a": 0.1}),
            tree=_tree(10), seed=1)


def test_a_family_rate_needs_a_per_copy_run_rate():
    with pytest.raises(ValueError, match="counted per copy"):
        joint.simulate(
            genomes.genome(loss=PerLineage(0.1),
                           families=[family("A", loss=PerCopy(0.2).scaled_by("trait", {"a": 2.0, "b": 1.0}))]),
            traits.discrete(states=["a", "b"], switch=0.2),
            tree=_tree(10), seed=1)
