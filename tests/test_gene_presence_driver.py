"""Gene presence as a conditioning driver — a genome making a trait's rate faster.

The other direction of a relation ZOMBI2 already had. A trait could drive a genome rate; now a gene
family's presence can drive anything a driver reaches, a trait's switch rate included.

The two tests worth reading are the cross-checks. The trajectory is derived from the family's gene
tree, so both of them grade it against something derived a different way: `has_family`, which reads
the genome at a node, and the recorded loss events, which say when a copy actually went.
"""

from __future__ import annotations

import pytest

from zombi2.genomes import simulate_genomes_family, simulate_genomes_ordered
from zombi2.rates import modifiers as mod
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_discrete


def _run(loss=0.15, seed=9, n_extant=40):
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=n_extant, seed=4)
    g = simulate_genomes_family(sp.complete_tree, initial_families=20, family_names=["tox"],
                                duplication=0.1, loss=loss, seed=seed)
    return sp.complete_tree, g


def test_the_trajectory_agrees_with_the_genome_at_every_branch_end():
    """Graded against `has_family`, which reads a node's genome rather than the family's gene tree —
    so agreement is two derivations of the same fact meeting, not one of them checking itself."""
    tree, g = _run()
    traj = g.presence("tox").as_driver_trajectory(tree)
    for node in tree.nodes.values():
        just_before_the_end = node.end_time - 1e-9
        if just_before_the_end < node.birth_time:
            continue                                  # a zero-length branch has no interior
        expected = "present" if g.has_family(node.id, "tox") else "absent"
        assert traj.value(node.id, just_before_the_end) == expected


def test_the_family_goes_absent_exactly_when_a_copy_was_lost():
    """The mid-branch case, and the one a per-branch snapshot cannot see. Every instant the family
    turns absent has to be a recorded loss — nothing else takes the last copy away mid-branch."""
    tree, g = _run()
    traj = g.presence("tox").as_driver_trajectory(tree)
    losses = {round(e.time, 9) for e in g.events
              if e.kind == "loss" and e.family == g.family_names["tox"]}
    changes = [(node_id, t) for node_id, starts in traj._starts.items()
               for t, state in zip(starts, traj._states[node_id])
               if state == "absent" and t > tree.nodes[node_id].birth_time]
    assert changes, "this run is meant to lose the family somewhere mid-branch"
    for node_id, t in changes:
        assert round(t, 9) in losses, f"lineage {node_id} went absent at {t}, which is not a loss"


def test_every_lineage_gets_an_answer():
    """A driver is asked about whatever branch the engine is on, so a lineage that never held the
    family has to answer 'absent' rather than raise."""
    tree, g = _run()
    traj = g.presence("tox").as_driver_trajectory(tree)
    assert set(traj._starts) == set(tree.nodes)
    for node in tree.nodes.values():
        assert traj.value(node.id, node.birth_time) in ("present", "absent")


def test_a_gene_can_drive_a_trait():
    tree, g = _run()
    kw = dict(states=["harmless", "pathogenic"], start="harmless", seed=1)
    plain = simulate_discrete(tree, switch=0.1, **kw)
    driven = simulate_discrete(tree, switch=0.1 * mod.DrivenBy(
        g.presence("tox"), {"present": 8.0, "absent": 1.0}), **kw)
    switches = lambda r: sum(1 for e in r.events if e.kind == "on_branch")
    assert switches(driven) > switches(plain)


def test_a_driver_that_never_changes_leaves_the_run_byte_identical():
    """The same contract a trait driver has: what perturbs the stream is the *breakpoints* a driver
    adds to the Gillespie horizon, so a family present on every branch for all of its life adds none.
    """
    tree, g = _run(loss=0.0)
    traj = g.presence("tox").as_driver_trajectory(tree)
    assert traj.states() == {"present"}
    kw = dict(states=["x", "y"], start="x", seed=2)
    plain = simulate_discrete(tree, switch=0.2, **kw)
    driven = simulate_discrete(tree, switch=0.2 * mod.DrivenBy(g.presence("tox"), {"present": 1.0}),
                               **kw)
    assert plain.values == driven.values


def test_ordered_runs_carry_the_same_signal():
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=12, seed=4)
    g = simulate_genomes_ordered(sp.complete_tree, initial_families=12, family_names=["tox"],
                                 duplication=0.1, loss=0.2, inversion=0.2, seed=7)
    traj = g.presence("tox").as_driver_trajectory(sp.complete_tree)
    assert traj.states() <= {"present", "absent"}
    for node in sp.complete_tree.extant_leaves():
        expected = "present" if g.has_family(node.id, "tox") else "absent"
        assert traj.value(node.id, node.end_time - 1e-9) == expected


def test_it_refuses_what_it_cannot_answer():
    tree, g = _run()
    with pytest.raises(KeyError, match="no named family"):
        g.presence("nope")
    other = simulate_species_tree(birth=1.0, death=0.2, n_extant=9, seed=77).complete_tree
    with pytest.raises(ValueError, match="different species trees"):
        g.presence("tox").as_driver_trajectory(other)


# --- modules: a named group of families, and how much of it a lineage carries ------------------

_FLG = [f"flg{i}" for i in range(6)]


def _module_run(loss=0.2, seed=9, n_extant=60):
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=n_extant, seed=4)
    g = simulate_genomes_family(sp.complete_tree, initial_families=20, family_names=_FLG,
                                modules={"flagellum": _FLG}, duplication=0.05, loss=loss, seed=seed)
    return sp.complete_tree, g


def test_completion_is_the_fraction_of_the_module_a_lineage_carries():
    """Graded against a naive per-tip count of `has_family`, which reads each node's genome rather
    than the families' gene trees — two derivations of the same fact meeting."""
    tree, g = _module_run()
    traj = g.completion("flagellum").as_driver_trajectory(tree)
    for node in tree.nodes.values():
        just_before = node.end_time - 1e-9
        if just_before < node.birth_time:
            continue
        naive = sum(g.has_family(node.id, m) for m in _FLG) / len(_FLG)
        assert traj.value(node.id, just_before) == naive


def test_completion_takes_every_step_between_none_and_all():
    tree, g = _module_run()
    states = g.completion("flagellum").as_driver_trajectory(tree).states()
    assert states == {k / len(_FLG) for k in range(len(_FLG) + 1)}


def test_a_module_of_one_family_is_that_family_s_presence():
    """The two readers share a span builder, so this is the seam between them: a one-family module
    must reduce to presence, 1.0 where present and 0.0 where absent."""
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=30, seed=4)
    tree = sp.complete_tree
    g = simulate_genomes_family(tree, initial_families=15, family_names=["tox"],
                                modules={"just_tox": ["tox"]}, duplication=0.1, loss=0.15, seed=9)
    one = g.completion("just_tox").as_driver_trajectory(tree)
    pres = g.presence("tox").as_driver_trajectory(tree)
    for node in tree.nodes.values():
        t = node.end_time - 1e-9
        if t < node.birth_time:
            continue
        assert one.value(node.id, t) == (1.0 if pres.value(node.id, t) == "present" else 0.0)


def test_every_lineage_gets_a_completion():
    tree, g = _module_run()
    traj = g.completion("flagellum").as_driver_trajectory(tree)
    assert set(traj._starts) == set(tree.nodes)


def test_completion_drives_a_trait():
    from zombi2.rates.mapping import Curve
    tree, g = _module_run()
    kw = dict(states=["sessile", "motile"], start="sessile", seed=2)
    plain = simulate_discrete(tree, switch=0.05, **kw)
    driven = simulate_discrete(tree, switch=0.05 * mod.DrivenBy(
        g.completion("flagellum"), Curve(lambda f: 0.05 + 30.0 * f ** 4)), **kw)
    switches = lambda r: sum(1 for e in r.events if e.kind == "on_branch")
    assert switches(driven) != switches(plain)


def test_the_module_history_sums_to_the_branch_lengths():
    tree, g = _module_run()
    for node_id, segs in g.completion("flagellum").history(tree).items():
        node = tree.nodes[node_id]
        assert abs(sum(d for _, d in segs) - (node.end_time - node.birth_time)) < 1e-9


def test_a_module_refuses_what_it_cannot_mean():
    tree, g = _module_run()
    with pytest.raises(KeyError, match="no module"):
        g.completion("nope")
    kw = dict(initial_families=4, family_names=["a", "b"], seed=1)
    for bad, fragment in (({"m": []}, "no families in it"),
                          ({"m": ["a", "z"]}, "not declared families"),
                          ({"m": ["a", "a"]}, "names a family twice"),
                          ({"": ["a"]}, "non-empty name")):
        with pytest.raises(ValueError, match=fragment):
            simulate_genomes_family(tree, modules=bad, **kw)
    with pytest.raises(TypeError, match="must be a dict"):
        simulate_genomes_family(tree, modules=["a"], **kw)
