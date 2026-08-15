"""Per-family rates — one named family running at rates of its own (design note §5).

A run's rates apply to every family in it. `family()` declares one family and gives it its own, and
what it leaves out falls back to the run's. These cover the declaration, the two older spellings it
absorbs, the arithmetic reaching both the totals and the copy pick, and the four refusals.
"""

import collections

import pytest

from zombi2.genomes import GeneFamily, family, simulate_genomes_family
from zombi2.params import LogNormal, PerCopy, PerLineage
from zombi2.species import simulate_species_tree


def _tree(n=25, seed=1):
    return simulate_species_tree(birth=1.0, n_extant=n, seed=seed).complete_tree


def _by_family(result, kind):
    """How many events of ``kind`` each family had. A transfer writes two edges, so the recipient
    side is counted to keep one event as one."""
    out = collections.Counter()
    for e in result.edges:
        if e.kind != kind:
            continue
        if kind == "transfer" and e.recipient is None:
            continue
        out[e.family] += 1
    return out


# --- the declaration ------------------------------------------------------------------------------

def test_family_declares_a_name_and_its_rates():
    f = family("IS1", transfer=PerCopy(1.5), loss=0.02, module="mobile")
    assert isinstance(f, GeneFamily)
    assert f.name == "IS1" and f.module == "mobile"
    assert f.written() == {"transfer": PerCopy(1.5), "loss": 0.02}


def test_family_needs_a_name():
    # the whole-genome spec is `genomes.genome(...)`; calling this without a name is that mistake
    with pytest.raises(ValueError, match="genomes.genome"):
        family(duplication=0.2, loss=0.1)


def test_a_family_with_no_rates_is_just_a_name():
    """A declaration with no rates is a name and nothing else — the family runs at the run's rates,
    and is only there so something can refer to it."""
    tree = _tree()
    g = simulate_genomes_family(tree, initial_families=20, duplication=0.1, loss=0.12,
                                families=[family("tox")], seed=5)
    assert set(g.family_names) == {"tox"}
    plain = simulate_genomes_family(tree, initial_families=21, duplication=0.1, loss=0.12, seed=5)
    # the named family is an ordinary one: 20 anonymous + 1 named runs as 21 anonymous would
    assert len(g.events) == len(plain.events)


def test_a_family_names_the_group_it_belongs_to():
    tree = _tree()
    g = simulate_genomes_family(tree, initial_families=10, duplication=0.1, loss=0.1, seed=5,
                                families=[family("nuoA", module="aerobic"),
                                          family("nuoB", module="aerobic")])
    assert set(g.modules) == {"aerobic"}
    assert g.completion("aerobic").history(tree)


@pytest.mark.parametrize("kw, expect", [
    (dict(family_names=["tox"]), "families=\\[family"),
    (dict(origins=[("n3", 0.1)]), "origin="),
    (dict(modules={"m": ["tox"]}), "module="),
])
def test_the_three_older_spellings_name_their_replacement(kw, expect):
    """Three keywords said one thing between them. Removing them without a sentence would leave a
    port with "unexpected keyword argument", which names the mistake and not the fix."""
    tree = _tree(10)
    with pytest.raises(TypeError, match=expect):
        simulate_genomes_family(tree, initial_families=5, loss=0.1, seed=1, **kw)


def test_a_declared_family_can_be_planted():
    """``origin=`` puts a *named* family at a chosen point, which ``origins=`` can only do
    anonymously — today a caller has to work the placed family's id out by arithmetic."""
    tree = _tree()
    lineage = sorted(tree.nodes)[3]
    node = tree.nodes[lineage]
    when = (node.birth_time + node.end_time) / 2
    g = simulate_genomes_family(tree, initial_families=5, duplication=0.05, loss=0.05, seed=2,
                                families=[family("late", origin=(lineage, when))])
    fid = g.family_names["late"]
    born = [e for e in g.edges if e.family == fid and e.kind == "origination"]
    assert len(born) == 1
    assert born[0].time == pytest.approx(when)
    assert born[0].lineage == lineage
    # it is not in the root genome, because it did not exist at the origin
    assert all(c.family != fid for c in g.initial_genome)


# --- the rates actually apply ---------------------------------------------------------------------

def test_a_family_transfers_at_its_own_rate():
    tree = _tree(20)
    g = simulate_genomes_family(tree, initial_families=30, duplication=0.1, loss=0.15,
                                transfer=0.02, seed=3,
                                families=[family("IS1", transfer=PerCopy(1.5), loss=0.02)])
    moved = _by_family(g, "transfer")
    is1 = g.family_names["IS1"]
    others = sum(v for k, v in moved.items() if k != is1)
    # 1.5 against 0.02 per copy, so the element should dominate the transfers outright
    assert moved[is1] > 10 * max(others, 1), f"IS1 {moved[is1]}, the other 30 families {others}"


def test_a_family_can_be_almost_never_lost():
    tree = _tree(20)
    g = simulate_genomes_family(tree, initial_families=40, duplication=0.05, loss=0.4, seed=4,
                                families=[family("core", loss=0.001)])
    core = g.family_names["core"]
    lost = _by_family(g, "loss")
    assert lost[core] == 0
    assert sum(v for k, v in lost.items() if k != core) > 50
    # and it survives to every tip, where the rest of the genome is being stripped
    assert all(any(c.family == core for c in gen) for gen in g.genomes.values())


def test_the_copy_pick_follows_the_same_weights():
    """The totals and the pick have to agree. A family lost 100x faster than the rest must take that
    share of the losses, not merely raise the total."""
    tree = _tree(18)
    g = simulate_genomes_family(tree, initial_families=20, duplication=0.25, loss=0.01, seed=6,
                                families=[family("fragile", loss=1.0)])
    fragile = g.family_names["fragile"]
    lost = _by_family(g, "loss")
    assert lost[fragile] > sum(v for k, v in lost.items() if k != fragile)


def test_a_written_rate_ignores_the_runs_per_family_draw():
    """A written rate is the rate. A draw meant to vary the run's rate among families must not
    multiply it, or the number written would not be the number run."""
    tree = _tree(15)
    spread = PerCopy(0.3).varying_among("families", LogNormal(0.0, 1.2))
    runs = [simulate_genomes_family(tree, initial_families=8, duplication=0.05, loss=spread,
                                    seed=s, families=[family("steady", loss=0.0)])
            for s in (1, 2, 3)]
    for g in runs:
        steady = g.family_names["steady"]
        assert _by_family(g, "loss")[steady] == 0, "a loss rate written as 0 lost copies anyway"


# --- what is refused ------------------------------------------------------------------------------

def test_a_family_rate_takes_no_verb_yet():
    tree = _tree(10)
    with pytest.raises(ValueError, match="does not take yet"):
        simulate_genomes_family(tree, initial_families=5, loss=0.1, seed=1,
                                families=[family("x", loss=PerCopy(0.1).changing_at({0: 1.0, 2: 0.5}))])


def test_a_family_rate_is_per_copy():
    tree = _tree(10)
    with pytest.raises(ValueError, match="per copy"):
        simulate_genomes_family(tree, initial_families=5, loss=0.1, seed=1,
                                families=[family("x", loss=PerLineage(0.1))])


def test_a_family_rate_composes_with_a_driven_run_rate():
    """The driver scales the run's rate. A family that states its own rate has no run rate in it to
    scale, so the two are added rather than multiplied — and the family's rate is the same number in
    caves and on the surface, which is what "the rate itself" means."""
    from zombi2.traits import simulate_discrete

    tree = _tree(20)
    habitat = simulate_discrete(tree, states=["cave", "surface"], switch=0.2, seed=1)
    driven = PerCopy(0.2).scaled_by(habitat, {"cave": 6.0, "surface": 1.0})
    g = simulate_genomes_family(tree, initial_families=20, duplication=0.1, loss=driven, seed=1,
                                families=[family("steady", loss=0.0)])
    steady = g.family_names["steady"]
    # the run's loss is driven and fires; the family that wrote 0.0 is never lost, in either habitat
    lost = _by_family(g, "loss")
    assert sum(v for k, v in lost.items() if k != steady) > 20
    assert lost[steady] == 0


def test_the_per_family_engine_refuses_a_family_rate():
    tree = _tree(10)
    with pytest.raises(ValueError, match="per-family engine"):
        simulate_genomes_family(tree, initial_families=5, loss=0.1, seed=1, parallel=True,
                                families=[family("x", loss=0.05)])


def test_names_must_be_unique():
    tree = _tree(10)
    with pytest.raises(ValueError, match="unique"):
        simulate_genomes_family(tree, initial_families=5, loss=0.1, seed=1,
                                families=[family("x"), family("x", loss=0.05)])


def test_families_takes_declarations_not_names():
    tree = _tree(10)
    with pytest.raises(TypeError, match=r"family\('IS1'\)"):
        simulate_genomes_family(tree, initial_families=5, loss=0.1, seed=1, families=["IS1"])
