"""A module's completeness as a live driver in joint runs (issue #425).

``"genomes:module:<group>"`` reads how complete a declared module is in a lineage while the run grows:
the share of the module's families the lineage carries, from 0 to 1. The tests use factors of 0 and
1, so each check is exact, and each positive test first requires at least one event.
"""

import pytest

from zombi2.genomes import family, simulate_genomes_family
from zombi2.params import PerCopy, Recipients
from zombi2.params.conditioned import resolve_driver
from zombi2.species import simulate_species_tree

MEMBERS = ("trpA", "trpB", "trpC")


def _only_when_complete(f):
    return 1.0 if f == 1.0 else 0.0


def _only_when_incomplete(f):
    return 0.0 if f == 1.0 else 1.0


def _always(f):
    return 1.0


@pytest.fixture(scope="module")
def tree():
    return simulate_species_tree(birth=1.0, death=0.2, n_extant=40, seed=1).complete_tree


@pytest.fixture(scope="module")
def clade(tree):
    """The lineage whose subtree holds closest to half of the extant tips."""
    extant = set(tree.extant_leaves())
    under = dict.fromkeys(tree.nodes, 0)
    for leaf in extant:
        x = leaf
        while x is not None:
            under[x] += 1
            x = tree.nodes[x].parent
    return min((n for n in tree.nodes if n != tree.root), key=lambda n: abs(under[n] - len(extant) / 2))


def _run(tree, **kw):
    base = dict(duplication=0.1, loss=0.3, transfer=0.5, origination=0.1, initial_families=8, seed=7)
    return simulate_genomes_family(tree, **{**base, **kw})


def _completion(result, tree):
    return resolve_driver(result.completion("trp"), tree, step=None, level="genomes.family")


def _just_before(tree, edge):
    """The instant before an event on its branch, so a loss reads the module as it was."""
    return max(edge.time - 1e-9, tree.nodes[edge.lineage].birth_time)


def _guarded(factor):
    return PerCopy(1.0).scaled_by("genomes:module:trp", factor)


# --- a module's completeness read during the run -------------------------------------------------

@pytest.mark.parametrize("factor, members_lost", [(_only_when_incomplete, False), (_always, True)],
                         ids=["guarded", "unguarded"])
def test_a_complete_module_protects_its_families(tree, factor, members_lost):
    """The module starts complete, and its families can only become absent by being lost, which a
    complete module forbids. Without the guard they are lost."""
    g = _run(tree, joint=True, families=[family(n, module="trp", loss=_guarded(factor)) for n in MEMBERS])
    ids = {g.family_names[n] for n in MEMBERS}
    assert bool([e for e in g.edges if e.kind == "loss" and e.family in ids]) == members_lost
    assert any(e.kind == "loss" and e.family not in ids for e in g.edges)


def test_a_family_is_lost_only_while_its_module_is_incomplete(tree, clade):
    """trpC is planted in one clade, so the module is incomplete everywhere else, and complete inside
    that clade once trpC is there. trpA and trpB are lost only where it is incomplete."""
    rule = _guarded(_only_when_incomplete)
    g = _run(tree, joint=True, families=[
        family("trpA", module="trp", loss=rule), family("trpB", module="trp", loss=rule),
        family("trpC", module="trp", loss=0.0, transfer=0.0, origin=clade)])
    done = _completion(g, tree)
    ids = {g.family_names["trpA"], g.family_names["trpB"]}
    lost = [e for e in g.edges if e.kind == "loss" and e.family in ids]
    assert lost, "no family of the module was lost, so the test shows nothing"
    assert all(done.value(e.lineage, _just_before(tree, e)) < 1.0 for e in lost)


@pytest.mark.parametrize("factor, anything_lost", [(_only_when_incomplete, False), (_always, True)],
                         ids=["guarded", "unguarded"])
def test_a_run_rate_reads_a_module(tree, factor, anything_lost):
    """The run's loss reads the module, so while it is complete no family at all is lost."""
    g = _run(tree, joint=True, loss=_guarded(factor), families=[family(n, module="trp") for n in MEMBERS])
    assert any(e.kind == "loss" for e in g.edges) == anything_lost


def test_a_transfer_to_reads_a_module(tree, clade):
    """X is received only by lineages whose module is complete: inside the clade trpC was planted in."""
    g = _run(tree, joint=True, families=[
        family("trpA", module="trp", loss=0.0), family("trpB", module="trp", loss=0.0),
        family("trpC", module="trp", loss=0.0, transfer=0.0, origin=clade),
        family("X", loss=0.0, transfer_to=Recipients().weighted_by("genomes:module:trp", _only_when_complete))])
    done = _completion(g, tree)
    x = g.family_names["X"]
    arrived = [e for e in g.edges if e.kind == "transfer" and e.recipient is not None and e.family == x]
    assert arrived, "X never arrived, so the test shows nothing"
    assert all(done.value(e.recipient, e.time) == 1.0 for e in arrived)


# --- what is refused -----------------------------------------------------------------------------

def test_a_module_driver_names_a_declared_module(tree):
    with pytest.raises(ValueError, match="does not declare"):
        _run(tree, joint=True, families=[
            family("trpA", module="trp", loss=PerCopy(0.3).scaled_by("genomes:module:his", _only_when_incomplete))])


def test_a_module_driver_takes_a_function_not_a_table(tree):
    with pytest.raises(ValueError, match="CONTINUOUS"):
        _run(tree, joint=True, families=[
            family(n, module="trp", loss=PerCopy(0.3).scaled_by("genomes:module:trp", {"present": 1.0, "absent": 20.0}))
            for n in MEMBERS])


def test_reading_a_module_needs_joint(tree):
    with pytest.raises(ValueError, match="joint=True"):
        _run(tree, families=[family(n, module="trp", loss=_guarded(_only_when_incomplete)) for n in MEMBERS])


def test_a_name_that_is_both_a_module_and_a_family_is_refused(tree):
    with pytest.raises(ValueError, match="declares both"):
        _run(tree, joint=True, families=[
            family("trpA", module="trp"),
            family("module:trp", loss=PerCopy(0.3).scaled_by("genomes:module:trp", _only_when_incomplete))])


def test_a_family_named_like_a_module_is_still_a_family(tree):
    """With no module of that name, "genomes:module:x" reads the family literally named "module:x",
    as it did before a module could be read."""
    g = _run(tree, joint=True, families=[
        family("module:x"), family("B", loss=PerCopy(0.3).scaled_by("genomes:module:x", {"present": 1.0, "absent": 2.0}))])
    assert "module:x" in g.family_names


def test_the_per_family_engine_refuses_a_module_driver(tree):
    with pytest.raises(ValueError, match="per-family engine"):
        _run(tree, joint=True, parallel=True,
             families=[family(n, module="trp", loss=_guarded(_only_when_incomplete)) for n in MEMBERS])
