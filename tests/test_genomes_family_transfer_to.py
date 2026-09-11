"""A family's own ``transfer_to``, and a ``transfer_to`` that reads a family's presence (issue #423).

Any recipient rule the run's ``transfer_to`` takes can be given to one declared family, for that
family's copies alone. In a joint run the rule can also be weighted by gene content the run is
building. The tests use weights of 1 and 0, so each check is exact: a copy never arrives where its
weight is 0. Family A is planted in one clade and never lost or transferred, so it is present in
exactly that subtree; B is never lost where a test needs it to keep moving.
"""

import pytest

from zombi2.genomes import family, genome, simulate_genomes_family, simulate_genomes_ordered
from zombi2.params import Between, PerCopy, Recipients
from zombi2.params.conditioned import resolve_driver
from zombi2.species import simulate_species_tree

ONLY_WITH = {"present": 1.0, "absent": 0.0}
ONLY_WITHOUT = {"present": 0.0, "absent": 1.0}


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


def _a_in(clade):
    """Family A, planted at ``clade`` and never lost or transferred."""
    return family("A", loss=0.0, transfer=0.0, origin=clade)


def _weighted(driver, table):
    return Recipients().weighted_by(driver, table)


def _timeline(driver, tree):
    return resolve_driver(driver, tree, step=None, level="genomes.family")


def _arrivals(result, name=None):
    fid = None if name is None else result.family_names[name]
    return [e for e in result.edges if e.kind == "transfer" and e.recipient is not None
            and (fid is None or e.family == fid)]


def _state(timeline, edge, attr="recipient"):
    return timeline.value(getattr(edge, attr), edge.time)


# --- a family's own rule, read off another family ------------------------------------------------

@pytest.mark.parametrize("table, allowed", [(ONLY_WITH, "present"), (ONLY_WITHOUT, "absent")],
                         ids=["only-with-A", "only-without-A"])
def test_a_family_rule_reads_another_family(tree, clade, table, allowed):
    g = _run(tree, joint=True, families=[
        _a_in(clade), family("B", loss=0.0, transfer_to=_weighted("genomes:A", table))])
    a = _timeline(g.presence("A"), tree)
    arrived = _arrivals(g, "B")
    assert arrived, "B never arrived, so the test shows nothing"
    assert all(_state(a, e) == allowed for e in arrived)


def test_the_other_families_keep_the_run_rule(tree, clade):
    g = _run(tree, joint=True, families=[
        _a_in(clade), family("B", loss=0.0, transfer_to=_weighted("genomes:A", ONLY_WITH))])
    a = _timeline(g.presence("A"), tree)
    others = [e for e in _arrivals(g) if e.family not in (g.family_names["A"], g.family_names["B"])]
    assert any(_state(a, e) == "absent" for e in others)


def test_the_run_rule_reads_a_family(tree, clade):
    """The run's own transfer_to weighted by live presence is a joint run, with no rate reading it."""
    g = _run(tree, joint=True, families=[_a_in(clade)], transfer_to=_weighted("genomes:A", ONLY_WITH))
    a = _timeline(g.presence("A"), tree)
    arrived = _arrivals(g)
    assert arrived
    assert all(_state(a, e) == "present" for e in arrived)


def test_two_families_each_arrive_where_the_other_is(tree, clade):
    g = _run(tree, joint=True, families=[
        family("A", loss=0.0, transfer=PerCopy(0.5), origin=clade,
               transfer_to=_weighted("genomes:B", ONLY_WITH)),
        family("B", loss=0.0, transfer_to=_weighted("genomes:A", ONLY_WITH))])
    a, b = _timeline(g.presence("A"), tree), _timeline(g.presence("B"), tree)
    a_arrived, b_arrived = _arrivals(g, "A"), _arrivals(g, "B")
    assert a_arrived and b_arrived
    assert all(_state(b, e) == "present" for e in a_arrived)
    assert all(_state(a, e) == "present" for e in b_arrived)


def test_gain_and_loss_read_the_same_family(tree, clade):
    """B arrives only where A is and is lost only where A is not. B is not lost before A exists."""
    start = tree.nodes[clade].birth_time
    g = _run(tree, joint=True, families=[
        _a_in(clade),
        family("B", transfer_to=_weighted("genomes:A", ONLY_WITH),
               loss=PerCopy(0.5).changing_at({0: 0.0, start: 1.0}).scaled_by("genomes:A", ONLY_WITHOUT))])
    a = _timeline(g.presence("A"), tree)
    b = g.family_names["B"]
    arrived = _arrivals(g, "B")
    lost = [e for e in g.edges if e.kind == "loss" and e.family == b]
    assert arrived and lost
    assert all(_state(a, e) == "present" for e in arrived)
    assert all(_state(a, e, attr="lineage") == "absent" for e in lost)


def test_between_reads_the_donor_and_the_receiver(tree, clade):
    g = _run(tree, joint=True, families=[
        _a_in(clade),
        family("B", loss=0.0, transfer_to=_weighted(
            "genomes:A", Between({("present", "present"): 1.0}, default=0.0)))])
    a = _timeline(g.presence("A"), tree)
    arrived = _arrivals(g, "B")
    assert arrived
    assert all(_state(a, e, attr="donor") == "present" and _state(a, e) == "present" for e in arrived)


def test_when_every_candidate_weighs_zero_nothing_arrives(tree):
    g = _run(tree, joint=True, families=[
        family("A"), family("B", loss=0.0, transfer_to=_weighted("genomes:A", {"present": 0.0, "absent": 0.0}))])
    assert not _arrivals(g, "B")
    assert _arrivals(g), "the other families still transfer"


def test_a_family_takes_the_other_rules_without_joint(tree):
    """A trait weighting on one family and "distance" on another: no family is read, so no joint run."""
    from zombi2.traits import simulate_discrete

    habitat = simulate_discrete(tree, states=["soil", "water"], switch=0.3, seed=2)
    g = _run(tree, families=[
        family("B", loss=0.0, transfer_to=_weighted(habitat, {"soil": 1.0, "water": 0.0})),
        family("C", loss=0.0, transfer_to="distance")])
    where = _timeline(habitat, tree)
    b = g.family_names["B"]
    arrived = _arrivals(g, "B")
    others_in_water = [e for e in _arrivals(g) if e.family != b and _state(where, e) == "water"]
    assert arrived and others_in_water and _arrivals(g, "C")
    assert all(_state(where, e) == "soil" for e in arrived)


# --- what is refused -----------------------------------------------------------------------------

def test_reading_a_family_needs_joint(tree):
    with pytest.raises(ValueError, match="joint=True"):
        _run(tree, families=[family("A"), family("B", transfer_to=_weighted("genomes:A", ONLY_WITH))])


def test_a_rule_reads_only_declared_families(tree):
    with pytest.raises(ValueError, match="does not declare"):
        _run(tree, joint=True, families=[family("B", transfer_to=_weighted("genomes:Z", ONLY_WITH))])


def test_a_rule_table_names_real_states(tree):
    with pytest.raises(ValueError, match="not among the driver's states"):
        _run(tree, joint=True, families=[
            family("A"), family("B", transfer_to=_weighted("genomes:A", {"presnt": 1.0, "absent": 0.0}))])


def test_a_family_rule_is_checked_like_the_run_rule(tree):
    with pytest.raises(ValueError, match="family 'B''s transfer_to"):
        _run(tree, families=[family("B", transfer_to=PerCopy(0.1))])


@pytest.mark.parametrize("kw", [
    dict(joint=True, families=[family("A")], transfer_to=Recipients().weighted_by("genomes:A", ONLY_WITH)),
    dict(families=[family("B", transfer_to="distance")]),
], ids=["run-rule-reads-a-family", "family-rule"])
def test_the_per_family_engine_refuses(tree, kw):
    with pytest.raises(ValueError, match="per-family engine"):
        _run(tree, parallel=True, **kw)


def test_the_ordered_engine_refuses_a_family_rule(tree):
    with pytest.raises(ValueError, match="own transfer_to"):
        simulate_genomes_ordered(tree, families=[family("B", transfer_to="distance")], seed=1)


def test_a_joint_genome_spec_refuses_a_family_rule():
    with pytest.raises(ValueError, match="transfer_to"):
        genome(families=[family("B", transfer_to="distance")])
