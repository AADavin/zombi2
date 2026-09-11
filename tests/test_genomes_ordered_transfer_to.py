"""A family's own ``transfer_to``, and a ``transfer_to`` that reads gene content, for ordered genomes
(issue #430).

The rules are the family resolution's (issue #423). What is new is the segment: a transfer moves a
block of neighbouring genes, and the block arrives whole, so a lineage can receive it only where
every rule among its genes allows. The tests use weights of 1 and 0, so each check is exact, and each
positive test first requires at least one transfer. A rule reads a genome as it was just before the
transfer, so presence is read an instant earlier.
"""

import collections

import numpy as np
import pytest

import zombi2.genomes.ordered as ordered
from zombi2.genomes import family, simulate_genomes_ordered
from zombi2.genomes._transfer import recipient_index_all, resolve_transfer_to
from zombi2.genomes.links import Link
from zombi2.params import Between, Fixed, Recipients
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
    return simulate_genomes_ordered(tree, **{**base, **kw})


def _a_in(clade):
    """Family A, planted at ``clade``, never lost and with no transfer rate of its own."""
    return family("A", loss=0.0, transfer=0.0, origin=clade)


def _weighted(driver, table):
    return Recipients().weighted_by(driver, table)


def _presence(result, name, tree):
    return resolve_driver(result.presence(name), tree, step=None, level="genomes.ordered")


def _state(timeline, tree, edge, attr="recipient"):
    """The driver on the edge's ``attr`` lineage an instant before the transfer."""
    node = getattr(edge, attr)
    return timeline.value(node, max(edge.time - 1e-9, tree.nodes[node].birth_time))


def _blocks(result):
    """The arriving copies of each transfer, grouped by (time, donor, recipient): one block each."""
    out = collections.defaultdict(list)
    for e in result.edges:
        if e.kind == "transfer" and e.recipient is not None:
            out[(e.time, e.donor, e.recipient)].append(e)
    return list(out.values())


def _carries(block, result, name):
    return any(e.family == result.family_names[name] for e in block)


# --- one gene at a time: the family resolution's rules ------------------------------------------

@pytest.mark.parametrize("table, allowed", [(ONLY_WITH, "present"), (ONLY_WITHOUT, "absent")],
                         ids=["only-with-A", "only-without-A"])
def test_a_family_rule_reads_another_family(tree, clade, table, allowed):
    g = _run(tree, joint=True, families=[
        _a_in(clade), family("B", loss=0.0, transfer_to=_weighted("genomes:A", table))])
    a = _presence(g, "A", tree)
    arrived = [b[0] for b in _blocks(g) if _carries(b, g, "B")]
    assert arrived, "B never arrived, so the test shows nothing"
    assert all(_state(a, tree, e) == allowed for e in arrived)


def test_the_other_families_keep_the_run_rule(tree, clade):
    g = _run(tree, joint=True, families=[
        _a_in(clade), family("B", loss=0.0, transfer_to=_weighted("genomes:A", ONLY_WITH))])
    a = _presence(g, "A", tree)
    others = [b[0] for b in _blocks(g) if not _carries(b, g, "B")]
    assert any(_state(a, tree, e) == "absent" for e in others)


def test_the_run_rule_reads_a_family(tree, clade):
    g = _run(tree, joint=True, families=[_a_in(clade)], transfer_to=_weighted("genomes:A", ONLY_WITH))
    a = _presence(g, "A", tree)
    arrived = [b[0] for b in _blocks(g)]
    assert arrived
    assert all(_state(a, tree, e) == "present" for e in arrived)


def test_the_run_rule_reads_a_module(tree, clade):
    """x is in one clade and y everywhere, so the module is complete exactly where x is."""
    g = _run(tree, joint=True, families=[
        family("x", loss=0.0, transfer=0.0, origin=clade, module="m"),
        family("y", loss=0.0, transfer=0.0, module="m")],
        transfer_to=Recipients().weighted_by("genomes:module:m", lambda share: 1.0 if share == 1.0 else 0.0))
    x = _presence(g, "x", tree)
    arrived = [b[0] for b in _blocks(g)]
    assert arrived
    assert all(_state(x, tree, e) == "present" for e in arrived)


def test_between_reads_the_donor_and_the_receiver(tree, clade):
    g = _run(tree, joint=True, families=[
        _a_in(clade),
        family("B", loss=0.0, transfer_to=_weighted(
            "genomes:A", Between({("present", "present"): 1.0}, default=0.0)))])
    a = _presence(g, "A", tree)
    arrived = [b[0] for b in _blocks(g) if _carries(b, g, "B")]
    assert arrived
    assert all(_state(a, tree, e, attr="donor") == "present" and _state(a, tree, e) == "present"
               for e in arrived)


def test_a_family_takes_the_other_rules_without_joint(tree):
    """A trait weighting on one family and "distance" on another: no family is read, so no joint run."""
    from zombi2.traits import simulate_discrete

    habitat = simulate_discrete(tree, states=["soil", "water"], switch=0.3, seed=2)
    g = _run(tree, families=[
        family("B", loss=0.0, transfer_to=_weighted(habitat, {"soil": 1.0, "water": 0.0})),
        family("C", loss=0.0, transfer_to="distance")])
    where = resolve_driver(habitat, tree, step=None, level="genomes.ordered")
    blocks = _blocks(g)
    with_b = [b[0] for b in blocks if _carries(b, g, "B")]
    others_in_water = [b[0] for b in blocks if not _carries(b, g, "B")
                       and _state(where, tree, b[0]) == "water"]
    assert with_b and others_in_water and any(_carries(b, g, "C") for b in blocks)
    assert all(_state(where, tree, e) == "soil" for e in with_b)


# --- a segment of several genes ------------------------------------------------------------------

def test_a_segment_carrying_the_family_follows_its_rule(tree, clade):
    """Segments of three genes: every one that carries B lands where A is, and the others often do
    not. A travels as a neighbour too, so the seed is one where many lineages still lack it."""
    g = _run(tree, joint=True, transfer_extent=Fixed(3), seed=1, families=[
        _a_in(clade), family("B", loss=0.0, transfer_to=_weighted("genomes:A", ONLY_WITH))])
    a = _presence(g, "A", tree)
    blocks = _blocks(g)
    with_b = [b for b in blocks if _carries(b, g, "B")]
    assert any(len(b) > 1 for b in with_b), "no segment carried B with a neighbour"
    assert all(_state(a, tree, b[0]) == "present" for b in with_b)
    assert sum(_state(a, tree, b[0]) == "absent" for b in blocks if not _carries(b, g, "B")) >= 20


@pytest.mark.parametrize("d_table, together", [(ONLY_WITH, True), (ONLY_WITHOUT, False)],
                         ids=["rules-agree", "rules-exclude-each-other"])
def test_a_segment_goes_only_where_every_rule_allows(tree, clade, d_table, together):
    """B and D start side by side. B goes only where A is. When D goes there too, segments carrying
    both arrive. When D goes only where A is not, no lineage allows both, and none arrives."""
    g = _run(tree, joint=True, transfer_extent=Fixed(3), initial_families=4, seed=3, families=[
        family("B", loss=0.0, transfer_to=_weighted("genomes:A", ONLY_WITH)),
        family("D", loss=0.0, transfer_to=_weighted("genomes:A", d_table)),
        _a_in(clade)])
    a = _presence(g, "A", tree)
    blocks = _blocks(g)
    both = [b for b in blocks if _carries(b, g, "B") and _carries(b, g, "D")]
    only_d = [b for b in blocks if _carries(b, g, "D") and not _carries(b, g, "B")]
    assert bool(both) == together
    assert only_d
    assert all(_state(a, tree, b[0]) == ("present" if d_table is ONLY_WITH else "absent") for b in only_d)


def test_the_rules_of_a_segment_multiply():
    """Two rules over three lineages: the first allows lineages 0 and 1, the second 1 and 2, so only
    1 is ever picked. With a third rule that allows only 2 beside the first, none is."""
    class States:
        def __init__(self, states):
            self._states = states

        def value(self, node, time):
            return self._states[node]

    rule = resolve_transfer_to(_weighted("genomes:A", {"yes": 1.0, "no": 0.0}))
    first = States({10: "yes", 11: "yes", 12: "no"})
    second = States({10: "no", 11: "yes", 12: "yes"})
    third = States({10: "no", 11: "no", 12: "yes"})
    rng = np.random.default_rng(0)
    alive, cand = [10, 11, 12], [0, 1, 2]
    picks = {recipient_index_all(rng, None, alive, cand, 99, 0.0, [(rule, None, first), (rule, None, second)], 1.0)
             for _ in range(200)}
    assert picks == {1}
    assert recipient_index_all(rng, None, alive, cand, 99, 0.0,
                               [(rule, None, first), (rule, None, third)], 1.0) is None


@pytest.mark.parametrize("run_rule, b_rule, c_rule, several", [
    ("distance", "distance", "distance", False),
    ("shared", "shared", None, False),
    ("uniform", "distance", None, True),
], ids=["equal-names", "same-object", "different-rules"])
def test_a_rule_is_counted_once(tree, monkeypatch, run_rule, b_rule, c_rule, several):
    """A family's rule equal to the run's, or the same object, is one rule, so a segment carrying
    both is picked by that rule alone. Two different rules go through the product."""
    calls = []
    real = ordered.recipient_index_all

    def spy(rng, tree_, alive, cand, donor, t, rules, depth):
        calls.append(len(rules))
        return real(rng, tree_, alive, cand, donor, t, rules, depth)

    monkeypatch.setattr(ordered, "recipient_index_all", spy)
    from zombi2.traits import simulate_discrete

    habitat = simulate_discrete(tree, states=["soil", "water"], switch=0.3, seed=2)
    rules = {"shared": _weighted(habitat, {"soil": 2.0, "water": 1.0})}
    families = [family("B", loss=0.0, transfer_to=rules.get(b_rule, b_rule))]
    if c_rule is not None:
        families.append(family("C", loss=0.0, transfer_to=rules.get(c_rule, c_rule)))
    g = _run(tree, transfer_extent=Fixed(3), families=families, transfer_to=rules.get(run_rule, run_rule))
    b = g.family_names["B"]
    assert any(any(e.family == b for e in block) and any(e.family != b for e in block)
               for block in _blocks(g)), "no segment carried B with another family"
    assert bool(calls) == several
    assert set(calls) <= {2}


# --- what is refused, and what is written ---------------------------------------------------------

def test_reading_gene_content_needs_joint(tree):
    with pytest.raises(ValueError, match="joint=True"):
        _run(tree, families=[family("A"), family("B", transfer_to=_weighted("genomes:A", ONLY_WITH))])


def test_a_rule_reads_only_declared_families(tree):
    with pytest.raises(ValueError, match="does not declare"):
        _run(tree, joint=True, transfer_to=_weighted("genomes:Z", ONLY_WITH))


def test_a_family_rule_is_checked_like_the_run_rule(tree):
    with pytest.raises(ValueError, match="family 'B''s transfer_to"):
        _run(tree, families=[family("B", transfer_to=5.0)])


def test_the_links_of_the_rules(tree):
    g = _run(tree, joint=True, transfer_to=_weighted("genomes:A", {"present": 2.0, "absent": 1.0}),
             families=[family("A"), family("B", transfer_to=_weighted("genomes:A", ONLY_WITH))])
    assert g.links == (
        Link("run", "transfer_to", "genomes:A", "weighted_by", "present=2.0;absent=1.0"),
        Link("B", "transfer_to", "genomes:A", "weighted_by", "present=1.0;absent=0.0"),
    )
