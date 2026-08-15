"""The genome level joined to itself — a rate reading gene content this same run is producing.

The design note's §6 loop. The driver is a live **name** (``"genomes:IS1"``, ``"genomes:count"``)
rather than a finished result, so neither the driver nor what it drives can be simulated first, and
``joint=True`` is how a run says so. These cover the signal, the exactness, and the refusals.
"""

import collections

import pytest

from zombi2.genomes import family, simulate_genomes_family
from zombi2.params import Curve, PerCopy
from zombi2.species import simulate_species_tree


def _tree(n=25, seed=4):
    return simulate_species_tree(birth=1.0, n_extant=n, seed=seed).complete_tree


def _donations(result, exclude=None):
    """Transfers **given** by each lineage. A driven ``transfer`` weights the donor, and every
    transfer writes two edges, so the recipient-side one is counted and grouped by its donor."""
    return collections.Counter(e.donor for e in result.edges
                               if e.kind == "transfer" and e.recipient is not None
                               and e.family != exclude)


# --- the signal -----------------------------------------------------------------------------------

def _mobile(seed=7, factor=30.0, **kw):
    tree = _tree(30)
    g = simulate_genomes_family(
        tree, initial_families=25, duplication=0.05, loss=0.12, seed=seed, joint=True,
        max_family_size=8,
        families=[family("IS1", origin=(11, None), transfer=PerCopy(0.30), loss=0.08)],
        transfer=PerCopy(0.025).scaled_by("genomes:IS1", {"present": factor, "absent": 1.0}), **kw)
    return tree, g


def test_carrying_the_element_makes_a_lineage_donate():
    tree, g = _mobile()
    is1 = g.family_names["IS1"]
    given = _donations(g, exclude=is1)
    tips = list(tree.extant_leaves())
    carriers = [n for n in tips if g.family_counts(n)[is1] > 0]
    others = [n for n in tips if n not in carriers]
    assert carriers and others, "the run needs both kinds of tip for the comparison to mean anything"
    per = lambda ns: sum(given.get(n, 0) for n in ns) / len(ns)
    assert per(carriers) > 5 * per(others), f"{per(carriers):.2f} against {per(others):.2f}"


def test_a_flat_mapping_is_the_undriven_model():
    """The same run with the factor at 1.0 must lose the signal. Otherwise the difference above could
    be about which lineages carry the element rather than about the driving."""
    tree, flat = _mobile(factor=1.0)
    is1 = flat.family_names["IS1"]
    given = _donations(flat, exclude=is1)
    tips = list(tree.extant_leaves())
    carriers = [n for n in tips if flat.family_counts(n)[is1] > 0]
    others = [n for n in tips if n not in carriers]
    per = lambda ns: sum(given.get(n, 0) for n in ns) / max(len(ns), 1)
    assert per(carriers) < 2 * max(per(others), 0.05)


def test_the_element_reaches_lineages_that_never_inherited_it():
    """The loop, not just the driving: the element moves itself, so it turns up in clades outside the
    one it was planted on, and those become donors in their turn."""
    tree, g = _mobile()
    is1 = g.family_names["IS1"]
    planted = 11
    below = set()
    stack = [planted]
    while stack:
        i = stack.pop()
        below.add(i)
        if tree.nodes[i].children:
            stack.extend(tree.nodes[i].children)
    carriers = {n for n in tree.extant_leaves() if g.family_counts(n)[is1] > 0}
    assert carriers - below, "the element never left the clade it was planted on"


def test_a_count_driver_reads_the_whole_genome():
    tree = _tree(20)
    g = simulate_genomes_family(
        tree, initial_families=10, loss=0.05, seed=2, joint=True,
        duplication=PerCopy(0.05).scaled_by("genomes:count", Curve(lambda n: 1.0 + n / 40.0)))
    assert len(g.events) > 10


def test_it_is_deterministic():
    a, b = _mobile()[1], _mobile()[1]
    assert [(e.time, e.kind, e.family, e.copy) for e in a.edges] == \
           [(e.time, e.kind, e.family, e.copy) for e in b.edges]


# --- what is refused ------------------------------------------------------------------------------

def test_a_live_driver_needs_joint_true():
    tree = _tree(10)
    with pytest.raises(ValueError, match="the run is joint"):
        simulate_genomes_family(tree, initial_families=5, loss=0.1, seed=1,
                                families=[family("IS1")],
                                transfer=PerCopy(0.05).scaled_by("genomes:IS1",
                                                                 {"present": 3.0, "absent": 1.0}))


def test_joint_true_needs_something_to_drive():
    tree = _tree(10)
    with pytest.raises(ValueError, match="no rate reads live gene content"):
        simulate_genomes_family(tree, initial_families=5, loss=0.1, transfer=0.05, seed=1,
                                joint=True, families=[family("IS1")])


def test_a_live_driver_must_name_a_declared_family():
    tree = _tree(10)
    with pytest.raises(ValueError, match="does not declare"):
        simulate_genomes_family(tree, initial_families=5, loss=0.1, seed=1, joint=True,
                                families=[family("IS1")],
                                transfer=PerCopy(0.05).scaled_by("genomes:IS7",
                                                                 {"present": 3.0, "absent": 1.0}))


def test_a_mapping_that_can_never_fire_is_refused():
    """The alphabet is known before the run starts — a family is present or absent — so the check is
    exhaustive. A typo would leave every lineage at the default factor: a run that reports it was
    driven and was not."""
    tree = _tree(10)
    with pytest.raises(ValueError, match="not among the driver's states"):
        simulate_genomes_family(tree, initial_families=5, loss=0.1, seed=1, joint=True,
                                families=[family("IS1")],
                                transfer=PerCopy(0.05).scaled_by("genomes:IS1",
                                                                 {"presnt": 3.0, "absent": 1.0}))


def test_a_table_on_the_count_driver_is_refused():
    tree = _tree(10)
    with pytest.raises(ValueError, match="CONTINUOUS"):
        simulate_genomes_family(tree, initial_families=5, loss=0.1, seed=1, joint=True,
                                duplication=PerCopy(0.05).scaled_by("genomes:count", {"many": 3.0}))


def test_the_per_family_engine_refuses_a_joint_run():
    tree = _tree(10)
    with pytest.raises(ValueError, match="per-family engine"):
        simulate_genomes_family(tree, initial_families=5, loss=0.1, seed=1, joint=True, parallel=True,
                                families=[family("IS1")],
                                transfer=PerCopy(0.05).scaled_by("genomes:IS1",
                                                                 {"present": 3.0, "absent": 1.0}))


def test_a_runaway_joint_run_raises():
    """A rate reading the genome's own content can feed itself. It raises rather than stopping early,
    for the reason the species engine gives: a run cut off at a size is no longer a sample from the
    process asked for."""
    tree = _tree(12)
    with pytest.raises(RuntimeError, match="feeding itself"):
        simulate_genomes_family(
            tree, initial_families=20, loss=0.0, seed=1, joint=True, max_family_size=None,
            duplication=PerCopy(0.5).scaled_by("genomes:count", Curve(lambda n: 1.0 + n / 5.0)))
