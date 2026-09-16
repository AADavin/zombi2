"""``lineage_multipliers.tsv``: each species branch's drawn rate multipliers.

A rate written with ``varying_among("lineages", law)`` draws a number for each branch of the species
tree before the run starts, and every family passing through that branch is multiplied by it. The run
keeps the numbers on the result, writes them with its other files, and ``read_run`` reads them back.
A run whose rates do not vary among lineages writes the header alone.
"""

import collections

import pytest

from zombi2.genomes import read_run, simulate_genomes_family, simulate_genomes_ordered
from zombi2.genomes.multipliers import (LINEAGE_TARGETS, ORDERED_LINEAGE_TARGETS,
                                        lineage_multipliers_from_tsv, lineage_multipliers_header)
from zombi2.params import Drift, LogNormal, PerChromosome, PerCopy, PerLineage, Random
from zombi2.species import simulate_species_tree


@pytest.fixture(scope="module")
def tree():
    return simulate_species_tree(birth=1.0, death=0.3, n_extant=15, seed=1).complete_tree


def _header():
    return lineage_multipliers_header() + "\n"


def _duplications_by_branch(run):
    """How many duplications happened on each species branch, from the gene trees — the event log
    records the family and the genes, and the branch is what the gene tree knows."""
    got: collections.Counter = collections.Counter()
    for gt in run.gene_trees.values():
        stack = [gt.complete]
        while stack:
            node = stack.pop()
            if node.kind == "duplication":
                got[node.species] += 1
            stack.extend(node.children)
    return got


def test_a_run_whose_rates_do_not_vary_writes_the_header_alone(tree, tmp_path):
    run = simulate_genomes_family(tree, duplication=0.2, transfer=0.3, loss=0.25,
                                  initial_families=5, seed=2)
    assert run.lineage_multipliers == {}
    run.write(tmp_path)
    assert (tmp_path / "lineage_multipliers.tsv").read_text(encoding="utf-8") == _header()


def test_every_branch_has_a_row_extinct_ones_included(tree):
    run = simulate_genomes_family(
        tree, duplication=PerCopy(0.2).varying_among("lineages", LogNormal(0.0, 0.5)),
        loss=0.2, initial_families=5, seed=2)
    assert set(run.lineage_multipliers) == set(tree.nodes)
    extinct = [i for i, nd in tree.nodes.items() if nd.fate == "extinct"]
    assert extinct and all(i in run.lineage_multipliers for i in extinct)
    assert all(set(row) == set(LINEAGE_TARGETS) for row in run.lineage_multipliers.values())


def test_a_rate_that_does_not_vary_gives_every_branch_one(tree):
    run = simulate_genomes_family(
        tree, duplication=PerCopy(0.2).varying_among("lineages", LogNormal(0.0, 0.5)),
        loss=0.2, initial_families=5, seed=2)
    assert {row["loss"] for row in run.lineage_multipliers.values()} == {1.0}
    assert len({row["duplication"] for row in run.lineage_multipliers.values()}) > 1


def test_one_random_on_two_rates_is_one_draw_per_branch(tree):
    """Writing one object or two is what decides whether a branch is fast at both (SPEC §5)."""
    speed = Random("lineages", LogNormal(0.0, 0.5))
    shared = simulate_genomes_family(tree, duplication=PerCopy(0.2).varying_among(speed),
                                     loss=PerCopy(0.2).varying_among(speed),
                                     initial_families=5, seed=2)
    assert all(row["duplication"] == row["loss"] for row in shared.lineage_multipliers.values())

    apart = simulate_genomes_family(
        tree, duplication=PerCopy(0.2).varying_among("lineages", LogNormal(0.0, 0.5)),
        loss=PerCopy(0.2).varying_among("lineages", LogNormal(0.0, 0.5)),
        initial_families=5, seed=2)
    assert any(row["duplication"] != row["loss"] for row in apart.lineage_multipliers.values())


def test_an_inherited_draw_keeps_a_daughter_near_its_parent(tree):
    """``Drift`` is the autocorrelated form: a daughter starts from its parent's value and is
    nudged. The independent draw is the null — the same spread, none of it heritable."""
    def gaps(rate):
        run = simulate_genomes_family(tree, duplication=rate, loss=0.2, initial_families=5, seed=2)
        table = run.lineage_multipliers
        return [abs(table[c]["duplication"] - table[i]["duplication"])
                for i, node in tree.nodes.items() for c in node.children]

    drifted = gaps(PerCopy(0.2).varying_among("lineages", Drift(LogNormal(0.0, 0.5))))
    drawn = gaps(PerCopy(0.2).varying_among("lineages", LogNormal(0.0, 0.5)))
    assert sum(drifted) / len(drifted) < sum(drawn) / len(drawn)


def test_a_branch_with_a_larger_multiplier_gets_more_duplications():
    """The number has to reach the events. duplication is PerLineage, so a branch expects
    ``base × multiplier × length`` duplications whatever its genome holds, and the cap is off so
    none is discarded."""
    tree = simulate_species_tree(birth=1.0, death=0.0, n_extant=16, seed=5).complete_tree
    length = {i: nd.end_time - nd.birth_time for i, nd in tree.nodes.items()}
    base, reps = 0.5, 60
    got: collections.Counter = collections.Counter()
    want: dict = collections.defaultdict(float)
    for rep in range(reps):
        run = simulate_genomes_family(
            tree, duplication=PerLineage(base).varying_among("lineages", LogNormal(0.0, 0.8)),
            loss=0.0, transfer=0.0, origination=0.0, initial_families=2,
            max_family_size=None, seed=1000 + rep)
        got += _duplications_by_branch(run)
        for i, row in run.lineage_multipliers.items():
            want[i] += base * row["duplication"] * length[i]
    rows = sorted(((want[i], got[i]) for i in tree.nodes if want[i] > 20), key=lambda p: p[0])
    assert len(rows) >= 8
    # both halves, so a run that ignored the multipliers (every branch at the base rate) fails:
    # the slow half would overshoot and the fast half fall short
    half = len(rows) // 2
    for part in (rows[:half], rows[half:]):
        expected = sum(w for w, _ in part)
        counted = sum(g for _, g in part)
        assert abs(counted - expected) < 4.0 * expected ** 0.5, (expected, counted)


def test_the_parallel_engine_gives_the_same_table_for_any_worker_count(tree):
    kw = dict(duplication=PerCopy(0.2).varying_among("lineages", LogNormal(0.0, 0.5)),
              origination=PerLineage(0.3).varying_among("lineages", Drift(LogNormal(0.0, 0.3))),
              loss=0.2, initial_families=5, seed=3)
    one = simulate_genomes_family(tree, parallel=1, **kw)
    four = simulate_genomes_family(tree, parallel=4, **kw)
    assert one.lineage_multipliers == four.lineage_multipliers
    assert [(e.time, e.kind) for e in one.events] == [(e.time, e.kind) for e in four.events]


def test_a_streamed_family_run_writes_the_table_the_run_in_memory_writes(tree, tmp_path):
    kw = dict(duplication=PerCopy(0.2).varying_among("lineages", LogNormal(0.0, 0.5)),
              loss=0.2, initial_families=5, seed=3)
    simulate_genomes_family(tree, stream_to=tmp_path / "streamed", **kw)
    in_memory = simulate_genomes_family(tree, parallel=1, **kw)
    in_memory.write(tmp_path / "memory")
    assert ((tmp_path / "streamed" / "lineage_multipliers.tsv").read_text(encoding="utf-8")
            == (tmp_path / "memory" / "lineage_multipliers.tsv").read_text(encoding="utf-8"))


def test_read_run_gives_the_multipliers_back(tree, tmp_path):
    run = simulate_genomes_family(
        tree, duplication=PerCopy(0.2).varying_among("lineages", LogNormal(0.0, 0.5)),
        loss=0.2, initial_families=5, seed=2)
    run.write(tmp_path)
    assert read_run(tmp_path).lineage_multipliers == run.lineage_multipliers


def test_the_file_names_branches_the_way_the_other_files_do(tree, tmp_path):
    run = simulate_genomes_family(
        tree, duplication=PerCopy(0.2).varying_among("lineages", LogNormal(0.0, 0.5)),
        loss=0.2, initial_families=5, seed=2)
    run.write(tmp_path)
    text = (tmp_path / "lineage_multipliers.tsv").read_text(encoding="utf-8")
    labels = set(tree.labels().values())
    named = {line.split("\t")[0] for line in text.splitlines()[1:]}
    assert named == labels
    assert lineage_multipliers_from_tsv(text, tree.labels()) == run.lineage_multipliers


def test_lineage_multipliers_is_an_output_that_can_be_left_out(tree, tmp_path):
    run = simulate_genomes_family(
        tree, duplication=PerCopy(0.2).varying_among("lineages", LogNormal(0.0, 0.5)),
        loss=0.2, initial_families=5, seed=2)
    run.write(tmp_path, outputs=("events",))
    assert not (tmp_path / "lineage_multipliers.tsv").exists()


def test_a_table_with_another_header_is_refused():
    with pytest.raises(ValueError, match="lineage_multipliers.tsv must start with"):
        lineage_multipliers_from_tsv("family\tduplication\n0\t1.0\n")


def test_the_ordered_engine_draws_for_every_rate_it_has(tree):
    """The ordered resolution counts rearrangements and chromosome events too, and a multiplier per
    branch scales whatever the event acts on, so every one of its rates takes the draw."""
    run = simulate_genomes_ordered(
        tree,
        duplication=PerCopy(0.2).varying_among("lineages", LogNormal(0.0, 0.5)),
        inversion=PerCopy(0.1).varying_among("lineages", LogNormal(0.0, 0.5)),
        fission=PerChromosome(0.05).varying_among("lineages", LogNormal(0.0, 0.5)),
        chromosomes=2, loss=0.1, initial_families=8, seed=2)
    assert set(run.lineage_multipliers) == set(tree.nodes)
    assert set(next(iter(run.lineage_multipliers.values()))) == set(ORDERED_LINEAGE_TARGETS)
    for target in ("duplication", "inversion", "fission"):
        assert len({row[target] for row in run.lineage_multipliers.values()}) > 1
    # a rate carrying no draw reads 1.0 on every branch, as at the family resolution
    assert {row["loss"] for row in run.lineage_multipliers.values()} == {1.0}


def test_an_ordered_run_writes_the_wider_table(tree, tmp_path):
    run = simulate_genomes_ordered(
        tree, duplication=PerCopy(0.2).varying_among("lineages", LogNormal(0.0, 0.5)),
        loss=0.1, initial_families=8, seed=2)
    run.write(tmp_path)
    text = (tmp_path / "lineage_multipliers.tsv").read_text(encoding="utf-8")
    assert text.splitlines()[0] == lineage_multipliers_header(ORDERED_LINEAGE_TARGETS)
    assert lineage_multipliers_from_tsv(text, tree.labels()) == run.lineage_multipliers


def test_an_ordered_run_that_does_not_vary_writes_the_header_alone(tree, tmp_path):
    run = simulate_genomes_ordered(tree, duplication=0.2, loss=0.1, initial_families=8, seed=2)
    assert run.lineage_multipliers == {}
    run.write(tmp_path)
    assert ((tmp_path / "lineage_multipliers.tsv").read_text(encoding="utf-8")
            == lineage_multipliers_header(ORDERED_LINEAGE_TARGETS) + "\n")


def test_a_flat_draw_beside_a_family_draw_changes_nothing(tree):
    """The two draws compose, and the ordered engine keeps the product in one weight per lineage.
    A per-lineage draw with no spread multiplies every branch by exactly 1.0 and consumes no
    randomness, so the run through the combined path has to be the run without it, to the last bit —
    which is what pins that weight against the one the per-family draw alone builds."""
    kw = dict(loss=0.15, transfer=0.05, inversion=0.1, initial_families=8, seed=11,
              max_family_size=None)
    fam = PerCopy(0.3).varying_among("families", LogNormal(0.0, 0.6))
    flat = (PerCopy(0.3).varying_among("families", LogNormal(0.0, 0.6))
            .varying_among("lineages", LogNormal(0.0, 0.0)))
    a = simulate_genomes_ordered(tree, duplication=fam, **kw)
    b = simulate_genomes_ordered(tree, duplication=flat, **kw)
    assert [(e.time, e.kind, e.family, e.parents, e.children) for e in a.events] == \
           [(e.time, e.kind, e.family, e.parents, e.children) for e in b.events]
    assert {v for row in b.lineage_multipliers.values() for v in row.values()} == {1.0}


def test_an_ordered_run_can_carry_a_draw_on_every_axis_at_once(tree):
    """A driver, a per-family draw and a per-lineage draw: the first two are refused together, the
    other pairs run, and the engine's own row checks hold for each."""
    from zombi2.params import Clade

    clade = Clade({"fast": ["n3"]})
    both_draws = simulate_genomes_ordered(
        tree, duplication=(PerCopy(0.2).varying_among("lineages", LogNormal(0.0, 0.5))
                           .varying_among("families", LogNormal(0.0, 0.5))),
        loss=0.1, initial_families=8, seed=2)
    assert both_draws.lineage_multipliers and both_draws.family_multipliers

    driven = simulate_genomes_ordered(
        tree, loss=(PerCopy(0.15).scaled_by(clade, {"fast": 3.0, "rest": 1.0})
                    .varying_among("lineages", LogNormal(0.0, 0.5))),
        duplication=0.1, initial_families=8, seed=2)
    assert driven.lineage_multipliers


def test_a_rate_cannot_be_both_drawn_and_inherited_among_lineages(tree):
    speed = Random("lineages", LogNormal(0.0, 0.5))
    drift = Random("lineages", Drift(LogNormal(0.0, 0.5)))
    with pytest.raises(ValueError, match="drawn and an inherited value among lineages"):
        simulate_genomes_family(tree, duplication=PerCopy(0.2).varying_among(speed).varying_among(drift),
                                loss=0.2, initial_families=5, seed=2)


def test_a_table_read_beside_another_tree_is_refused(tree):
    """The table and the tree are written together, so a branch the tree does not hold means the two
    are from different runs — said, rather than kept under a key that looks like a node id."""
    text = lineage_multipliers_header() + "\nnot_a_branch\t1.0\t1.0\t1.0\t1.0\n"
    with pytest.raises(ValueError, match="not in the species tree read with it"):
        lineage_multipliers_from_tsv(text, tree.labels())
    # without a tree the same table reads, keyed by the name it carries
    assert lineage_multipliers_from_tsv(text) == {
        "not_a_branch": dict.fromkeys(LINEAGE_TARGETS, 1.0)}


def test_a_family_that_writes_its_own_rate_is_not_scaled_by_the_branch():
    """A written rate replaces the run's rate, so the branch multiplier, which was written on the
    run's rate, has nothing there to scale. The same reading the per-family draw gets, and the
    weights are where it is decided, so they are what the test reads."""
    from zombi2.genomes.family import _family_weights

    # one lineage: the families running at the run's rate weigh 2.0 between them, and a family that
    # wrote its own rate contributes 5.0. A branch multiplier of 3.0 scales the first, not the second
    assert _family_weights(1.0, [2.0], [5.0], [3.0]) == [1.0 * 3.0 * 2.0 + 5.0]
    # and with no draw on the rate it is the expression it always was
    assert _family_weights(1.0, [2.0], [5.0]) == [1.0 * 2.0 + 5.0]


def test_a_family_that_writes_its_own_rate_still_lets_the_run_vary(tree):
    """The run's own families follow the branch multipliers while the declared one does not, so both
    kinds can sit in one run."""
    from zombi2.genomes import family as declare

    run = simulate_genomes_family(
        tree, duplication=PerCopy(0.4).varying_among("lineages", LogNormal(0.0, 0.9)),
        loss=0.2, initial_families=5, families=[declare("steady", duplication=PerCopy(0.4))],
        seed=2)
    assert len({row["duplication"] for row in run.lineage_multipliers.values()}) > 1
    assert "steady" in run.family_names


def test_a_family_rate_refuses_the_draw_and_says_where_it_goes(tree):
    from zombi2.genomes import family as declare

    with pytest.raises(ValueError, match="belongs on the run's"):
        simulate_genomes_family(
            tree, loss=0.2, initial_families=5,
            families=[declare("fast",
                              loss=PerCopy(0.4).varying_among("lineages", LogNormal(0.0, 0.5)))],
            seed=2)


def test_an_ordered_branch_with_a_larger_multiplier_gets_more_inversions():
    """The same check the family resolution gets, on a rearrangement: inversion is PerLineage, so a
    branch expects ``base × multiplier × length`` inversions whatever its genome holds."""
    tree = simulate_species_tree(birth=1.0, death=0.0, n_extant=16, seed=5).complete_tree
    length = {i: nd.end_time - nd.birth_time for i, nd in tree.nodes.items()}
    base, reps = 0.5, 60
    got: collections.Counter = collections.Counter()
    want: dict = collections.defaultdict(float)
    for rep in range(reps):
        run = simulate_genomes_ordered(
            tree, duplication=0.0, loss=0.0, transfer=0.0, origination=0.0,
            inversion=PerLineage(base).varying_among("lineages", LogNormal(0.0, 0.8)),
            initial_families=2, max_family_size=None, seed=1000 + rep)
        got += collections.Counter(x.lineage for x in run.rearrangements
                                   if type(x).__name__ == "Inversion")
        for i, row in run.lineage_multipliers.items():
            want[i] += base * row["inversion"] * length[i]
    rows = sorted(((want[i], got[i]) for i in tree.nodes if want[i] > 20), key=lambda p: p[0])
    assert len(rows) >= 8
    half = len(rows) // 2
    for part in (rows[:half], rows[half:]):
        expected = sum(w for w, _ in part)
        counted = sum(g for _, g in part)
        assert abs(counted - expected) < 4.0 * expected ** 0.5, (expected, counted)


def test_a_streamed_ordered_run_writes_the_table(tree, tmp_path):
    kw = dict(duplication=PerCopy(0.2).varying_among("lineages", LogNormal(0.0, 0.5)),
              loss=0.1, initial_families=8, seed=3)
    simulate_genomes_ordered(tree, stream_to=tmp_path / "streamed", **kw)
    in_memory = simulate_genomes_ordered(tree, **kw)
    in_memory.write(tmp_path / "memory")
    assert ((tmp_path / "streamed" / "lineage_multipliers.tsv").read_text(encoding="utf-8")
            == (tmp_path / "memory" / "lineage_multipliers.tsv").read_text(encoding="utf-8"))
