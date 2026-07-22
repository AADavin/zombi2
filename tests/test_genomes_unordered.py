"""Tests for the unordered D/L/O gene-family core (zombi2.genomes_unordered)."""

import collections
import pathlib
import tempfile

import pytest

from zombi2.rates import modifiers as mod
from zombi2.species import simulate_species_tree
from zombi2.genomes import simulate_genomes_unordered


def _tree(seed=1, n_extant=12, death=0.3):
    return simulate_species_tree(birth=1.0, death=death, n_extant=n_extant, seed=seed)


def _transfers(events):
    """Reconstruct each transfer as ``(donor_lineage, recipient_lineage, time)`` from its two rows
    (same time + parent: the continuation on the donor, and the copy on the recipient — the one that
    names a recipient)."""
    rows = collections.defaultdict(list)
    for e in events:
        if e.kind == "transfer":
            rows[(e.time, e.parent)].append(e)
    out = []
    for pair in rows.values():
        xfer = next(r for r in pair if r.recipient is not None)     # the transferred copy, on the recipient
        cont = next(r for r in pair if r.recipient is None)         # the donor's continuation
        out.append((cont.lineage, xfer.lineage, xfer.time))
    return out


# --- the walk covers the whole complete tree -------------------------------

def test_genomes_on_every_node_including_extinct():
    sp = _tree(seed=3, death=0.5)
    g = simulate_genomes_unordered(sp, duplication=0.2, loss=0.2, origination=0.5, initial_families=4, seed=1)
    assert set(g.genomes) == set(sp.complete_tree.nodes)          # every node has a genome
    extinct = {n.id for n in sp.complete_tree.extinct()}
    assert extinct and extinct <= set(g.genomes)                 # extinct lineages included


def test_accepts_a_result_or_a_bare_tree():
    sp = _tree(seed=7)
    a = simulate_genomes_unordered(sp, origination=0.5, initial_families=3, seed=1)
    b = simulate_genomes_unordered(sp.complete_tree, origination=0.5, initial_families=3, seed=1)
    assert [(e.time, e.kind, e.copy) for e in a.events] == [(e.time, e.kind, e.copy) for e in b.events]


# --- determinism -----------------------------------------------------------

def test_deterministic_given_seed():
    sp = _tree(seed=2)
    kw = dict(duplication=0.3, loss=0.2, origination=0.5, initial_families=4, seed=9)
    a, b = simulate_genomes_unordered(sp, **kw), simulate_genomes_unordered(sp, **kw)
    assert [(e.time, e.kind, e.lineage, e.copy, e.parent) for e in a.events] == \
           [(e.time, e.kind, e.lineage, e.copy, e.parent) for e in b.events]
    assert a.genomes == b.genomes


def test_different_seeds_differ():
    sp = _tree(seed=2)
    a = simulate_genomes_unordered(sp, duplication=0.3, loss=0.2, origination=0.5, initial_families=4, seed=1)
    b = simulate_genomes_unordered(sp, duplication=0.3, loss=0.2, origination=0.5, initial_families=4, seed=2)
    assert len(a.events) != len(b.events) or \
        [e.copy for e in a.events] != [e.copy for e in b.events]


# --- the three events behave --------------------------------------------------

def test_initial_families_seed_originations_at_the_crown():
    sp = _tree(seed=1)
    root = sp.complete_tree.root
    t0 = sp.complete_tree.nodes[root].birth_time
    g = simulate_genomes_unordered(sp, initial_families=6, seed=1)          # no D/L/O rates
    crown = [e for e in g.events if e.time == t0]
    assert len(crown) == 6 and all(e.kind == "origination" for e in crown)
    assert len(g.genomes[root]) == 6                               # the root carries all 6


def test_no_rates_means_pure_inheritance():
    # with every rate 0, nothing happens beyond the crown seeding + the speciation re-ids: every node
    # carries the root's families (per-node ids differ, but the family multiset is unchanged)
    sp = _tree(seed=4)
    g = simulate_genomes_unordered(sp, initial_families=3, seed=1)
    root_counts = g.family_counts(sp.complete_tree.root)
    assert all(g.family_counts(i) == root_counts for i in g.genomes)   # families inherited unchanged
    assert {e.kind for e in g.events} <= {"origination", "speciation"}  # only crown births + splits
    assert all(e.time == 0.0 for e in g.events if e.kind == "origination")


def test_origination_only_families_never_exceed_one_copy():
    sp = _tree(seed=6)
    g = simulate_genomes_unordered(sp, origination=0.8, initial_families=2, seed=1)  # no duplication
    for node_id in g.genomes:
        assert all(count == 1 for count in g.family_counts(node_id).values())


def test_duplication_grows_a_family():
    sp = _tree(seed=6)
    g = simulate_genomes_unordered(sp, duplication=0.8, initial_families=3, seed=1)  # no loss, no origination
    biggest = max((max(g.family_counts(i).values(), default=0)) for i in g.genomes)
    assert biggest > 1                                            # some family reached >1 copy
    # duplication never introduces a new family (only origination does); speciation re-ids at splits
    assert {e.kind for e in g.events} <= {"origination", "duplication", "speciation"}


def test_loss_can_shrink_and_empty_a_genome():
    # high loss, no origination/duplication except the seeded families -> some lineage loses all
    sp = _tree(seed=6)
    g = simulate_genomes_unordered(sp, loss=2.0, initial_families=4, seed=1)
    sizes = [len(g.genomes[i]) for i in g.genomes]
    assert min(sizes) < 4                                          # at least one node shrank
    assert any(e.kind == "loss" for e in g.events)


def test_duplication_bifurcates_into_two_same_family_children():
    # ZOMBI1 model: a duplication ends the gene and starts two fresh ids descending from it
    sp = _tree(seed=8)
    g = simulate_genomes_unordered(sp, duplication=0.6, initial_families=3, seed=1)
    fam_of = {e.copy: e.family for e in g.events if e.kind != "loss"}   # every gene's birth family
    kids = collections.defaultdict(list)
    for e in g.events:
        if e.kind == "duplication":
            assert e.parent in fam_of and e.family == fam_of[e.parent]  # parent is a real, same-family gene
            kids[e.parent].append(e.copy)
    assert kids and all(len(cs) == 2 for cs in kids.values())          # each duplication has two descendants


def test_every_born_copy_id_is_unique():
    sp = _tree(seed=8, death=0.5)
    g = simulate_genomes_unordered(sp, duplication=0.4, loss=0.3, origination=0.6, initial_families=5, seed=2)
    born = [e.copy for e in g.events if e.kind in ("origination", "duplication")]
    assert len(born) == len(set(born))


def test_family_counts_matches_the_genome():
    sp = _tree(seed=9)
    g = simulate_genomes_unordered(sp, duplication=0.3, loss=0.2, origination=0.5, initial_families=4, seed=1)
    for node_id in g.genomes:
        assert sum(g.family_counts(node_id).values()) == len(g.genomes[node_id])


# --- empty / validation ----------------------------------------------------

def test_empty_run_has_no_events_or_content():
    sp = _tree(seed=1)
    g = simulate_genomes_unordered(sp, seed=1)                            # no families, no rates
    assert g.events == []
    assert all(genome == () for genome in g.genomes.values())


def test_validation():
    sp = _tree(seed=1)
    with pytest.raises(ValueError):
        simulate_genomes_unordered(sp, initial_families=-1, seed=1)
    with pytest.raises(ValueError):
        simulate_genomes_unordered(sp, initial_families=True, seed=1)     # bool is not a valid count
    with pytest.raises(ValueError):
        simulate_genomes_unordered(sp, origination=-1.0, initial_families=1, seed=1)   # negative rate (via scope)


# --- modifiers: OnTime (skyline) is wired; the rest are rejected, not silently dropped ---

def test_time_skyline_modifier_is_supported():
    # OnTime reads only `time`, which the walk supplies, so a skyline origination works: the rate
    # drops to 0 at t=1.5, so no family originates after it
    sp = simulate_species_tree(birth=1.0, death=0.2, total_time=4.0, seed=3)
    r = simulate_genomes_unordered(sp, origination=1.0 * mod.OnTime({0: 1.0, 1.5: 0.0}), seed=1)
    orig_times = [e.time for e in r.events if e.kind == "origination"]
    assert orig_times and max(orig_times) < 1.5


def test_unsupported_modifiers_are_rejected_not_silently_dropped():
    sp = _tree(seed=1)
    # clade drift would need per-lineage threading the walk doesn't do → reject, don't no-op
    with pytest.raises(ValueError, match="does not support"):
        simulate_genomes_unordered(sp, duplication=0.5 * mod.FromParent(spread=0.8), initial_families=3, seed=1)
    # OnTotalDiversity reads a `diversity` context the genome walk doesn't supply → reject, don't crash raw
    with pytest.raises(ValueError, match="does not support"):
        simulate_genomes_unordered(sp, loss=0.25 * mod.OnTotalDiversity(cap=100), initial_families=3, seed=1)


def test_non_default_scope_is_rejected_not_silently_mismatched():
    # a non-default scope sets the total rate one way while the engine still picks the affected
    # copy/lineage the default way — reject it (a PerCopy origination would be base×0 copies, a no-op)
    from zombi2.rates import scope
    sp = _tree(seed=1)
    with pytest.raises(ValueError, match="scope overrides are a later slice"):
        simulate_genomes_unordered(sp, origination=scope.PerCopy(2.0), seed=1)
    with pytest.raises(ValueError, match="scope overrides are a later slice"):
        simulate_genomes_unordered(sp, duplication=scope.PerLineage(0.5), initial_families=3, seed=1)
    with pytest.raises(ValueError, match="scope overrides are a later slice"):
        simulate_genomes_unordered(sp, loss=scope.Global(0.3), initial_families=3, seed=1)
    # the defaults — bare number and the explicit default scope — are accepted
    simulate_genomes_unordered(sp, origination=scope.PerLineage(0.5), duplication=scope.PerCopy(0.5),
                               initial_families=1, seed=1)


# --- transfer: horizontal moves between contemporaneous lineages -----------

def _alive_at(tree, node_id, t):
    n = tree.nodes[node_id]
    return n.birth_time <= t <= n.end_time


def test_transfer_events_are_contemporaneous_donor_to_recipient():
    sp = _tree(seed=3, death=0.4)
    g = simulate_genomes_unordered(sp, transfer=0.4, origination=0.5, initial_families=4, seed=1)
    xfers = _transfers(g.events)
    assert xfers
    for donor, recipient, t in xfers:
        assert donor != recipient                               # a different lineage (default)
        assert _alive_at(sp.complete_tree, donor, t)            # donor alive at t
        assert _alive_at(sp.complete_tree, recipient, t)        # recipient alive at t


def test_transfer_copy_descends_from_a_real_donor_copy():
    sp = _tree(seed=8)
    g = simulate_genomes_unordered(sp, transfer=0.5, initial_families=5, seed=2)
    born = {e.copy for e in g.events if e.kind != "loss"}       # every gene id that was born
    xfer_rows = [e for e in g.events if e.kind == "transfer"]
    assert xfer_rows
    for e in xfer_rows:
        assert e.parent in born and e.copy in born             # the donor gene and the new copy are real


def test_only_transfer_events_carry_a_recipient():
    sp = _tree(seed=2)
    g = simulate_genomes_unordered(sp, duplication=0.2, transfer=0.3, loss=0.2, origination=0.4,
                                   initial_families=4, seed=1)
    for e in g.events:
        if e.recipient is not None:
            assert e.kind == "transfer"                        # only a transfer names a recipient


def test_no_transfer_at_zero_rate():
    sp = _tree(seed=1)
    g = simulate_genomes_unordered(sp, duplication=0.3, loss=0.2, origination=0.4, initial_families=5, seed=1)
    assert all(e.kind != "transfer" for e in g.events)


def test_transfer_is_deterministic():
    sp = _tree(seed=2)
    kw = dict(duplication=0.2, transfer=0.4, loss=0.2, origination=0.4, initial_families=5, seed=9)
    a, b = simulate_genomes_unordered(sp, **kw), simulate_genomes_unordered(sp, **kw)
    assert [str(e) for e in a.events] == [str(e) for e in b.events]
    assert a.genomes == b.genomes


def test_replacement_can_displace_a_resident():
    # with loss=0 the only way a copy is lost is a replacement transfer overwriting a homologous copy;
    # each such loss sits at the same instant as a transfer
    sp = _tree(seed=6)
    g = simulate_genomes_unordered(sp, transfer=1.0, loss=0.0, replacement=True, initial_families=8, seed=1)
    losses = [e for e in g.events if e.kind == "loss"]
    xfer_times = {e.time for e in g.events if e.kind == "transfer"}
    assert losses                                               # replacement did displace some copies
    assert all(e.time in xfer_times for e in losses)           # every loss co-occurs with a transfer


def test_additive_transfer_never_loses():
    sp = _tree(seed=6)
    g = simulate_genomes_unordered(sp, transfer=1.0, loss=0.0, replacement=False, initial_families=8, seed=1)
    assert all(e.kind != "loss" for e in g.events)             # additive transfer only ever adds


def test_default_transfer_is_never_self_but_self_transfer_runs():
    sp = _tree(seed=4)
    g = simulate_genomes_unordered(sp, transfer=0.6, initial_families=5, seed=1)
    assert _transfers(g.events)
    assert all(donor != recipient for donor, recipient, _ in _transfers(g.events))
    s = simulate_genomes_unordered(sp, transfer=0.6, self_transfer=True, initial_families=5, seed=1)
    assert any(e.kind == "transfer" for e in s.events)         # runs (self donor==recipient now allowed)


def test_distance_mode_runs_and_is_deterministic():
    sp = _tree(seed=7, death=0.4)
    from zombi2.genomes import Distance
    a = simulate_genomes_unordered(sp, transfer=0.5, transfer_to="distance", initial_families=5, seed=3)
    b = simulate_genomes_unordered(sp, transfer=0.5, transfer_to=Distance(decay=1.0), initial_families=5, seed=3)
    assert [str(e) for e in a.events] == [str(e) for e in b.events]  # "distance" == Distance(decay=1.0)
    assert any(e.kind == "transfer" for e in a.events)


def test_transfer_can_come_from_the_dead():
    # high death + transfer: some transfers are donated by a lineage that later goes extinct
    sp = _tree(seed=3, death=0.7)
    g = simulate_genomes_unordered(sp, transfer=0.8, origination=0.6, initial_families=4, seed=2)
    donor_fates = {sp.complete_tree.nodes[e.lineage].fate for e in g.events if e.kind == "transfer"}
    assert "extinct" in donor_fates


def test_transfer_to_validation():
    sp = _tree(seed=1)
    with pytest.raises(ValueError, match="transfer_to"):
        simulate_genomes_unordered(sp, transfer=0.3, transfer_to="bogus", initial_families=3, seed=1)


def test_distance_decay_validation():
    from zombi2.genomes import Distance
    Distance(decay=0.0)                       # zero is fine — the uniform limit
    for bad in (-1.0, float("inf"), float("nan"), True):
        with pytest.raises(ValueError, match="Distance decay"):
            Distance(decay=bad)


# --- the written outputs ---------------------------------------------------

def _rows(path):
    lines = path.read_text().splitlines()
    cols = lines[0].split("\t")
    return cols, [dict(zip(cols, ln.split("\t"))) for ln in lines[1:] if ln]


def test_written_node_columns_carry_the_n_label():
    # species_events.tsv and trait_values.tsv have always written n<id>; the genome tables used to
    # write bare ints, so the same node read two ways in one output directory.
    sp = _tree(seed=2)
    g = simulate_genomes_unordered(sp, duplication=0.3, transfer=0.3, loss=0.2, origination=0.5,
                                   initial_families=4, seed=7)
    with tempfile.TemporaryDirectory() as d:
        out = pathlib.Path(d)
        g.write(out, outputs=("events", "genomes"))
        _, events = _rows(out / "genome_events.tsv")
        assert events and all(r["lineage"].startswith("n") for r in events)
        assert any(r["recipient"].startswith("n") for r in events if r["recipient"])
        # gene-copy columns are NOT lineages and stay bare
        assert all(not r["copy"].startswith("n") for r in events)
        assert all(not r["parent"].startswith("n") for r in events if r["parent"])
        _, genomes = _rows(out / "genomes.tsv")
        assert genomes and all(r["lineage"].startswith("n") for r in genomes)


def test_written_log_round_trips_through_the_reader():
    from zombi2.genomes.events import events_from_tsv, events_tsv
    sp = _tree(seed=2)
    g = simulate_genomes_unordered(sp, duplication=0.3, transfer=0.3, loss=0.2, origination=0.5,
                                   initial_families=4, seed=7)
    back = events_from_tsv(events_tsv(g.events))
    assert [(e.time, e.kind, e.lineage, e.recipient) for e in back] == \
           [(e.time, e.kind, e.lineage, e.recipient) for e in g.events]


def test_write_genomes_covers_every_node_where_profiles_covers_only_tips():
    sp = _tree(seed=5, death=0.6)
    g = simulate_genomes_unordered(sp, duplication=0.3, loss=0.3, origination=0.6,
                                   initial_families=5, seed=3)
    with tempfile.TemporaryDirectory() as d:
        out = pathlib.Path(d)
        g.write(out, outputs=("genomes", "profiles"))
        _, rows = _rows(out / "genomes.tsv")
        written = {r["lineage"] for r in rows}
        assert written == {f"n{s}" for s in g.genomes if g.genomes[s]}
        internal = {f"n{n.id}" for n in sp.complete_tree.nodes.values() if n.children is not None}
        assert written & internal, "ancestral genomes must be in there, not just the tips"
        # profiles is the extant-only view
        tips = {f"n{n.id}" for n in sp.complete_tree.extant()}
        assert set((out / "profiles.tsv").read_text().splitlines()[0].split("\t")[1:]) == tips


def test_write_gene_trees_emits_one_newick_per_family():
    sp = _tree(seed=5, death=0.6)
    g = simulate_genomes_unordered(sp, duplication=0.3, loss=0.3, origination=0.6,
                                   initial_families=5, seed=3)
    with tempfile.TemporaryDirectory() as d:
        out = pathlib.Path(d)
        g.write(out, outputs=("gene_trees",))
        for fam, gt in g.gene_trees.items():
            complete = out / f"gene_tree_fam{fam}_complete.nwk"
            assert complete.read_text().strip() == gt.to_newick("complete")
            extant = out / f"gene_tree_fam{fam}_extant.nwk"
            # a family with no survivor has no extant tree, and writes no file for it
            assert extant.exists() == (gt.to_newick("extant") is not None)


# --- ByFamily: per-family rate heterogeneity -------------------------------

def _dup_per_family(g, n_families):
    counts = collections.Counter(e.family for e in g.events if e.kind == "duplication")
    return [counts.get(f, 0) for f in range(n_families)]


def test_by_family_spreads_the_rates_without_moving_their_mean():
    # the point of the modifier: families stop being interchangeable. The draw is mean-corrected,
    # so widening the spread must widen the spread of outcomes without inflating the average.
    sp = _tree(seed=1, n_extant=20, death=0.0)
    flat = simulate_genomes_unordered(sp, duplication=0.25, loss=0.25, initial_families=150, seed=3)
    varied = simulate_genomes_unordered(sp, duplication=0.25 * mod.ByFamily(spread=0.5),
                                        loss=0.25, initial_families=150, seed=3)
    f, v = _dup_per_family(flat, 150), _dup_per_family(varied, 150)
    import statistics
    assert statistics.pstdev(v) > 1.5 * statistics.pstdev(f)      # families genuinely differ
    assert statistics.mean(v) == pytest.approx(statistics.mean(f), rel=0.35)   # mean is held


def test_a_run_with_no_by_family_is_untouched():
    # the weighted path costs something, so it must only be taken when it is asked for
    sp = _tree(seed=2, n_extant=12)
    a = simulate_genomes_unordered(sp, duplication=0.2, transfer=0.1, loss=0.2,
                                   initial_families=10, seed=5)
    b = simulate_genomes_unordered(sp, duplication=0.2, transfer=0.1, loss=0.2,
                                   initial_families=10, seed=5)
    assert [(e.time, e.kind, e.copy) for e in a.events] == [(e.time, e.kind, e.copy) for e in b.events]


def test_by_family_is_deterministic_given_the_seed():
    sp = _tree(seed=2, n_extant=12)
    kw = dict(duplication=0.2 * mod.ByFamily(spread=0.6), loss=0.2, initial_families=20, seed=5)
    a = simulate_genomes_unordered(sp, **kw)
    b = simulate_genomes_unordered(sp, **kw)
    assert [(e.time, e.kind, e.copy) for e in a.events] == [(e.time, e.kind, e.copy) for e in b.events]


def test_family_speed_moves_every_rate_of_a_family_together():
    # the other placement: one draw per family, scaling all its rates. A family that duplicates a
    # lot should also be losing a lot — which is exactly what a per-rate ByFamily does NOT give.
    sp = _tree(seed=1, n_extant=20, death=0.0)
    g = simulate_genomes_unordered(sp, duplication=0.25, loss=0.25,
                                   family_speed=mod.ByFamily(spread=0.6),
                                   initial_families=150, seed=3)
    dup = collections.Counter(e.family for e in g.events if e.kind == "duplication")
    los = collections.Counter(e.family for e in g.events if e.kind == "loss")
    fams = [f for f in range(150) if dup.get(f, 0) + los.get(f, 0) > 4]
    assert len(fams) > 20
    xs = [dup.get(f, 0) for f in fams]
    ys = [los.get(f, 0) for f in fams]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    cov = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sx = sum((a - mx) ** 2 for a in xs) ** 0.5
    sy = sum((b - my) ** 2 for b in ys) ** 0.5
    assert cov / (sx * sy) > 0.3          # a fast family is fast at everything


def test_by_family_is_refused_on_origination():
    sp = _tree(seed=1, n_extant=8)
    with pytest.raises(ValueError, match="families are CREATED"):
        simulate_genomes_unordered(sp, origination=0.5 * mod.ByFamily(spread=0.3), seed=1)


def test_family_speed_takes_a_by_family_modifier():
    sp = _tree(seed=1, n_extant=8)
    with pytest.raises(ValueError, match="family_speed takes a ByFamily"):
        simulate_genomes_unordered(sp, loss=0.2, family_speed=0.5, seed=1)


def test_by_family_with_driven_by_is_refused_for_now():
    sp = _tree(seed=1, n_extant=8)
    with pytest.raises(ValueError, match="later slice"):
        simulate_genomes_unordered(sp, loss=0.2 * mod.ByFamily(spread=0.3),
                                   duplication=0.2 * mod.DrivenBy("x.tsv", {"a": 2.0}), seed=1)


# --- max_family_size: a per-genome ceiling on a family's copies ------------

def _biggest_family(g):
    return max((collections.Counter(c.family for c in gen).most_common(1) or [(None, 0)])[0][1]
               for gen in g.genomes.values() if gen)


def test_max_family_size_binds_exactly():
    # duplication far above loss compounds without bound; the cap is what stops it, and it stops it
    # AT the number given rather than somewhere near it
    sp = _tree(seed=1, n_extant=12, death=0.0)
    for cap in (2, 5, 9):
        g = simulate_genomes_unordered(sp, duplication=0.9, loss=0.05, initial_families=6,
                                       max_family_size=cap, seed=2)
        assert _biggest_family(g) == cap


def test_max_family_size_none_lifts_the_ceiling():
    sp = _tree(seed=1, n_extant=12, death=0.0)
    capped = simulate_genomes_unordered(sp, duplication=0.9, loss=0.05, initial_families=6,
                                        max_family_size=5, seed=2)
    free = simulate_genomes_unordered(sp, duplication=0.9, loss=0.05, initial_families=6,
                                      max_family_size=None, seed=2)
    assert _biggest_family(free) > _biggest_family(capped)


def test_a_float_cap_scales_with_the_tree():
    # the point of the float form: the bound travels with the size of the run
    from zombi2.genomes import resolve_max_family_size
    small = simulate_species_tree(birth=1.0, death=0.0, n_extant=6, seed=1)
    big = simulate_species_tree(birth=1.0, death=0.0, n_extant=40, seed=1)
    a = resolve_max_family_size(2.0, len(small.complete_tree.nodes))
    b = resolve_max_family_size(2.0, len(big.complete_tree.nodes))
    assert b > a
    assert resolve_max_family_size(7, 999) == 7          # an int is absolute, tree size ignored
    assert resolve_max_family_size(None, 999) is None


def test_max_family_size_is_validated():
    from zombi2.genomes import resolve_max_family_size
    for bad in (0, -1, -2.0, True, "big"):
        with pytest.raises(ValueError):
            resolve_max_family_size(bad, 100)


def test_the_cap_also_holds_when_a_transfer_arrives():
    # a transfer adds a copy to the recipient, so the ceiling has to hold there too or a family
    # could be pushed past it sideways
    sp = _tree(seed=3, n_extant=15, death=0.0)
    g = simulate_genomes_unordered(sp, duplication=0.4, transfer=0.6, loss=0.05,
                                   initial_families=8, max_family_size=4, seed=1)
    assert _biggest_family(g) <= 4
