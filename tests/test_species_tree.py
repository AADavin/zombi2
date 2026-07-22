"""Tests for the forward birth-death engine (species)."""

import re

import pytest

from zombi2.rates import modifiers as mod
from zombi2.rates import scope
from zombi2.species import Event, simulate_species_tree


def test_yule_reaches_n_extant_with_no_extinction():
    r = simulate_species_tree(birth=1.0, death=0.0, n_extant=50, seed=1)
    assert r.n_extant == 50
    assert all(e.kind == "speciation" for e in r.events)  # Yule → no deaths
    assert all(n.fate == "extant" for n in r.complete_tree.leaves())


def test_birth_death_has_extinctions_and_survivors():
    r = simulate_species_tree(birth=1.0, death=0.4, n_extant=60, seed=7)
    assert r.n_extant == 60
    assert "extinction" in {e.kind for e in r.events}
    assert len(r.complete_tree.extinct()) > 0


def test_total_time_stop_ends_extant_lineages_at_the_present():
    r = simulate_species_tree(birth=1.0, death=0.2, total_time=4.0, seed=3)
    for n in r.complete_tree.extant():
        assert n.end_time == pytest.approx(4.0)


def test_deterministic_given_seed():
    a = simulate_species_tree(birth=1.0, death=0.3, n_extant=40, seed=42)
    b = simulate_species_tree(birth=1.0, death=0.3, n_extant=40, seed=42)
    assert len(a.complete_tree.nodes) == len(b.complete_tree.nodes)
    assert [(e.time, e.kind, e.node) for e in a.events] == [(e.time, e.kind, e.node) for e in b.events]


def test_different_seeds_differ():
    a = simulate_species_tree(birth=1.0, death=0.3, n_extant=40, seed=1)
    b = simulate_species_tree(birth=1.0, death=0.3, n_extant=40, seed=2)
    assert len(a.complete_tree.nodes) != len(b.complete_tree.nodes) or \
        [e.time for e in a.events] != [e.time for e in b.events]


def test_tree_structure_invariants():
    r = simulate_species_tree(birth=1.0, death=0.2, n_extant=30, seed=5)
    for node in r.complete_tree.nodes.values():
        if node.children is None:
            assert node.fate in ("extant", "extinct")
        else:
            c1, c2 = node.children
            assert c1 != c2
            assert r.complete_tree.nodes[c1].parent == node.id
            assert r.complete_tree.nodes[c2].parent == node.id
            assert node.fate == "speciation"


def test_events_are_time_ordered():
    r = simulate_species_tree(birth=1.0, death=0.3, n_extant=50, seed=9)
    times = [e.time for e in r.events]
    assert times == sorted(times)


def test_n_extant_conditions_on_survival():
    # seed 7's first attempt dies on the first event; conditioning restarts until 60 living
    r = simulate_species_tree(birth=1.0, death=0.4, n_extant=60, seed=7)
    assert r.n_extant == 60


def test_raises_when_n_extant_unreachable():
    with pytest.raises(RuntimeError):
        simulate_species_tree(birth=0.1, death=3.0, n_extant=200, seed=1)


# --- the scope × modifiers grammar actually drives the engine -------------

def test_global_grows_slower_than_per_lineage():
    # per-lineage birth compounds (exponential); Global is a constant tree-wide budget (linear)
    per = simulate_species_tree(birth=1.0, death=0.0, total_time=4.0, seed=1)
    glob = simulate_species_tree(birth=scope.Global(2.0), death=0.0, total_time=4.0, seed=1)
    assert per.n_extant > glob.n_extant


def test_diversity_caps_growth():
    # OnTotalDiversity(cap=20): the birth factor falls to 0 at 20 lineages, so the tree saturates
    r = simulate_species_tree(birth=1.0 * mod.OnTotalDiversity(cap=20), death=0.0, total_time=100.0, seed=1)
    assert r.n_extant <= 20   # never exceeds the cap
    assert r.n_extant >= 15   # but it grew toward it


def test_no_scope_default_is_per_lineage():
    # a bare number and an explicit PerLineage wrapper give the same tree
    a = simulate_species_tree(birth=1.0, death=0.2, n_extant=40, seed=11)
    b = simulate_species_tree(birth=scope.PerLineage(1.0), death=scope.PerLineage(0.2), n_extant=40, seed=11)
    assert [(e.time, e.kind) for e in a.events] == [(e.time, e.kind) for e in b.events]


# --- validation -----------------------------------------------------------

def test_validation():
    with pytest.raises(ValueError):
        simulate_species_tree(birth=1.0)                          # neither n_extant nor total_time
    with pytest.raises(ValueError):
        simulate_species_tree(birth=1.0, n_extant=10, total_time=5.0)    # both
    with pytest.raises(ValueError):
        simulate_species_tree(birth=-1.0, n_extant=10)            # negative rate (via scope)
    with pytest.raises(TypeError):
        simulate_species_tree(birth="fast", n_extant=10)          # non-numeric rate
    with pytest.raises(ValueError):
        simulate_species_tree(birth=1.0, n_extant=0)              # non-positive n_extant
    with pytest.raises(ValueError):
        simulate_species_tree(birth=1.0, total_time=-2.0)               # non-positive total_time


def test_event_is_frozen_record():
    e = Event(1.0, "extinction", 3)
    with pytest.raises(Exception):
        e.time = 2.0  # type: ignore[misc]


def test_skyline_stops_births_after_a_zero_breakpoint():
    # birth 1.0 on [0, 2), then 0 → the interval-aware sampler must forbid births at/after t=2
    r = simulate_species_tree(birth=1.0 * mod.OnTime({0: 1.0, 2.0: 0.0}), death=0.0, total_time=10.0, seed=1)
    spec_times = [e.time for e in r.events if e.kind == "speciation"]
    assert spec_times                 # growth happened before the breakpoint
    assert max(spec_times) < 2.0      # and nothing after the rate dropped to 0


def test_skyline_is_deterministic():
    kw = dict(birth=1.0 * mod.OnTime({0: 2.0, 3: 0.2}), death=0.1, total_time=6.0, seed=4)
    a = simulate_species_tree(**kw)
    b = simulate_species_tree(**kw)
    assert [(e.time, e.kind) for e in a.events] == [(e.time, e.kind) for e in b.events]


# --- the extant tree + Newick output --------------------------------------

def test_extant_tree_prunes_to_survivors():
    r = simulate_species_tree(birth=1.0, death=0.4, n_extant=40, seed=3)
    ext = r.extant_tree
    assert len(ext.leaves()) == r.n_extant                       # exactly the survivors
    assert all(n.fate == "extant" for n in ext.leaves())         # no extinct tips
    assert all(n.fate != "extinct" for n in ext.nodes.values())  # no extinct nodes remain
    for n in ext.nodes.values():                                 # bifurcating
        assert n.children is None or len(n.children) == 2
    # every branch is now strictly positive — the present sits after the last split (n_extant fix)
    for n in ext.nodes.values():
        assert n.end_time - n.birth_time > 0


def test_yule_extant_equals_complete_leaves():
    r = simulate_species_tree(birth=1.0, death=0.0, n_extant=30, seed=1)
    assert len(r.extant_tree.leaves()) == len(r.complete_tree.leaves()) == 30


def test_newick_is_wellformed():
    r = simulate_species_tree(birth=1.0, death=0.3, n_extant=25, seed=2)
    nwk = r.extant_tree.to_newick()
    assert nwk.endswith(";")
    assert nwk.count("(") == nwk.count(")")   # balanced parens
    assert nwk.count(",") == 25 - 1           # a bifurcating tree of 25 tips has 24 joins


def test_newick_root_carries_the_stem():
    # a forward run starts from one lineage, so the root's branch is real simulated time: the stem,
    # origin to first split. Writing it is what keeps the tree's full height in the file.
    r = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=2)
    root = r.complete_tree.nodes[r.complete_tree.root]
    stem = root.end_time - root.birth_time
    first_split = min(e.time for e in r.events if e.kind == "speciation")

    assert stem == pytest.approx(first_split)
    assert re.search(rf"\)n{r.complete_tree.root}:[0-9.e+-]+;$", r.complete_tree.to_newick())
    assert float(r.complete_tree.to_newick().rsplit(":", 1)[1].rstrip(";")) == pytest.approx(
        stem, rel=1e-5)


def test_write_produces_newick_files(tmp_path):
    r = simulate_species_tree(birth=1.0, death=0.2, n_extant=20, seed=5)
    r.write(tmp_path)
    assert (tmp_path / "species_complete.nwk").read_text().strip().endswith(";")
    assert (tmp_path / "species_extant.nwk").read_text().strip().endswith(";")


def test_write_records_the_event_log(tmp_path):
    r = simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=5)
    r.write(tmp_path)                                            # events are always written
    lines = (tmp_path / "species_events.tsv").read_text().splitlines()
    assert lines[0] == "time\tkind\tlineage\tchildren"
    assert len(lines) == 1 + len(r.events)                      # one row per recorded event
    speciation = next(ln for ln in lines[1:] if "\tspeciation\t" in ln)
    kids = speciation.split("\t")[-1]
    assert ";" in kids and kids.count("n") == 2                 # a speciation lists its two children


def test_write_is_selective(tmp_path):
    r = simulate_species_tree(birth=1.0, death=0.3, n_extant=20, fossils=0.5, seed=3)
    r.write(tmp_path, outputs=["extant", "events"])
    assert {p.name for p in tmp_path.iterdir()} == {"species_extant.nwk", "species_events.tsv"}  # only what was asked


def test_write_rejects_unknown_output(tmp_path):
    r = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=1)
    with pytest.raises(ValueError, match="unknown write outputs"):
        r.write(tmp_path, outputs=["complete", "bogus"])


def test_extant_tree_is_deterministic():
    a = simulate_species_tree(birth=1.0, death=0.3, n_extant=30, seed=8)
    b = simulate_species_tree(birth=1.0, death=0.3, n_extant=30, seed=8)
    assert a.extant_tree.to_newick() == b.extant_tree.to_newick()


def test_dead_tree_has_no_extant_tree():
    r = simulate_species_tree(birth=0.1, death=10.0, total_time=5.0, seed=1)
    assert r.n_extant == 0
    assert r.extant_tree is None


# --- FromParent (clade drift): rates drift down the tree, picking is rate-weighted ---

def test_weighted_index_respects_weights():
    import numpy as np

    from zombi2.species import _weighted_index
    rng = np.random.default_rng(0)
    weights = [1.0, 1.0, 8.0]           # index 2 carries 80% of the total rate
    counts = [0, 0, 0]
    for _ in range(20000):
        counts[_weighted_index(rng, weights, sum(weights))] += 1
    assert 0.77 < counts[2] / 20000 < 0.83   # ≈ 0.8
    assert counts[0] > 0 and counts[1] > 0    # the light lineages still get picked sometimes


def test_clade_drift_is_deterministic_given_seed():
    kw = dict(birth=1.0 * mod.FromParent(spread=0.5), death=0.1, n_extant=40, seed=3)
    a = simulate_species_tree(**kw)
    b = simulate_species_tree(**kw)
    assert [(e.time, e.kind, e.node) for e in a.events] == [(e.time, e.kind, e.node) for e in b.events]


def test_inherited_zero_spread_reaches_target():
    # spread 0 → every step is ×1, so no drift; still a valid birth-death tree
    r = simulate_species_tree(birth=1.0 * mod.FromParent(spread=0.0), death=0.2, n_extant=40, seed=5)
    assert r.n_extant == 40


def test_death_can_drift_independently():
    # drift lives on death, not birth; birth and death are bent independently
    r = simulate_species_tree(birth=1.0, death=0.4 * mod.FromParent(spread=0.5), n_extant=50, seed=4)
    assert r.n_extant == 50
    assert len(r.complete_tree.extinct()) > 0


def test_clade_drift_composes_with_diversity_cap():
    # clade drift × diversity-dependence: the cap still bounds the tree
    r = simulate_species_tree(
        birth=1.0 * mod.FromParent(spread=0.4) * mod.OnTotalDiversity(cap=25), death=0.0, total_time=100.0, seed=1)
    assert r.n_extant <= 25          # the cap is a hard ceiling even with drift
    assert r.n_extant >= 12          # and the tree grew toward it


def test_inherited_requires_per_lineage_scope():
    # per-lineage drift on a Global (tree-wide) budget is contradictory — reject it clearly
    with pytest.raises(ValueError, match="per lineage"):
        simulate_species_tree(birth=scope.Global(1.0) * mod.FromParent(spread=0.2), n_extant=10, seed=1)


def test_drifting_birth_with_non_drifting_global_death_is_allowed():
    # only the drifting rate must be per lineage; a Global death budget alongside it is fine
    r = simulate_species_tree(
        birth=1.0 * mod.FromParent(spread=0.3), death=scope.Global(0.2), n_extant=30, seed=2)
    assert r.n_extant == 30


# --- the level rejects what it does not wire (SPEC §5) --------------------

@pytest.mark.parametrize("modifier", [mod.ByLineage(spread=0.3),
                                      mod.DrivenBy("habitat.tsv", {"a": 2.0})])
def test_an_unwired_modifier_raises_rather_than_being_ignored(modifier):
    # an unthreaded modifier returns its default factor of 1.0, so silently accepting one would run
    # a model the user did not ask for — the whole point of declaring WIRED_MODIFIERS
    with pytest.raises(ValueError, match="does not support"):
        simulate_species_tree(birth=1.0 * modifier, n_extant=10, seed=1)
    with pytest.raises(ValueError, match="does not support"):
        simulate_species_tree(birth=1.0, death=0.1 * modifier, n_extant=10, seed=1)


@pytest.mark.parametrize("wrapper", [scope.PerCopy, scope.PerSite, scope.PerChromosome])
def test_a_foreign_scope_raises(wrapper):
    # the species engine counts lineages; any other unit would ask it for a count it has no idea of
    with pytest.raises(ValueError, match="counts lineages"):
        simulate_species_tree(birth=wrapper(1.0), n_extant=10, seed=1)


def _colless(result):
    """Colless imbalance of the extant tree: Σ over internal nodes of |left tips − right tips|.
    Higher = more lopsided. A pure Yule tree is comparatively balanced; heritable rate drift
    concentrates tips in the fast clades, so it climbs."""
    tree = result.extant_tree
    size = {}
    for i in sorted(tree.nodes, reverse=True):        # children (higher ids) before parents
        nd = tree.nodes[i]
        size[i] = 1 if nd.children is None else sum(size[c] for c in nd.children)
    return sum(abs(size[nd.children[0]] - size[nd.children[1]])
               for nd in tree.nodes.values() if nd.children is not None)


def test_clade_drift_is_more_imbalanced_than_yule():
    # the signature of heritable rate drift, at a fixed tip count so it is shape not size: fast
    # clades are inherited, so they hoard the tips and the tree is far more lopsided than Yule
    import statistics
    seeds = range(40)
    yule = [_colless(simulate_species_tree(birth=1.0, death=0.0, n_extant=64, seed=s)) for s in seeds]
    drift = [_colless(simulate_species_tree(birth=1.0 * mod.FromParent(spread=0.9), death=0.0, n_extant=64, seed=s))
             for s in seeds]
    assert statistics.mean(drift) > 1.5 * statistics.mean(yule)   # observed ≈ 2.7× (margin to spare)


# --- mass extinctions: (time, fraction_lost) survival pulses, time forward from the crown, total_time mode ---

def test_mass_extinction_culls_diversity():
    import statistics
    seeds = range(30)
    no_pulse = [simulate_species_tree(birth=1.0, death=0.2, total_time=5.0, seed=s).n_extant for s in seeds]
    culled = [simulate_species_tree(birth=1.0, death=0.2, total_time=5.0, mass_extinctions=[(3.0, 0.9)], seed=s).n_extant
              for s in seeds]
    # a 90% cull at time 3.0 leaves far fewer survivors even after some regrowth to the present
    assert statistics.mean(culled) < 0.5 * statistics.mean(no_pulse)


def test_total_mass_extinction_wipes_the_tree():
    r = simulate_species_tree(birth=1.0, death=0.2, total_time=5.0, mass_extinctions=[(2.5, 1.0)], seed=1)
    assert r.n_extant == 0            # fraction lost = 1.0 kills every standing lineage
    assert r.extant_tree is None


def test_mass_extinction_deaths_land_at_the_pulse_instant():
    # a pulse at time 2.0 (forward from the crown) puts its deaths exactly there
    r = simulate_species_tree(birth=1.0, death=0.2, total_time=5.0, mass_extinctions=[(2.0, 0.75)], seed=3)
    culled = [e for e in r.events if e.kind == "extinction" and e.time == pytest.approx(2.0)]
    assert len(culled) > 0
    assert all(r.complete_tree.nodes[e.node].fate == "extinct" for e in culled)


def test_multiple_mass_extinctions_each_fire():
    r = simulate_species_tree(birth=1.2, death=0.1, total_time=6.0,
                              mass_extinctions=[(4.0, 0.5), (2.0, 0.5)], seed=2)
    times = {e.time for e in r.events if e.kind == "extinction"}
    assert any(t == pytest.approx(2.0) for t in times)   # each pulse fires at its own time
    assert any(t == pytest.approx(4.0) for t in times)


def test_zero_fraction_pulse_kills_nobody():
    r = simulate_species_tree(birth=1.0, death=0.2, total_time=5.0, mass_extinctions=[(2.0, 0.0)], seed=3)
    at_instant = [e for e in r.events if e.kind == "extinction" and e.time == pytest.approx(2.0)]
    assert at_instant == []          # survival 1.0 → the pulse removes no one


def test_mass_extinction_is_deterministic():
    kw = dict(birth=1.0, death=0.3, total_time=5.0, mass_extinctions=[(2.0, 0.6)], seed=7)
    a = simulate_species_tree(**kw)
    b = simulate_species_tree(**kw)
    assert [(e.time, e.kind, e.node) for e in a.events] == [(e.time, e.kind, e.node) for e in b.events]


def test_mass_extinction_requires_total_time():
    with pytest.raises(ValueError, match="total_time"):
        simulate_species_tree(birth=1.0, n_extant=10, mass_extinctions=[(1.0, 0.5)], seed=1)


def test_mass_extinction_validation():
    with pytest.raises(ValueError):
        simulate_species_tree(birth=1.0, total_time=5.0, mass_extinctions=[(5.0, 0.5)], seed=1)   # time not < total_time
    with pytest.raises(ValueError):
        simulate_species_tree(birth=1.0, total_time=5.0, mass_extinctions=[(0.0, 0.5)], seed=1)   # time not > 0
    with pytest.raises(ValueError):
        simulate_species_tree(birth=1.0, total_time=5.0, mass_extinctions=[(2.0, 1.5)], seed=1)   # fraction > 1
    with pytest.raises(ValueError):
        simulate_species_tree(birth=1.0, total_time=5.0, mass_extinctions=[(2.0, -0.1)], seed=1)  # fraction < 0


# --- fossils: Poisson(rate × branch length) side-output along the complete tree ---

def test_fossils_scale_with_rate():
    low = simulate_species_tree(birth=1.0, death=0.4, n_extant=40, fossils=0.3, seed=3)
    high = simulate_species_tree(birth=1.0, death=0.4, n_extant=40, fossils=1.0, seed=3)
    assert len(high.fossils) > len(low.fossils) > 0


def test_no_fossils_at_zero_rate():
    r = simulate_species_tree(birth=1.0, death=0.4, n_extant=40, fossils=0.0, seed=3)
    assert r.fossils == []


def test_fossil_times_lie_within_their_branch():
    r = simulate_species_tree(birth=1.0, death=0.4, n_extant=40, fossils=0.5, seed=3)
    assert r.fossils
    for lineage, t in r.fossils:
        node = r.complete_tree.nodes[lineage]
        assert node.birth_time <= t <= node.end_time


def test_fossils_are_sorted_by_time():
    r = simulate_species_tree(birth=1.0, death=0.4, n_extant=40, fossils=0.5, seed=3)
    assert [t for _, t in r.fossils] == sorted(t for _, t in r.fossils)


def test_fossils_do_not_change_the_tree():
    # a pure side output: the tree is grown before fossils are drawn, so it is bit-identical
    plain = simulate_species_tree(birth=1.0, death=0.4, n_extant=40, seed=3)
    withf = simulate_species_tree(birth=1.0, death=0.4, n_extant=40, fossils=0.7, seed=3)
    assert withf.n_extant == plain.n_extant
    assert withf.extant_tree.to_newick() == plain.extant_tree.to_newick()


def test_fossils_are_deterministic():
    kw = dict(birth=1.0, death=0.4, n_extant=40, fossils=0.5, seed=9)
    assert simulate_species_tree(**kw).fossils == simulate_species_tree(**kw).fossils


def test_fossils_recovered_along_extinct_branches_too():
    # fossils fall along ALL branches of the complete tree, so extinct lineages get fossils
    r = simulate_species_tree(birth=1.0, death=0.4, n_extant=40, fossils=1.0, seed=3)
    fates = {r.complete_tree.nodes[i].fate for i, _ in r.fossils}
    assert "extinct" in fates


def test_fossils_validation():
    with pytest.raises(ValueError):
        simulate_species_tree(birth=1.0, n_extant=10, fossils=-0.1, seed=1)


def test_fossils_write_tsv(tmp_path):
    r = simulate_species_tree(birth=1.0, death=0.4, n_extant=30, fossils=0.5, seed=3)
    r.write(tmp_path)
    lines = (tmp_path / "species_fossils.tsv").read_text().splitlines()
    assert lines[0] == "lineage\ttime"
    assert len(lines) == 1 + len(r.fossils)
    # no fossils file when none were recovered
    simulate_species_tree(birth=1.0, death=0.4, n_extant=30, seed=3).write(tmp_path / "nof")
    assert not (tmp_path / "nof" / "species_fossils.tsv").exists()


# --- incomplete sampling (rho): observe a fraction of the survivors ---

def _survivor_ids(result):
    return {n.id for n in result.complete_tree.extant()} | {n.id for n in result.complete_tree.unsampled()}


def test_sampling_relabels_not_removes():
    # n_extant stops at 40 SURVIVORS; sampling then splits them into extant + unsampled
    r = simulate_species_tree(birth=1.0, death=0.3, n_extant=40, sampling=0.5, seed=3)
    assert len(r.complete_tree.extant()) + len(r.complete_tree.unsampled()) == 40
    assert 0 < r.n_extant < 40                       # some observed, some not


def test_extant_tree_is_the_sampled_survivors():
    r = simulate_species_tree(birth=1.0, death=0.3, n_extant=40, sampling=0.5, seed=3)
    assert len(r.extant_tree.leaves()) == r.n_extant                 # the extant tree is the observed one
    assert all(n.fate == "extant" for n in r.extant_tree.nodes.values() if n.children is None)
    for n in r.extant_tree.nodes.values():                           # still bifurcating after pruning
        assert n.children is None or len(n.children) == 2


def test_sampling_one_observes_everyone():
    r = simulate_species_tree(birth=1.0, death=0.3, n_extant=40, sampling=1.0, seed=3)
    assert r.complete_tree.unsampled() == []
    assert r.n_extant == 40


def test_sampling_does_not_change_the_grown_tree():
    # sampling only relabels survivors after growth, so the grown survivor set is identical
    full = simulate_species_tree(birth=1.0, death=0.3, n_extant=40, sampling=1.0, seed=3)
    half = simulate_species_tree(birth=1.0, death=0.3, n_extant=40, sampling=0.5, seed=3)
    assert _survivor_ids(full) == _survivor_ids(half)


def test_sampling_fraction_matches_rho():
    import statistics
    seeds = range(30)
    fracs = []
    for s in seeds:
        r = simulate_species_tree(birth=1.0, death=0.3, n_extant=40, sampling=0.5, seed=s)
        fracs.append(r.n_extant / 40)
    assert 0.45 < statistics.mean(fracs) < 0.55       # observed fraction ≈ ρ = 0.5


def test_sampling_is_deterministic():
    kw = dict(birth=1.0, death=0.3, n_extant=40, sampling=0.5, seed=9)
    a = {n.id for n in simulate_species_tree(**kw).complete_tree.extant()}
    b = {n.id for n in simulate_species_tree(**kw).complete_tree.extant()}
    assert a == b


def test_sampling_validation():
    for bad in (0.0, 1.5, -0.1):
        with pytest.raises(ValueError):
            simulate_species_tree(birth=1.0, n_extant=10, sampling=bad, seed=1)
