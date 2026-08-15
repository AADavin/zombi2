"""A continuously diffusing trait driving speciation — QuaSSE (design note §8).

This is the one joint model that is **not** exact. Every other driver here changes only at events,
and an event ends the Gillespie step, so the rate is constant in between. A diffusion moves at every
instant. The run therefore slices: the value is held fixed across `step` and released at the
boundary, where each lineage takes the exact transition law of its own diffusion. So the trait is
exact and only its coupling to speciation is approximated.
"""

import math
import statistics

import pytest

from zombi2 import joint, species, traits
from zombi2.params import Curve, PerLineage, Scalar
from zombi2.traits import TraitsResult


def _run(slope=0.5, step=0.05, seed=1, n_extant=50, birth=0.4, death=0.05, **trait_kw):
    """Bigger animals speciate faster, and body size diffuses."""
    return joint.simulate(
        species.birth_death(
            birth=PerLineage(birth).scaled_by(
                "trait", Curve(lambda x, a=slope: math.exp(a * x)), step=step),
            death=death, n_extant=n_extant),
        traits.continuous(start=0.0, rate=1.0, **trait_kw), seed=seed)


def _tip_mean(r):
    tips = sorted(r.complete_tree.extant_leaves())
    return statistics.fmean(r.trait.node_values[i] for i in tips)


# --- what comes back ------------------------------------------------------------------------------

def test_it_returns_a_tree_and_an_ordinary_continuous_trait():
    r = _run()
    ct = r.complete_tree
    assert sum(1 for n in ct.nodes.values() if n.fate == "extant") == 50
    assert isinstance(r.trait, TraitsResult) and r.trait.kind == "continuous"
    assert set(r.trait.node_values) == set(ct.nodes)
    assert all(isinstance(v, float) and math.isfinite(v) for v in r.trait.node_values.values())


def test_the_log_opens_with_the_value_at_time_zero():
    """A diffusion cannot be rebuilt from events, so the log carries what the trait level's own log
    carries — the origin, and any jump at a split. The values ride in `node_values`."""
    r = _run()
    first = r.trait.events[0]
    assert (first.kind, first.lineage, first.from_state, first.to_state) == \
           ("initial", r.complete_tree.root, None, 0.0)
    assert {c.kind for c in r.trait.events} == {"initial"}


def test_it_is_deterministic():
    a, b = _run(seed=7), _run(seed=7)
    assert a.trait.node_values == b.trait.node_values
    assert a.complete_tree.to_newick() == b.complete_tree.to_newick()


def test_total_time_stops_the_run_there():
    r = joint.simulate(
        species.birth_death(
            birth=PerLineage(0.5).scaled_by("trait", Curve(lambda x: math.exp(0.3 * x)), step=0.05),
            death=0.1, total_time=6.0),
        traits.continuous(rate=1.0), seed=3)
    ends = [n.end_time for n in r.complete_tree.nodes.values() if n.fate == "extant"]
    assert ends and all(e == pytest.approx(6.0) for e in ends)


def test_both_levels_write():
    """A joint run's claim is that each level lands on disk exactly as its own command writes it."""
    import pathlib
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        _run(seed=2).write(d)
        names = {p.name for p in pathlib.Path(d).rglob("*") if p.is_file()}
    assert "trait_values.tsv" in names and "species_complete.nwk" in names


# --- the model actually does what it says ---------------------------------------------------------

def test_a_bigger_value_speciating_faster_leaves_the_tips_bigger():
    """The QuaSSE signature, and the thing that makes the run joint at all.

    Under plain Brownian motion the tips average the root value, 0. Let the birth rate rise with the
    value and the standing lineages are a biased sample of the diffusion — the fast ones left more
    descendants — so the tips average well above where they started. The flat curve is the same run
    with only the exponent changed, so it controls for everything else.
    """
    driven = statistics.fmean(_tip_mean(_run(slope=0.5, seed=s)) for s in range(1, 21))
    flat = statistics.fmean(_tip_mean(_run(slope=0.0, seed=s)) for s in range(1, 21))
    assert abs(flat) < 1.0, f"an undriven run should sit near the root value, got {flat:+.3f}"
    assert driven > flat + 1.5, f"driven {driven:+.3f} against flat {flat:+.3f}"


def test_a_scalar_log_link_drives_it_too():
    """`Curve` is not the only continuous mapping: a `Scalar` is the log-link the SSE literature
    writes, and it reaches the same engine by the same route."""
    r = joint.simulate(
        species.birth_death(
            birth=PerLineage(0.4).scaled_by("trait", Scalar(0.5), step=0.05),
            death=0.05, n_extant=50),
        traits.continuous(rate=1.0), seed=4)
    assert sum(1 for n in r.complete_tree.nodes.values() if n.fate == "extant") == 50


def test_the_answer_holds_up_when_the_step_is_halved():
    """The check the manual asks for. Slicing is the approximation, so a run whose answer moves when
    `step` is halved was read at too coarse a resolution."""
    coarse = statistics.fmean(_tip_mean(_run(step=0.05, seed=s)) for s in range(1, 31))
    fine = statistics.fmean(_tip_mean(_run(step=0.025, seed=s)) for s in range(1, 31))
    assert abs(coarse - fine) < 0.6, f"step=0.05 gave {coarse:+.3f}, step=0.025 gave {fine:+.3f}"


def test_reverting_to_an_optimum_holds_the_trait_near_it():
    """Ornstein–Uhlenbeck, over slices: each boundary applies the exact OU transition, so the pull
    is not something the slicing approximates."""
    vals = []
    for s in range(1, 11):
        r = _run(slope=0.3, seed=s, reverts_to=4.0, pull=1.5)
        tips = sorted(r.complete_tree.extant_leaves())
        vals.extend(r.trait.node_values[i] for i in tips)
    assert statistics.fmean(vals) == pytest.approx(4.0, abs=1.0)
    assert statistics.stdev(vals) < 1.5, "a strong pull should hold the spread in"


def test_a_jump_at_the_split_is_logged_as_such():
    r = _run(seed=5, at_speciation=0.5)
    jumps = [c for c in r.trait.events if c.kind == "on_speciation"]
    assert len(jumps) > 10
    assert all(c.from_state != c.to_state for c in jumps)


def test_a_skyline_on_the_variance_rate_is_taken():
    """σ² may read the clock, because a growing tree has one."""
    from zombi2.params import PerLineage as PL
    r = joint.simulate(
        species.birth_death(
            birth=PerLineage(0.4).scaled_by("trait", Curve(lambda x: math.exp(0.3 * x)), step=0.05),
            death=0.05, n_extant=40),
        traits.continuous(rate=PL(2.0).changing_at({3.0: 0.05})), seed=6)
    assert sum(1 for n in r.complete_tree.nodes.values() if n.fate == "extant") == 40


# --- what is refused ------------------------------------------------------------------------------

def test_a_diffusing_driver_needs_a_step():
    with pytest.raises(ValueError, match="needs a step"):
        joint.simulate(
            species.birth_death(birth=PerLineage(0.4).scaled_by("trait", Curve(math.exp)),
                                n_extant=20),
            traits.continuous(rate=1.0), seed=1)


def test_two_readings_of_one_live_trait_must_agree_on_the_step():
    with pytest.raises(ValueError, match="agree on step"):
        joint.simulate(
            species.birth_death(
                birth=PerLineage(0.4).scaled_by("trait", Curve(math.exp), step=0.05),
                death=PerLineage(0.1).scaled_by("trait", Curve(math.exp), step=0.01),
                n_extant=20),
            traits.continuous(rate=1.0), seed=1)


def test_a_state_table_cannot_map_a_diffusing_value():
    with pytest.raises(ValueError, match="never equals a state name"):
        joint.simulate(
            species.birth_death(
                birth=PerLineage(0.4).scaled_by("trait", {"small": 1.0, "large": 2.0}, step=0.05),
                n_extant=20),
            traits.continuous(rate=1.0), seed=1)


def test_a_variance_rate_that_reads_a_finished_tree_is_refused():
    """σ² scaled by the standing diversity reads the whole lineages-through-time curve, which only a
    finished tree has. Refused rather than quietly read as 1.0 (SPEC §5)."""
    from zombi2.params import TotalDiversity
    with pytest.raises(ValueError, match="still growing"):
        joint.simulate(
            species.birth_death(
                birth=PerLineage(0.4).scaled_by("trait", Curve(math.exp), step=0.05), n_extant=20),
            traits.continuous(rate=PerLineage(1.0).scaled_by(TotalDiversity(cap=100))), seed=1)


def test_the_optimum_and_the_pull_come_together():
    with pytest.raises(ValueError, match="both or neither"):
        traits.continuous(rate=1.0, reverts_to=2.0)


def test_a_continuous_trait_with_a_genome_on_a_given_tree_says_it_is_not_built():
    from zombi2 import genomes
    tree = species.simulate_species_tree(birth=1.0, n_extant=10, seed=1).complete_tree
    with pytest.raises(NotImplementedError, match="CONTINUOUS trait"):
        joint.simulate(genomes.genome(duplication=0.1, loss=0.1),
                       traits.continuous(rate=1.0), tree=tree, seed=1)


def test_a_run_that_can_never_reach_its_target_says_so_rather_than_walking_for_ever():
    """A sliced run cannot use the other engines' "nothing is scheduled, so stop" test: a rate of
    zero now says nothing about the rate one slice later, because the driver is still moving."""
    import zombi2.joint as J
    old = J._MAX_SLICES
    J._MAX_SLICES = 500
    try:
        with pytest.raises(RuntimeError, match="without reaching n_extant"):
            joint.simulate(
                species.birth_death(
                    birth=PerLineage(0.0).scaled_by("trait", Curve(lambda x: 1.0), step=0.05),
                    n_extant=20),
                traits.continuous(rate=1.0), seed=1)
    finally:
        J._MAX_SLICES = old


def test_the_runaway_guard_still_holds():
    """A driven birth rate feeds the driver's own sample: the fast lineages leave more descendants,
    which are fast too. Under a steep enough curve that has no realistic end, and the guard raises
    rather than truncating — a tree cut off at a size is no longer the process asked for."""
    with pytest.raises(RuntimeError, match="still growing"):
        joint.simulate(
            species.birth_death(
                birth=PerLineage(1.0).scaled_by("trait", Curve(lambda x: math.exp(3.0 * x)),
                                                step=0.05),
                total_time=20.0),
            traits.continuous(rate=4.0), seed=1, max_lineages=400)
