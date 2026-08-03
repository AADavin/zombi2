"""Joint slice 2 — a discrete trait drives speciation, grown jointly (BiSSE/MuSSE).

The joint half of the one mechanism: `mod.DrivenBy("trait", mapping)` with the live level name
"trait" instead of a file, grown by `joint.simulate`. Covers the process spec, the result shape,
determinism, the state-dependent-diversification signal, MuSSE, full BiSSE (λ and μ), and validation.
"""

from collections import Counter

import pytest

from zombi2 import traits
from zombi2 import joint
from zombi2.joint import JointResult
from zombi2.rates import modifiers as mod
from zombi2.traits import DiscreteTrait, TraitsResult


def _bisse(birth_large=4.0, death=0.2, switch=0.15, n_extant=200, seed=1):
    return joint.simulate_joint(
        birth=1.0 * mod.DrivenBy("trait", {"small": 1.0, "large": birth_large}),
        death=death,
        trait=traits.discrete(states=["small", "large"], switch=switch),
        n_extant=n_extant, seed=seed,
    )


def _fraction_large(res):
    c = Counter(res.trait.values.values())
    return c["large"] / sum(c.values())


# --- the process spec -----------------------------------------------------------------------------

def test_discrete_spec_is_unexecuted_bundle():
    spec = traits.discrete(states=["a", "b"], switch=0.1, start="a")
    assert isinstance(spec, DiscreteTrait)
    assert spec.states == ("a", "b") and spec.start == "a"


def test_discrete_spec_validates():
    with pytest.raises(ValueError):
        traits.discrete(states=["only"], switch=0.1)          # < 2 states
    with pytest.raises(ValueError):
        traits.discrete(states=["a", "a"], switch=0.1)        # duplicate
    with pytest.raises(ValueError):
        traits.discrete(states=["a", "b"])                    # no switch


# --- the result shape -----------------------------------------------------------------------------

def test_joint_result_carries_both_levels():
    res = _bisse(n_extant=120, seed=3)
    assert isinstance(res, JointResult)
    assert res.n_extant == 120
    assert isinstance(res.trait, TraitsResult) and res.trait.kind == "discrete"
    # trait state recorded at EVERY node (extant, extinct, internal), tips readable
    assert set(res.trait.node_values) == set(res.complete_tree.nodes)
    assert set(res.trait.values_by_id) == {n.id for n in res.complete_tree.extant_leaves()}
    # the derived stochastic map reconstructs from the switch log (durations sum to branch lengths)
    hist = res.trait.history
    for i, node in res.complete_tree.nodes.items():
        assert sum(d for _, d in hist[i]) == pytest.approx(node.end_time - node.birth_time)


def test_joint_writes_both_levels(tmp_path):
    res = _bisse(n_extant=80, seed=2)
    res.write(tmp_path)
    for f in ("species_complete.nwk", "species_extant.nwk", "species_events.tsv",
              "trait_values.tsv", "trait_events.tsv", "trait_tree.nwk"):
        assert (tmp_path / f).exists(), f"missing {f}"


# --- determinism ----------------------------------------------------------------------------------

def test_joint_is_deterministic():
    a, b = _bisse(seed=7), _bisse(seed=7)
    assert [(e.time, e.kind, e.node) for e in a.events] == [(e.time, e.kind, e.node) for e in b.events]
    assert [(c.time, c.kind, c.lineage, c.to_state) for c in a.trait.events] == \
           [(c.time, c.kind, c.lineage, c.to_state) for c in b.trait.events]


# --- the state-dependent-diversification signal ---------------------------------------------------

def test_fast_state_over_represented_at_tips():
    # "large" speciates 4× faster → it should dominate the tips, far above the 0.5 stationary split
    fracs = [_fraction_large(_bisse(seed=s)) for s in (1, 2, 3)]
    assert min(fracs) > 0.6, f"fast state not over-represented: {fracs}"


def test_asymmetry_beats_symmetry():
    # with equal birth the tip split has no diversification bias; asymmetry pushes it toward "large"
    asym = sum(_fraction_large(_bisse(birth_large=4.0, seed=s)) for s in (1, 2, 3)) / 3
    symm = sum(_fraction_large(_bisse(birth_large=1.0, seed=s)) for s in (1, 2, 3)) / 3
    assert asym > symm + 0.15, f"asymmetric {asym:.2f} not clearly above symmetric {symm:.2f}"


def test_musse_three_states_runs():
    res = joint.simulate_joint(
        birth=1.0 * mod.DrivenBy("trait", {"lo": 1.0, "mid": 2.0, "hi": 4.0}),
        death=0.1,
        trait=traits.discrete(states=["lo", "mid", "hi"], switch=0.2),
        n_extant=150, seed=4,
    )
    seen = set(res.trait.values.values())
    assert seen and seen <= {"lo", "mid", "hi"}


def test_full_bisse_drives_birth_and_death():
    # both λ and μ state-dependent: "large" speciates faster AND goes extinct slower
    res = joint.simulate_joint(
        birth=1.0 * mod.DrivenBy("trait", {"small": 1.0, "large": 3.0}),
        death=0.3 * mod.DrivenBy("trait", {"small": 2.0, "large": 0.5}),
        trait=traits.discrete(states=["small", "large"], switch=0.2),
        n_extant=150, seed=5,
    )
    assert _fraction_large(res) > 0.6


def test_total_time_mode():
    res = joint.simulate_joint(
        birth=1.0 * mod.DrivenBy("trait", {"small": 1.0, "large": 2.0}),
        death=0.1,
        trait=traits.discrete(states=["small", "large"], switch=0.3),
        total_time=4.0, seed=6,
    )
    # every extant lineage reaches the present at total_time
    assert all(n.end_time == pytest.approx(4.0) for n in res.complete_tree.extant_leaves())


# --- validation -----------------------------------------------------------------------------------

def test_trait_must_be_a_process_spec():
    with pytest.raises(TypeError, match="traits.discrete"):        # a dict is not a DiscreteTrait spec
        joint.simulate_joint(birth=1.0 * mod.DrivenBy("trait", {"a": 1.0}),
                       trait={"states": ["a", "b"]}, n_extant=10, seed=1)


def test_non_trait_source_rejected():
    with pytest.raises(ValueError, match="trait"):
        joint.simulate_joint(birth=1.0 * mod.DrivenBy("habitat.tsv", {"a": 1.0}),
                       trait=traits.discrete(states=["a", "b"], switch=0.1), n_extant=10, seed=1)


def test_must_actually_drive_something():
    with pytest.raises(ValueError, match="drive"):
        joint.simulate_joint(birth=1.0, death=0.1,
                       trait=traits.discrete(states=["a", "b"], switch=0.1), n_extant=10, seed=1)


def test_one_of_n_extant_or_total_time():
    with pytest.raises(ValueError, match="exactly one"):
        joint.simulate_joint(birth=1.0 * mod.DrivenBy("trait", {"a": 2.0}),
                       trait=traits.discrete(states=["a", "b"], switch=0.1),
                       n_extant=10, total_time=3.0, seed=1)


def test_fromparent_rejected():
    with pytest.raises(ValueError, match="FromParent"):
        joint.simulate_joint(birth=1.0 * mod.FromParent(spread=0.2) * mod.DrivenBy("trait", {"a": 2.0}),
                       trait=traits.discrete(states=["a", "b"], switch=0.1), n_extant=10, seed=1)


def test_mapping_matching_no_trait_state_is_refused():
    # a driving mapping whose states are none of the trait's would leave every lineage at the default
    # factor — a silently undriven run — so it is refused rather than run as if it were driven. The
    # trait's alphabet is known up front, so the exhaustive check names the offending key precisely.
    with pytest.raises(ValueError, match="not among the driver's states"):
        joint.simulate_joint(birth=1.0 * mod.DrivenBy("trait", {"tiny": 2.0}),   # trait is small/large
                       death=0.2, trait=traits.discrete(states=["small", "large"], switch=0.15),
                       n_extant=10, seed=1)


def test_mapping_naming_a_non_trait_state_is_refused():
    # the trait's declared states are known up front, so the joint check is exhaustive: even when some
    # keys match, a key outside the alphabet ('tiny' ∉ small/large) is a state that can never occur — a
    # typo whose factor would silently never apply — so the whole mapping is refused, not run as driven
    with pytest.raises(ValueError, match="not among the driver's states"):
        joint.simulate_joint(birth=1.0 * mod.DrivenBy("trait", {"small": 2.0, "tiny": 3.0}),
                       death=0.2, trait=traits.discrete(states=["small", "large"], switch=0.15),
                       n_extant=10, seed=1)


def test_a_gene_content_mapping_that_can_never_fire_is_refused_too():
    """The joint GENE-CONTENT path gets the same guard the joint trait path has.

    It did not have one. A typo'd ``{"presnt": 3.0}`` left every lineage at the default factor, so
    birth was never actually driven — the run was the plain birth-death model — and it completed in
    silence, reporting a coupled run. The failure is the worst kind for a simulator: the output is a
    well-formed tree and a well-formed genome, and only the association the run was made to measure
    is missing.

    Both gene-content drivers have an alphabet known before the race starts, so the check is
    exhaustive. A named family is ``present`` or ``absent``; a count is a number, which a
    ``{state: factor}`` table can never equal."""
    from zombi2 import genomes as g

    def run(mapping, **kw):
        return joint.simulate_joint(
            birth=1.0 * mod.DrivenBy(kw.pop("source"), mapping),
            genome=g.family(duplication=0.1, loss=0.1, origination=0.2, initial_families=3,
                            family_names=["toxin"]),
            n_extant=15, seed=1)

    with pytest.raises(ValueError, match="not among the driver's states"):
        run({"presnt": 3.0, "absent": 1.0}, source="genomes:toxin")      # typo: presnt
    with pytest.raises(ValueError, match="CONTINUOUS"):
        run({"many": 3.0}, source="genomes:count")                       # a table on a numeric driver

    result = run({"present": 3.0, "absent": 1.0}, source="genomes:toxin")  # the correct spelling runs
    assert len(result.complete_tree.nodes) > 15


# --- a trait that drives speciation and also jumps at the split ----------------------------------

def test_joint_trait_can_jump_at_speciation_while_driving_it():
    """The two options compose: `at_speciation=` gives the trait a chance of changing *at* a split,
    and the same trait drives how fast that split happens. Both kinds of change must appear in the
    log, distinguishable, or the model is only one of its two halves."""
    import collections
    from zombi2 import joint, traits
    from zombi2.rates import modifiers as mod

    res = joint.simulate_joint(
        birth=1.0 * mod.DrivenBy("trait", {"small": 1.0, "large": 2.5}),
        death=0.2,
        trait=traits.discrete(states=["small", "large"], switch=0.3, at_speciation=0.4),
        n_extant=30, seed=1)

    kinds = collections.Counter(e.kind for e in res.trait.events)
    assert kinds["on_speciation"] > 0, "no jump at a split — at_speciation was ignored"
    assert kinds["on_branch"] > 0, "no change along a branch — the CTMC was ignored"


def test_joint_speciation_jump_is_off_by_default():
    """Without `at_speciation`, a joint trait changes only along branches — the jump is opt-in, and
    a run that does not ask for it must not get it."""
    import collections
    from zombi2 import joint, traits
    from zombi2.rates import modifiers as mod

    res = joint.simulate_joint(
        birth=1.0 * mod.DrivenBy("trait", {"small": 1.0, "large": 2.5}),
        death=0.2,
        trait=traits.discrete(states=["small", "large"], switch=0.3),
        n_extant=30, seed=1)
    assert collections.Counter(e.kind for e in res.trait.events)["on_speciation"] == 0


def test_joint_refuses_per_lineage_rate_variation():
    """The species level takes ByLineage; the joint engine does not thread it.

    Accepting it here would run a model without the rate variation that was asked for — and because
    the same ``--birth`` expression works on ``zombi2 species``, silently ignoring it is a trap rather
    than merely a gap. It used to be accepted (SPEC §5's rejection rule was applied to FromParent and
    DrivenBy but not to this one)."""
    import pytest

    from zombi2 import joint, traits
    from zombi2.rates import modifiers as mod

    with pytest.raises(ValueError, match="ByLineage"):
        joint.simulate_joint(
            birth=1.0 * mod.ByLineage(spread=0.5) * mod.DrivenBy("trait", {"small": 1.0, "large": 2.0}),
            death=0.1, n_extant=8, seed=1,
            trait=traits.DiscreteTrait(states=("small", "large"), switch=0.3))


def test_joint_refuses_a_modifier_it_does_not_thread():
    """The gate was a negative list — it named FromParent and DrivenBy and let everything else
    through — where every other level declares what it takes. `ByFamily` was the one that slipped:
    accepted, then returning its default factor of 1.0, so the run was quietly not the model asked
    for. `OnTime` and `OnTotalDiversity` were never the problem; the loop threads both and steps at
    their breakpoints, which is why they are in `WIRED_MODIFIERS` rather than refused alongside."""
    import pytest

    from zombi2 import joint, traits
    from zombi2.rates import modifiers as mod

    trait = traits.DiscreteTrait(states=("small", "large"), switch=0.3)
    driven = mod.DrivenBy("trait", {"small": 1.0, "large": 2.0})
    with pytest.raises(ValueError, match="no gene families"):
        joint.simulate_joint(birth=1.0 * mod.ByFamily(spread=0.5) * driven, death=0.1,
                             n_extant=8, seed=1, trait=trait)
    # the two covariates are genuinely wired, so they must still run
    for modifier in (mod.OnTime({0: 1.0, 0.2: 0.4}), mod.OnTotalDiversity(cap=40)):
        assert joint.simulate_joint(birth=1.0 * modifier * driven, death=0.1, n_extant=8,
                                    seed=1, trait=trait).species.n_extant == 8
