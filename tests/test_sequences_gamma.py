"""Across-site rate variation: the discretised Gamma (``+Γ``) and the invariant class (``+I``).

Two things are being checked, and they pull in opposite directions. The **classes must be right** —
they are the only new maths in the feature, and a mis-discretised Gamma looks perfectly plausible
from the outside — so they are pinned against the values published for PAML, against scipy, and
against the mean-1 contract directly. And the **rest of the level must not move**: a model without
rate variation has to draw exactly what it drew before (there is a golden pin for that), and a
phylogram's branch length has to stay the mean substitutions per site, which is what lets a run with
variation and one without be compared at all.
"""

from __future__ import annotations

import numpy as np
import pytest

from zombi2.rates import LogNormal, PerSite
from zombi2.sequences import simulate_sequences
from zombi2.sequences._site_rates import _gamma_quantile, _incomplete_gamma, discrete_gamma
from zombi2.sequences.substitution_models import hky85, jc69

from test_sequences import _pair_run     # the hand-built gene tree the level's own tests use


# --- the discretisation itself --------------------------------------------------------------------

def test_discrete_gamma_matches_the_published_categories():
    """The standard four-category values (Yang 1994), as printed by PAML and every tool since.

    Pinned rather than derived: this is the one place the implementation could be subtly wrong in a
    way nothing downstream would reveal, since any monotone set of classes produces plausible-looking
    alignments."""
    assert discrete_gamma(0.5, 4) == pytest.approx((0.033388, 0.251916, 0.820268, 2.894428), abs=1e-5)
    assert discrete_gamma(1.0, 4) == pytest.approx((0.136954, 0.476752, 1.0, 2.386294), abs=1e-5)


def test_the_hand_rolled_numerics_match_scipy():
    """scipy is a test-only dependency here, so the incomplete gamma and its inverse are written out
    in `_site_rates`. This is the check that says the hand-rolled pair is the real thing."""
    special = pytest.importorskip("scipy.special")
    stats = pytest.importorskip("scipy.stats")
    for a in (0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0):
        for x in (1e-4, 0.01, 0.5, 1.0, 3.0, 10.0, 50.0, 200.0):
            assert _incomplete_gamma(a, x) == pytest.approx(float(special.gammainc(a, x)), abs=1e-12)
        for p in (0.001, 0.1, 0.25, 0.5, 0.75, 0.9, 0.999):
            assert _gamma_quantile(p, a) == pytest.approx(float(stats.gamma.ppf(p, a)), rel=1e-9)


@pytest.mark.parametrize("shape", [0.05, 0.5, 1.0, 5.0, 100.0])
@pytest.mark.parametrize("categories", [2, 4, 8, 16])
@pytest.mark.parametrize("invariant", [0.0, 0.3, 0.7])
def test_every_class_set_has_mean_one(shape, categories, invariant):
    # the invariant the phylograms rest on, asserted directly rather than inferred from an alignment
    m = jc69().across_sites(gamma_shape=shape, invariant=invariant, rate_categories=categories)
    assert sum(r * s for r, s in zip(m.site_rates, m.site_shares)) == pytest.approx(1.0, abs=1e-12)


def test_invariant_alone_still_averages_one():
    # +I without a Gamma is one never-changing class and one ordinary one; the ordinary sites have to
    # speed up to carry the whole mean, or +I would quietly slow the sequence down instead of
    # concentrating its change
    m = jc69().across_sites(invariant=0.25)
    assert m.site_rates == pytest.approx((0.0, 1 / 0.75))
    assert m.site_shares == pytest.approx((0.25, 0.75))


def test_the_model_name_records_the_decoration():
    # the name is what the run's summary line prints and what the log records, so it is contract
    assert hky85().across_sites(gamma_shape=0.5).name == "HKY85+G4"
    assert hky85().across_sites(gamma_shape=0.5, invariant=0.1).name == "HKY85+I+G4"
    assert hky85().across_sites(invariant=0.1).name == "HKY85+I"
    assert hky85().across_sites(gamma_shape=0.5, rate_categories=8).name == "HKY85+G8"


def test_decorating_twice_is_refused_rather_than_silently_replacing():
    """A second call replaces the classes instead of layering on them, and the name would claim both
    — `JC69+G4+I` over a model that is only `+I`. Every combination has one spelling: one call."""
    once = jc69().across_sites(gamma_shape=0.5)
    with pytest.raises(ValueError, match="already varies across sites"):
        once.across_sites(invariant=0.2)


@pytest.mark.parametrize("shape", [1e-3, 0.01, 1e3, 1e6])
def test_extreme_shapes_stay_normalised(shape):
    """The quantile is found by bracket-and-bisect, so the far tails are where it would give up.

    Only the contract is asserted here — non-negative classes averaging 1. Ascending order is
    checked separately, below: at a shape this small the classes are all within round-off of zero
    except the last, and which of several 1e-60 values comes out larger is noise rather than
    meaning."""
    rates = jc69().across_sites(gamma_shape=shape).site_rates
    assert sum(rates) / len(rates) == pytest.approx(1.0, abs=1e-12)
    assert all(r >= 0.0 for r in rates)


@pytest.mark.parametrize("shape", [0.1, 0.5, 1.0, 2.0, 5.0, 20.0])
def test_the_classes_come_out_in_ascending_order(shape):
    # over the range anyone actually fits (published analyses live between about 0.1 and 5)
    rates = jc69().across_sites(gamma_shape=shape, rate_categories=8).site_rates
    assert list(rates) == sorted(rates)


def test_the_original_model_is_left_alone():
    m = hky85(2.0)
    m.across_sites(gamma_shape=0.5)
    assert m.name == "HKY85" and m.site_rates == (1.0,)      # frozen: decorating returns a new model


@pytest.mark.parametrize("kwargs, message", [
    (dict(gamma_shape=0.0), "positive"),
    (dict(gamma_shape=-1.0), "positive"),
    (dict(gamma_shape=float("inf")), "finite"),
    (dict(invariant=1.0), r"\[0, 1\)"),
    (dict(invariant=-0.1), r"\[0, 1\)"),
    (dict(gamma_shape=0.5, rate_categories=1), "at least 2"),
    (dict(rate_categories=6), "needs.*gamma_shape|add gamma_shape"),
    (dict(), "nothing to vary"),
])
def test_across_sites_rejects_what_it_cannot_mean(kwargs, message):
    with pytest.raises((ValueError, TypeError), match=message):
        jc69().across_sites(**kwargs)


# --- what it does to a run -------------------------------------------------------------------------

def test_a_flat_model_still_draws_exactly_what_it_drew():
    """A golden pin, because the whole design rests on one guard: a model with a single rate class
    must not consume a single extra draw, so adding across-site variation leaves a flat run
    bit-identical. Re-pinned when each level moved to its own RNG stream (`zombi2.rng`) — that
    deliberately changed what a sequences seed draws, which is a different thing from this guard,
    and the guard still holds against the run beside it."""
    r = simulate_sequences(_pair_run(1.0, 2.0), model=jc69(), length=40, seed=7)
    assert r.alignments[0]["n1_g1"] == "GCAGATGCGGAGGTGGATCTTGTCCATTCAATTTACCGGG"
    assert r.founding[0] == "GCAGGAAAGGATGAAAATACTGGAGACTCGATAATGAACA"


def test_the_phylogram_is_unchanged_by_across_site_variation():
    """Branch lengths stay the **mean** substitutions per site.

    This is the statement that makes rate variation comparable to none: the classes are normalised to
    mean 1, so the tree behind a ``+I+G4`` alignment is the same tree as behind a flat one at the same
    rate, and only the spread of change across columns differs."""
    run = _pair_run(1.0, 2.0)
    flat = simulate_sequences(run, model=jc69(), length=200, substitution=0.4, seed=5)
    varied = simulate_sequences(run, model=jc69().across_sites(gamma_shape=0.4, invariant=0.2),
                                length=200, substitution=0.4, seed=5)
    assert flat.phylograms == varied.phylograms
    assert flat.species_phylogram == varied.species_phylogram


def test_the_mean_rate_is_unchanged_by_rate_variation():
    """Measured, not assumed: at a divergence low enough that a site is almost never hit twice, the
    observed p-distance from the founding sequence *is* the mean rate — so a mis-normalised class set
    shows up here as a number that has moved. A structural check on the classes would not catch a
    normalisation applied in the wrong place."""
    run = _pair_run(0.0, 1.0)

    def p_distance(model, seed):
        r = simulate_sequences(run, model=model, length=200_000, substitution=0.002, seed=seed)
        founding = np.frombuffer(r.founding[0].encode(), dtype=np.uint8)
        return np.mean([np.mean(np.frombuffer(s.encode(), dtype=np.uint8) != founding)
                        for s in r.alignments[0].values()])

    flat = np.mean([p_distance(jc69(), s) for s in (11, 12, 13)])
    for model in (jc69().across_sites(gamma_shape=0.3),
                  jc69().across_sites(gamma_shape=0.1, rate_categories=8),
                  jc69().across_sites(invariant=0.4),
                  jc69().across_sites(gamma_shape=0.5, invariant=0.25)):
        got = np.mean([p_distance(model, s) for s in (11, 12, 13)])
        assert got == pytest.approx(flat, rel=0.05), f"{model.name} shifted the mean rate"


def _constant_columns(result, family: int = 0) -> np.ndarray:
    """Per column, whether every sequence in the family agrees there."""
    rows = np.array([np.frombuffer(s.encode(), dtype=np.uint8)
                     for s in result.alignments[family].values()])
    return (rows == rows[0]).all(axis=0)


def test_gamma_leaves_more_columns_untouched():
    """The point of ``+Γ`` is the unevenness, so that is what is asserted.

    At one mean rate, concentrating the change into fewer sites leaves more columns completely
    unchanged than spreading it evenly does — which is exactly the pattern real alignments show and
    the reason the model exists."""
    run = _pair_run(0.0, 1.0)
    flat = simulate_sequences(run, model=jc69(), length=40_000, substitution=0.5, seed=4)
    varied = simulate_sequences(run, model=jc69().across_sites(gamma_shape=0.2),
                                length=40_000, substitution=0.5, seed=4)
    assert _constant_columns(varied).mean() > _constant_columns(flat).mean() + 0.1


def test_invariant_sites_never_change():
    """``invariant=0.4`` at a divergence that saturates everything else: at least 40% of columns must
    come through untouched, because those sites cannot change at all."""
    run = _pair_run(0.0, 1.0)
    saturated = simulate_sequences(run, model=jc69(), length=20_000, substitution=5.0, seed=6)
    with_inv = simulate_sequences(run, model=jc69().across_sites(invariant=0.4),
                                  length=20_000, substitution=5.0, seed=6)
    assert _constant_columns(with_inv).mean() >= 0.39      # 0.4 minus sampling tolerance
    assert _constant_columns(saturated).mean() < 0.35      # nothing is protected without +I


def test_the_invariant_sites_are_the_founding_letters():
    """A site in the rate-0 class keeps its founding state down the whole tree — not merely 'changes
    rarely'. Checked against the founding sequence rather than between tips, which is what says the
    class was skipped rather than sampled with a tiny probability."""
    run = _pair_run(0.0, 1.0)
    r = simulate_sequences(run, model=jc69().across_sites(invariant=0.5), length=20_000,
                           substitution=8.0, seed=8)
    founding = np.frombuffer(r.founding[0].encode(), dtype=np.uint8)
    matches_founding = np.array([np.frombuffer(s.encode(), dtype=np.uint8) == founding
                                 for s in r.alignments[0].values()]).all(axis=0)
    # at this divergence a variable site matches the founding letter about a quarter of the time by
    # chance, so ~0.5 protected + ~0.5*0.25 lucky ≈ 0.62; the floor below is the protected half
    assert matches_founding.mean() >= 0.49


def test_rate_variation_composes_with_the_lineage_clock():
    """The two axes are orthogonal — the clock says which lineages run fast, the model which sites do
    — and both apply at once, deterministically."""
    run = _pair_run(1.0, 2.0)
    model = jc69().across_sites(gamma_shape=0.5, invariant=0.1)
    kwargs = dict(model=model, length=300, substitution=PerSite(1.0).varying_among('lineages', LogNormal(0.0, 0.3)), seed=3)
    a = simulate_sequences(run, **kwargs)
    b = simulate_sequences(run, **kwargs)
    strict = simulate_sequences(run, model=model, length=300, seed=3)
    assert a.alignments == b.alignments                    # reproducible
    assert a.alignments != strict.alignments               # and the clock is doing something


def test_a_protein_model_takes_rate_variation_too():
    # the classes multiply a branch length, so they know nothing about the alphabet
    from zombi2.sequences.substitution_models import AMINO_ACIDS, lg

    m = lg().across_sites(gamma_shape=0.7)
    assert m.k == 20 and sum(r * s for r, s in zip(m.site_rates, m.site_shares)) == pytest.approx(1.0)
    r = simulate_sequences(_pair_run(1.0, 2.0), model=m, length=200, seed=2)
    assert set("".join(r.alignments[0].values())) <= set(AMINO_ACIDS)
