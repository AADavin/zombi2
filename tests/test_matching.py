"""Rejection-ABC profile matching (:func:`zombi2.match_profiles`).

The uniform (built-in) model that the plumbing and recovery tests exercise runs on the Rust
engine, so the whole module is skipped when ``zombi2_core`` isn't built. ``match_profiles``
picks the engine automatically (uniform -> Rust, family/callable -> Python); there is no
``engine`` argument.
"""

import numpy as np
import pytest

import zombi2 as z

pytestmark = pytest.mark.skipif(not z.rust_available(),
                                reason="zombi2_core (Rust extension) not built")


# --- summary statistics ----------------------------------------------------------

def _pm(matrix, species, families=None):
    matrix = np.array(matrix, dtype=int)
    families = families or [str(i + 1) for i in range(matrix.shape[0])]
    return z.ProfileMatrix(families=families, species=species, matrix=matrix)


def test_frequency_spectrum():
    pm = _pm([[1, 0, 0], [2, 1, 0], [1, 1, 3]], ["A", "B", "C"])
    # families present in 1, 2, 3 species respectively -> one each
    assert np.array_equal(z.frequency_spectrum(pm, 3), [1, 1, 1])


def test_frequency_spectrum_empty():
    pm = _pm(np.zeros((0, 4), dtype=int), ["A", "B", "C", "D"], families=[])
    assert np.array_equal(z.frequency_spectrum(pm, 4), [0, 0, 0, 0])


def test_genome_sizes_align_and_missing():
    pm = _pm([[1, 0, 0], [2, 1, 0], [1, 1, 3]], ["A", "B", "C"])
    # column sums are A=4, B=2, C=3; requested order reshuffles and adds an absent species
    assert np.array_equal(z.genome_sizes(pm, ["C", "A", "B"]), [3, 4, 2])
    assert np.array_equal(z.genome_sizes(pm, ["A", "D"]), [4, 0])


def test_copy_number_spectrum_lumps_tail():
    pm = _pm([[1, 0, 0], [2, 1, 0], [1, 1, 5]], ["A", "B", "C"])
    # present values: 1,2,1,1,1,5 -> copy 1:four, 2:one, 3:none, >=4:one
    assert np.array_equal(z.copy_number_spectrum(pm, max_copies=4), [4, 1, 0, 1])


def test_default_summary_length():
    pm = _pm([[1, 0, 0], [2, 1, 0]], ["A", "B", "C"])
    summarize = z.default_summary(["A", "B", "C"])
    # frequency spectrum (S) + genome sizes (S) + copy spectrum (4)
    assert summarize(pm).shape == (3 + 3 + 4,)


# --- TSV round-trip --------------------------------------------------------------

def test_tsv_round_trip_text_and_file(tmp_path):
    pm = _pm([[1, 0, 2], [0, 3, 1]], ["sp1", "sp2", "sp3"], families=["famA", "famB"])
    back = z.ProfileMatrix.from_tsv(pm.to_tsv())
    assert back.species == pm.species and back.families == pm.families
    assert np.array_equal(back.matrix, pm.matrix)

    path = tmp_path / "profiles.tsv"
    path.write_text(pm.to_tsv())
    from_file = z.ProfileMatrix.from_tsv(str(path))
    assert np.array_equal(from_file.matrix, pm.matrix)


# --- prior handling & validation -------------------------------------------------

def _small_tree(n=6, seed=1):
    return z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=n, age=4.0, seed=seed)


def _cheap_fit(**kw):
    """A tiny fit for plumbing checks (uniform model -> Rust engine)."""
    tree = _small_tree()
    emp = z.simulate_genomes(tree, duplication=0.1, loss=0.15, origination=0.6,
                             initial_families=10, seed=3).profiles
    params = dict(tree=tree, empirical=emp,
                  priors={"duplication": (0, 0.3), "loss": (0, 0.4), "origination": (0, 1.5)},
                  n_sims=15, accept=0.2, initial_families=10, seed=1)
    params.update(kw)
    return z.match_profiles(**params)


def test_unknown_parameter_rejected():
    tree = _small_tree()
    emp = z.simulate_genomes(tree, duplication=0.1, origination=0.5, seed=1).profiles
    with pytest.raises(ValueError):
        z.match_profiles(tree, emp, priors={"speciation": (0, 1)}, n_sims=4)


def test_empty_priors_rejected():
    tree = _small_tree()
    emp = z.simulate_genomes(tree, duplication=0.1, origination=0.5, seed=1).profiles
    with pytest.raises(ValueError):
        z.match_profiles(tree, emp, priors={}, n_sims=4)


def test_fixed_float_prior_is_constant():
    # a bare float prior is a fixed value -> that column never varies
    fit = _cheap_fit(priors={"duplication": 0.05, "loss": (0, 0.4), "origination": (0, 1.5)})
    dup = fit.samples[:, fit.param_names.index("duplication")]
    assert np.allclose(dup, 0.05)


def test_bad_accept():
    tree = _small_tree()
    emp = z.simulate_genomes(tree, duplication=0.1, origination=0.5, seed=1).profiles
    priors = {"duplication": (0, 0.3), "origination": (0, 1.5)}
    with pytest.raises(ValueError):
        z.match_profiles(tree, emp, priors=priors, accept=1.5, n_sims=4)


# --- ABCFit plumbing -------------------------------------------------------------

def test_accepted_are_the_closest_and_shapes():
    fit = _cheap_fit()  # n_sims=15, accept=0.2 -> k=3
    k = len(fit.accepted)
    assert k == round(0.2 * 15)
    assert fit.samples.shape == (15, len(fit.param_names))
    assert fit.distances.shape == (15,)
    # accepted are exactly the k smallest distances, and tolerance is their max
    expected = set(np.argsort(fit.distances)[:k].tolist())
    assert set(fit.accepted.tolist()) == expected
    assert fit.tolerance == pytest.approx(fit.distances[fit.accepted].max())
    assert (fit.distances[fit.accepted] <= fit.tolerance + 1e-12).all()


def test_best_is_argmin_and_posterior_keys():
    fit = _cheap_fit()
    assert fit.best == {n: float(fit.samples[np.argmin(fit.distances), j])
                        for j, n in enumerate(fit.param_names)}
    assert set(fit.posterior) == set(fit.param_names)
    assert all(v.shape == (len(fit.accepted),) for v in fit.posterior.values())


def test_accept_as_integer_count():
    fit = _cheap_fit(accept=4)
    assert len(fit.accepted) == 4


def test_determinism_same_seed():
    a, b = _cheap_fit(), _cheap_fit()
    assert np.array_equal(a.samples, b.samples)
    assert np.array_equal(a.distances, b.distances)
    assert np.array_equal(a.accepted, b.accepted)


# --- end-to-end recovery ---------------------------------------------------------

def test_recover_injected_rates_43_leaves():
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=43, age=5.0, seed=1)
    truth = dict(duplication=0.3, transfer=0.1, loss=0.6, origination=2.0)
    emp = z.simulate_genomes(tree, initial_families=20, seed=101, output="profiles", **truth)

    priors = {"duplication": (0, 1), "transfer": (0, 0.5),
              "loss": (0, 1.5), "origination": (0, 3)}
    fit = z.match_profiles(tree, emp, priors=priors, n_sims=800, accept=0.05,
                           initial_families=20, seed=7)

    s = fit.summary()
    # The well-identified rates: truth inside the 95% credible interval...
    assert s["origination"]["lo95"] <= 2.0 <= s["origination"]["hi95"]
    assert s["duplication"]["lo95"] <= 0.3 <= s["duplication"]["hi95"]
    # ...and the posterior median lands near the truth.
    assert abs(s["origination"]["median"] - 2.0) <= 0.8   # within ~40%
    assert abs(s["duplication"]["median"] - 0.3) <= 0.15

    # The data concentrate the posterior below the prior for the identifiable rates.
    prior_std = {"origination": 3 / np.sqrt(12), "duplication": 1 / np.sqrt(12),
                 "loss": 1.5 / np.sqrt(12)}
    post = fit.posterior
    for name, ps in prior_std.items():
        assert post[name].std() < ps

    # The search found a materially better fit than a typical prior draw.
    assert fit.distances.min() < 0.8 * np.median(fit.distances)

    # NOTE: loss and transfer are NOT asserted as point-recovered. From copy-number
    # profiles alone the gain/loss balance is a ridge (raising origination pulls loss up
    # to compensate), and transfer barely moves presence structure — exactly the
    # identifiability limit match_profiles is meant to expose, so we report a posterior,
    # not a point.


# --- parallel inner loop ---------------------------------------------------------

def test_parallel_matches_serial():
    # results depend only on the seed, not the process count (draws are pre-generated)
    tree = _small_tree()
    emp = z.simulate_genomes(tree, duplication=0.1, loss=0.15, origination=0.6,
                             initial_families=10, seed=3).profiles
    kw = dict(priors={"duplication": (0, 0.3), "loss": (0, 0.4), "origination": (0, 1.5)},
              n_sims=40, accept=0.2, initial_families=10, seed=1)
    serial = z.match_profiles(tree, emp, **kw)
    parallel = z.match_profiles(tree, emp, processes=2, **kw)
    assert np.array_equal(serial.samples, parallel.samples)
    assert np.allclose(serial.distances, parallel.distances)
    assert np.array_equal(serial.accepted, parallel.accepted)


# --- family-sampled-rates model --------------------------------------------------

def test_custom_model_allows_arbitrary_param_names():
    # a callable model params->RateModel may use any parameter names (runs on Python engine)
    tree = _small_tree()
    emp = z.simulate_genomes(tree, duplication=0.1, loss=0.15, origination=0.6,
                             initial_families=10, seed=3).profiles
    model = lambda p: z.UniformRates(duplication=p["d"], loss=p["l"], origination=p["o"])
    fit = z.match_profiles(tree, emp, priors={"d": (0, 0.3), "l": (0, 0.4), "o": (0, 1.5)},
                           model=model, n_sims=20, accept=0.2, initial_families=10, seed=1)
    assert set(fit.posterior) == {"d", "l", "o"}


def test_family_model_recovers_heterogeneous_rates():
    # per-family rates drawn around known means; a hard cap keeps Python sims bounded
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=20, age=5.0, seed=1)
    truth = z.FamilySampledRates(duplication=z.Gamma(2, 0.15), transfer=z.Gamma(2, 0.05),
                                 loss=z.Gamma(2, 0.30), origination=2.0)
    emp = z.simulate_genomes(tree, truth, initial_families=15, max_family_size=20, seed=101).profiles

    priors = {"duplication": (0, 0.6), "transfer": (0, 0.3),
              "loss": (0, 1.2), "origination": (0, 3)}
    fit = z.match_profiles(tree, emp, priors=priors, model="family", family_shape=2.0,
                           n_sims=150, accept=0.15, initial_families=15, max_family_size=20,
                           seed=7, processes=2)
    assert set(fit.posterior) == set(z.matching.RATE_PARAMS)
    s = fit.summary()
    # the identifiable means (duplication from copy number, origination from family count)
    assert s["duplication"]["lo95"] <= 0.30 <= s["duplication"]["hi95"]
    assert s["origination"]["lo95"] <= 2.0 <= s["origination"]["hi95"]


# --- spectrum diagnostic ---------------------------------------------------------

def test_spectra_data_shapes_and_custom_guard():
    fit = _cheap_fit()
    d = fit.spectra_data()
    s = fit.n_species
    assert len(d["k"]) == s and len(d["empirical"]) == s
    assert d["accepted"].shape == (len(fit.accepted), s)

    # a custom summary has no known spectrum slice -> the diagnostic refuses
    tree = _small_tree()
    emp = z.simulate_genomes(tree, duplication=0.1, origination=0.5, seed=1).profiles
    fit2 = z.match_profiles(tree, emp, priors={"origination": (0, 1)},
                            statistics=lambda pm: np.array([float(pm.matrix.sum())]),
                            n_sims=10, accept=0.3, seed=1)
    with pytest.raises(ValueError):
        fit2.spectra_data()


def _small_genomes(seed=3):
    tree = _small_tree()
    g = z.simulate_genomes(tree, duplication=0.12, transfer=0.05, loss=0.18, origination=0.7,
                           initial_families=10, seed=seed)
    return tree, g


def test_event_count_summary_matches_log():
    from collections import Counter
    _, g = _small_genomes()
    ev = z.event_count_summary(g)
    c = Counter(r.event for r in g.event_log)
    assert ev.shape == (3,)
    assert ev[0] == c[z.EventType.DUPLICATION]
    assert ev[1] == c[z.EventType.TRANSFER]
    assert ev[2] == c[z.EventType.LOSS]


def test_gene_tree_summary_layout_and_weights():
    tree, g = _small_genomes()
    species = [n.name for n in tree.extant_leaves()]
    s = len(species)
    summ = z.default_gene_tree_summary(species)
    vec = summ(g)
    assert vec.shape == (2 * s + 4 + 3,)
    # the frequency spectrum still leads the vector (so the spectrum diagnostic keeps working)
    assert np.array_equal(vec[:s], z.frequency_spectrum(g.profiles, s))
    # weights: 1 for the profile block, one equal value > 1 for the three event counts
    w = summ.feature_weights
    assert (w[:2 * s + 4] == 1).all()
    assert w[-1] > 1 and len(set(np.round(w[-3:], 6))) == 1


def test_gene_trees_path_runs_with_diagnostics():
    tree, g = _small_genomes()
    fit = z.match_profiles(tree, g, priors={"duplication": (0, 0.4), "loss": (0, 0.5),
                                            "origination": (0, 1.5)},
                           gene_trees=True, n_sims=40, accept=0.25, initial_families=10, seed=1)
    assert set(fit.posterior) == {"duplication", "loss", "origination"}
    assert fit.spectra_data()["accepted"].shape[1] == fit.n_species


def test_gene_trees_requires_genomes():
    tree, g = _small_genomes()
    pm = g.profiles
    with pytest.raises(TypeError):        # a bare profile lacks gene trees
        z.match_profiles(tree, pm, priors={"origination": (0, 1)}, gene_trees=True, n_sims=4)


def test_feature_weights_change_acceptance():
    tree = _small_tree()
    emp = z.simulate_genomes(tree, duplication=0.1, loss=0.15, origination=0.6,
                             initial_families=10, seed=3).profiles
    kw = dict(priors={"duplication": (0, 0.3), "loss": (0, 0.4), "origination": (0, 1.5)},
              n_sims=60, accept=0.25, initial_families=10, seed=1)
    base = z.match_profiles(tree, emp, **kw)
    s = len(emp.species)
    w = np.ones(2 * s + 4)
    w[:s] = 5.0                            # up-weight the frequency spectrum
    weighted = z.match_profiles(tree, emp, feature_weights=w, **kw)
    # same draws, different distance -> generally a different accepted set
    assert not np.array_equal(base.accepted, weighted.accepted)


def _adjust_fit():
    tree = _small_tree()
    emp = z.simulate_genomes(tree, duplication=0.12, loss=0.18, origination=0.7,
                             initial_families=10, seed=3).profiles
    return z.match_profiles(tree, emp,
                            priors={"duplication": (0, 0.4), "loss": (0, 0.5), "origination": (0, 2)},
                            n_sims=120, accept=0.25, initial_families=10, seed=1)


def test_regression_adjust_shapes_and_nonneg():
    fit = _adjust_fit()
    adj = fit.regression_adjust()
    assert set(adj) == set(fit.param_names)
    k = len(fit.accepted)
    for v in adj.values():
        assert v.shape == (k,) and (v >= 0).all()
    assert fit.regression_adjust() is adj  # cached


def test_adjusted_summary_is_ordered():
    fit = _adjust_fit()
    s = fit.summary(adjusted=True)
    assert set(s) == set(fit.param_names)
    for st in s.values():
        assert set(st) == {"mean", "median", "lo95", "hi95"}
        assert st["lo95"] <= st["median"] <= st["hi95"]


def test_weighted_quantile_matches_numpy_when_uniform():
    v = np.array([3.0, 1.0, 2.0, 5.0, 4.0])
    w = np.ones_like(v)
    assert z.matching._wquantile(v, w, 0.5) == pytest.approx(np.median(v))


def test_plot_spectra_returns_axes():
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fit = _cheap_fit()
    ax = fit.plot_spectra(draws=True)
    assert ax.get_ylabel() == "number of gene families"
    plt.close("all")


# --- ABC-SMC ---------------------------------------------------------------------

def _smc_setup():
    tree = _small_tree()
    emp = z.simulate_genomes(tree, duplication=0.1, loss=0.15, origination=0.6,
                             initial_families=10, seed=3).profiles
    return tree, emp


def test_smc_runs_and_is_weighted():
    tree, emp = _smc_setup()
    fit = z.match_profiles_smc(
        tree, emp, priors={"duplication": (0, 0.3), "loss": (0, 0.4), "origination": (0, 1.5)},
        rounds=2, n_particles=25, initial_families=10, seed=1, max_attempts_factor=40)
    assert fit.sample_weights is not None
    assert len(fit.accepted) == 25                 # population size stays fixed
    assert set(fit.posterior) == {"duplication", "loss", "origination"}
    assert fit.n_simulations >= 25
    o = fit.summary()["origination"]
    assert o["lo95"] <= o["median"] <= o["hi95"]   # weighted quantiles ordered


def test_smc_requires_uniform_priors():
    tree, emp = _smc_setup()
    with pytest.raises(ValueError):
        z.match_profiles_smc(tree, emp, priors={"duplication": z.Gamma(2, 0.1)},
                             rounds=2, n_particles=10, seed=1)


def test_smc_reproducible():
    tree, emp = _smc_setup()
    kw = dict(priors={"duplication": (0, 0.3), "origination": (0, 1.5)}, rounds=2,
              n_particles=20, initial_families=10, seed=1, max_attempts_factor=40)
    a = z.match_profiles_smc(tree, emp, **kw)
    b = z.match_profiles_smc(tree, emp, **kw)
    assert np.array_equal(a.samples, b.samples)
    assert np.allclose(a.sample_weights, b.sample_weights)
