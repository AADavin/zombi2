"""Codon site models — dN/dS varies among sites (M1a / M2a / M3 / M7 / M8).

Each model is a mixture of GY94/MG94 matrices that share one mutation process (hence one stationary
distribution) and differ only in ``ω``, normalised on a **single shared scale** so that purifying
classes genuinely evolve slower. The suite pins:

* the ``Beta`` discretisation (M7/M8) averages to the beta mean ``p/(p+q)``;
* the mixture's genome-wide ``dN/dS`` equals the proportion-weighted mean ``ω`` — checked both by an
  exact flux computation and by counting substitutions on evolved sequences;
* the shared scale preserves across-class rate heterogeneity (purifying < neutral exit rate);
* reductions (M1a ``p0=0`` → neutral, M8 ``p0=1`` → M7, single-class M3 → M0) and the engine
  integration (in-frame coding DNA, no stop codons, mutually exclusive with gamma).
"""

import numpy as np
import pytest

from zombi2.sequences.codon_models import (
    GENETIC_CODE, SENSE_CODONS, CodonSiteModel, beta_category_omegas, gy94,
    is_codon_site_model, m1a, m2a, m3, m7, m8, make_codon_site_model, _mean_rate, _syn_masks,
)
from zombi2.sequences.models import GammaRates, evolve_on_tree


class _N:
    def __init__(self, gid, children=()):
        self.gid, self.children = gid, list(children)


_FREQS = (0.3, 0.2, 0.25, 0.25)


def _flux_dnds(mix):
    """Exact mixture dN/dS from component fluxes vs a neutral (ω=1) twin — no simulation."""
    syn, non = _syn_masks()
    pi = mix.stationary
    n = sum(p * float((pi[:, None] * c.Q * non).sum()) for p, c in zip(mix.proportions, mix.components))
    s = sum(p * float((pi[:, None] * c.Q * syn).sum()) for p, c in zip(mix.proportions, mix.components))
    neu = gy94(2.0, 1.0, freqs=_FREQS)
    n0 = float((pi[:, None] * neu.Q * non).sum())
    s0 = float((pi[:, None] * neu.Q * syn).sum())
    return (n / n0) * (s0 / s)


def _count_syn_nonsyn(model, root_dna, t, seed):
    a = _N("a")
    ev = evolve_on_tree(_N("r", [a]), {a: t}, model, np.random.default_rng(seed), root_seq=root_dna)
    des = ev["a"]
    syn = non = 0
    for i in range(0, len(root_dna), 3):
        c1, c2 = root_dna[i:i + 3], des[i:i + 3]
        if c1 == c2:
            continue
        if GENETIC_CODE[c1] == GENETIC_CODE[c2]:
            syn += 1
        else:
            non += 1
    return syn, non


# --------------------------------------------------------------------------- #
# Beta discretisation (M7 / M8)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("p,q,ncat", [(2.0, 2.0, 4), (0.5, 2.0, 6), (3.0, 1.0, 8), (1.0, 1.0, 5)])
def test_beta_categories_average_to_beta_mean(p, q, ncat):
    om = beta_category_omegas(p, q, ncat)
    assert om.shape == (ncat,)
    assert (om >= 0).all() and (om <= 1).all()          # ω in [0, 1]
    assert np.all(np.diff(om) >= 0)                     # increasing
    assert om.mean() == pytest.approx(p / (p + q), abs=1e-6)   # discretisation is unbiased


# --------------------------------------------------------------------------- #
# mean_omega and the exact flux oracle
# --------------------------------------------------------------------------- #
def _models():
    return {
        "M1a": m1a(2.0, p0=0.6, omega0=0.1, freqs=_FREQS),
        "M2a": m2a(2.0, p0=0.5, omega0=0.1, p1=0.3, omega2=2.5, freqs=_FREQS),
        "M3": m3(2.0, omegas=[0.05, 0.5, 2.0], proportions=[0.5, 0.3, 0.2], freqs=_FREQS),
        "M7": m7(2.0, beta_p=0.6, beta_q=1.5, ncat=5, freqs=_FREQS),
        "M8": m8(2.0, beta_p=0.6, beta_q=1.5, p0=0.85, omega_s=3.0, ncat=5, freqs=_FREQS),
    }


@pytest.mark.parametrize("name", ["M1a", "M2a", "M3", "M7", "M8"])
def test_mean_omega_equals_weighted_mean(name):
    mix = _models()[name]
    assert mix.mean_omega == pytest.approx(float((mix.proportions * mix.omegas).sum()))
    assert mix.proportions.sum() == pytest.approx(1.0)


@pytest.mark.parametrize("name", ["M1a", "M2a", "M3", "M7", "M8"])
def test_flux_dnds_equals_mean_omega(name):
    """The shared-scale construction makes the exact genome-wide dN/dS equal the mean ω."""
    mix = _models()[name]
    assert _flux_dnds(mix) == pytest.approx(mix.mean_omega, abs=1e-9)


@pytest.mark.parametrize("name", ["M1a", "M2a", "M8"])
def test_omega_recovered_from_simulated_substitutions(name):
    """Counting syn/non-syn substitutions on evolved sequences recovers the mixture's mean dN/dS.

    Uses a short branch (few multiple hits) so net codon differences ≈ substitution counts; the
    high-ω classes would otherwise saturate and bias the estimate down."""
    mix = _models()[name]
    neu = gy94(2.0, 1.0, freqs=_FREQS)
    rng = np.random.default_rng(0)
    root = "".join(rng.choice(SENSE_CODONS, size=16000))
    est = []
    for seed in range(5):
        s1, n1 = _count_syn_nonsyn(mix, root, 0.1, seed)
        s0, n0 = _count_syn_nonsyn(neu, root, 0.1, 100 + seed)
        est.append((n1 / s1) / (n0 / s0))
    assert np.mean(est) == pytest.approx(mix.mean_omega, abs=0.06)


# --------------------------------------------------------------------------- #
# Shared-scale heterogeneity, positive selection, reductions
# --------------------------------------------------------------------------- #
def test_shared_scale_preserves_rate_heterogeneity():
    """A purifying class has a lower total substitution rate than the neutral class — the whole point
    of one shared scale (per-class normalisation would erase this and break dN/dS = mean ω)."""
    mix = m2a(2.0, p0=0.5, omega0=0.1, p1=0.3, omega2=2.5, freqs=_FREQS)
    pi = mix.stationary
    rates = [_mean_rate(c.Q, pi) for c in mix.components]        # classes: ω0<1, ω=1, ω2>1
    assert rates[0] < rates[1] < rates[2]                        # purifying slower, positive faster
    assert float((mix.proportions * np.array(rates)).sum()) == pytest.approx(1.0)   # mixture mean = 1


def test_positive_selection_class_present_only_where_expected():
    assert m2a(2.0).omegas.max() > 1.0 and m8(2.0, beta_p=1, beta_q=2).omegas.max() > 1.0
    assert m1a(2.0).omegas.max() <= 1.0
    assert m7(2.0, beta_p=1, beta_q=2).omegas.max() <= 1.0


def test_m1a_p0_zero_is_neutral():
    mix = m1a(2.0, p0=0.0, freqs=_FREQS)                          # all sites neutral
    assert mix.mean_omega == pytest.approx(1.0)
    assert _flux_dnds(mix) == pytest.approx(1.0, abs=1e-9)


def test_m8_p0_one_reduces_to_m7():
    a = m8(2.0, beta_p=1.2, beta_q=2.0, p0=1.0, ncat=6, freqs=_FREQS)
    b = m7(2.0, beta_p=1.2, beta_q=2.0, ncat=6, freqs=_FREQS)
    assert a.mean_omega == pytest.approx(b.mean_omega)           # the ω_s class carries proportion 0


def test_m3_single_class_matches_m0():
    mix = m3(3.0, omegas=[0.37], proportions=[1.0], freqs=_FREQS)
    m0 = gy94(3.0, 0.37, freqs=_FREQS)
    assert len(mix.components) == 1
    assert np.allclose(mix.components[0].Q, m0.Q)                # identical matrix (shared scale = own scale)
    assert np.allclose(mix.stationary, m0.stationary)


# --------------------------------------------------------------------------- #
# Engine integration
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name", ["M1a", "M2a", "M3", "M7", "M8"])
def test_evolve_in_frame_no_stops(name):
    mix = _models()[name]
    a, b = _N("a"), _N("b")
    seqs = evolve_on_tree(_N("r", [a, b]), {a: 1.5, b: 1.5}, mix, np.random.default_rng(3), length=400)
    for dna in seqs.values():
        assert len(dna) == 1200                                  # 400 codons -> 1200 nt
        codons = {dna[i:i + 3] for i in range(0, len(dna), 3)}
        assert codons <= set(SENSE_CODONS)                       # every codon is a sense codon


def test_mixture_and_gamma_are_mutually_exclusive():
    mix = m2a(2.0, freqs=_FREQS)
    a = _N("a")
    with pytest.raises(ValueError, match="per-site"):
        evolve_on_tree(_N("r", [a]), {a: 1.0}, mix, np.random.default_rng(1),
                       length=30, gamma=GammaRates(0.5))


def test_zero_branch_copies_parent():
    mix = m7(2.0, beta_p=1.0, beta_q=2.0, freqs=_FREQS)
    a = _N("a")
    rng = np.random.default_rng(2)
    root_dna = "".join(rng.choice(SENSE_CODONS, size=90))
    seqs = evolve_on_tree(_N("r", [a]), {a: 0.0}, mix, rng, root_seq=root_dna)
    assert seqs["a"] == root_dna


# --------------------------------------------------------------------------- #
# Construction / factory / validation
# --------------------------------------------------------------------------- #
def test_make_codon_site_model_routes_and_validates():
    assert is_codon_site_model("m1a") and is_codon_site_model("M8") and not is_codon_site_model("gy94")
    mix = make_codon_site_model("m2a", kappa=2.0, base="gy94", omega0=0.1, omega2=3.0)
    assert isinstance(mix, CodonSiteModel) and mix.name == "M2a"
    # mg94 base also works and shares its stationary across classes
    mg = make_codon_site_model("m7", base="mg94", beta_p=1.0, beta_q=2.0)
    assert all(np.allclose(c.stationary, mg.components[0].stationary) for c in mg.components)


def test_make_codon_site_model_missing_params_raise():
    with pytest.raises(ValueError, match="requires"):
        make_codon_site_model("m7")                              # needs beta_p / beta_q
    with pytest.raises(ValueError, match="requires"):
        make_codon_site_model("m3")                              # needs omegas / proportions


def test_codon_site_model_validation():
    good = gy94(2.0, 0.5, freqs=_FREQS)
    with pytest.raises(ValueError, match="sum to 1"):
        CodonSiteModel("bad", (good, good), [0.5, 0.6], [0.5, 1.0])
    with pytest.raises(ValueError, match="stationary"):
        CodonSiteModel("bad", (gy94(2.0, 0.5, freqs=_FREQS), gy94(2.0, 0.5, freqs=(0.4, 0.1, 0.2, 0.3))),
                       [0.5, 0.5], [0.5, 0.5])
    with pytest.raises(ValueError, match="p0"):
        m1a(2.0, p0=1.5)
    with pytest.raises(ValueError):
        m2a(2.0, p0=0.7, p1=0.7)                                 # p0 + p1 > 1
