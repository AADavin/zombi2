"""Codon substitution models (GY94 / MG94).

Validated against analytic properties and simulation:

* **structure** — the state space is exactly the 61 sense codons, ``Q`` connects only single-nucleotide
  neighbours, stop codons are unreachable, rows sum to zero and the mean rate is 1;
* **reversibility** — detailed balance ``π_i Q_ij = π_j Q_ji`` holds to machine precision, and
  ``P(t)=exp(Qt)`` matches a reference matrix exponential;
* **dN/dS** — :func:`expected_dnds` recovers the injected ``omega`` exactly from the matrices, and a
  synonymous/non-synonymous substitution count on evolved sequences recovers it statistically;
* **kappa** — under uniform codon frequencies the synonymous transition/transversion rate ratio is
  ``kappa``;
* **stationarity** — a long branch drives the base composition to the model's target frequencies.
"""

from collections import Counter

import numpy as np
import pytest

from zombi2.sequences.codon_models import (
    GENETIC_CODE, SENSE_CODONS, STOP_CODONS, expected_dnds, f3x4, gy94, mg94, translate,
)
from zombi2.sequences.models import evolve_on_tree, is_codon_model, make_model


class _N:
    """A minimal tree node (gid + children) for the engine tests."""

    def __init__(self, gid, children=()):
        self.gid, self.children = gid, list(children)


def _reference_expm(A):
    """A scipy-free reference ``exp(A)`` (scaling-and-squaring with a Taylor series) for the P(t) check."""
    n = 0
    M = A.astype(float).copy()
    while np.abs(M).sum(1).max() > 0.5:
        M /= 2.0
        n += 1
    E = np.eye(A.shape[0])
    term = np.eye(A.shape[0])
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        for k in range(1, 40):
            term = term @ M / k
            E = E + term
        for _ in range(n):
            E = E @ E
    return E


def _diff_pos(c1, c2):
    return [k for k in range(3) if c1[k] != c2[k]]


# --------------------------------------------------------------------------- #
# The genetic code / alphabet
# --------------------------------------------------------------------------- #
def test_sense_codon_alphabet():
    assert len(SENSE_CODONS) == 61
    assert STOP_CODONS == frozenset({"TAA", "TAG", "TGA"})
    # a few standard-code spot checks
    assert translate("ATG") == "M" and translate("TGG") == "W" and translate("TTT") == "F"
    assert GENETIC_CODE["TAA"] == "*"


@pytest.mark.parametrize("mk", [gy94, mg94])
def test_alphabet_and_shape(mk):
    m = mk(kappa=2.5, omega=0.4, freqs=(0.3, 0.2, 0.2, 0.3))
    assert m.k == 61
    assert m.unit == 3                     # three nucleotides per codon state
    assert tuple(m.alphabet) == SENSE_CODONS


# --------------------------------------------------------------------------- #
# Matrix structure and reversibility
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mk", [gy94, mg94])
def test_matrix_structure(mk):
    """Q connects only single-nucleotide neighbours; rows sum to 0; mean rate is 1; π > 0."""
    m = mk(kappa=3.0, omega=0.5, freqs=(0.3, 0.2, 0.2, 0.3))
    Q, pi = m.Q, m.stationary
    assert (pi > 0).all() and np.isclose(pi.sum(), 1.0)
    assert np.allclose(Q.sum(axis=1), 0.0, atol=1e-12)             # valid rate matrix
    assert np.isclose(-(pi * np.diag(Q)).sum(), 1.0)              # normalised to 1 subst/site
    off = ~np.eye(61, dtype=bool)
    assert (Q[off] >= 0).all()                                    # off-diagonals non-negative
    for i, c1 in enumerate(SENSE_CODONS):
        for j, c2 in enumerate(SENSE_CODONS):
            if i != j and len(_diff_pos(c1, c2)) != 1:
                assert Q[i, j] == 0.0                             # multi-step changes forbidden


@pytest.mark.parametrize("mk", [gy94, mg94])
def test_detailed_balance(mk):
    """Time-reversibility: π_i Q_ij == π_j Q_ji (so the reversible eigendecomposition is valid)."""
    m = mk(kappa=4.0, omega=0.3, freqs=f3x4([[0.4, 0.1, 0.2, 0.3],
                                             [0.2, 0.3, 0.3, 0.2],
                                             [0.1, 0.4, 0.2, 0.3]]) if mk is gy94
           else [[0.4, 0.1, 0.2, 0.3], [0.2, 0.3, 0.3, 0.2], [0.1, 0.4, 0.2, 0.3]])
    flux = m.stationary[:, None] * m.Q
    assert np.abs(flux - flux.T).max() < 1e-12


@pytest.mark.parametrize("mk", [gy94, mg94])
def test_p_matrix_matches_expm_and_is_stochastic(mk):
    m = mk(kappa=2.0, omega=0.6, freqs=(0.3, 0.2, 0.25, 0.25))
    P = m.p_matrix(0.37)
    assert np.allclose(P.sum(axis=1), 1.0, atol=1e-9)             # rows are distributions
    assert (P >= 0).all()
    assert np.allclose(m.p_matrix(0.0), np.eye(61), atol=1e-12)   # P(0) = I
    assert np.abs(P - _reference_expm(m.Q * 0.37)).max() < 1e-9   # matches exp(Qt)


def test_p_matrix_raises_no_floating_point_warnings():
    """The 61-state matmul must not leak the spurious BLAS divide/overflow RuntimeWarnings."""
    import warnings
    m = gy94(kappa=2.0, omega=0.5, freqs=(0.3, 0.2, 0.2, 0.3))
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        m.p_matrix(0.1)
        m.p_matrix(1.5)


# --------------------------------------------------------------------------- #
# dN/dS (omega)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mk", [gy94, mg94])
@pytest.mark.parametrize("omega", [0.1, 0.5, 1.0, 2.0])
def test_expected_dnds_recovers_omega(mk, omega):
    """The matrix-level dN/dS (vs the ω=1 twin) returns the injected omega exactly."""
    freqs = (0.3, 0.2, 0.25, 0.25)
    sel = mk(kappa=3.0, omega=omega, freqs=freqs)
    neu = mk(kappa=3.0, omega=1.0, freqs=freqs)
    assert expected_dnds(sel, neu) == pytest.approx(omega, abs=1e-9)


def _count_syn_nonsyn(model, root_dna, t, seed):
    """Evolve ``root_dna`` down a single branch of length ``t`` and count syn / non-syn codon changes."""
    a = _N("a")
    root = _N("r", [a])
    ev = evolve_on_tree(root, {a: t}, model, np.random.default_rng(seed), root_seq=root_dna)
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


@pytest.mark.parametrize("mk", [gy94, mg94])
def test_omega_recovered_from_simulated_substitutions(mk):
    """dN/dS = (N_sel/S_sel)/(N_neu/S_neu) on evolved sequences recovers the injected omega ≈ 0.2.

    The neutral (ω=1) run supplies the synonymous/non-synonymous *opportunity*; the selected run's
    depressed non-synonymous count divided by that opportunity is the realised dN/dS. Averaged over
    several seeds on a moderate branch (few multiple hits), it lands near the injected 0.2.
    """
    freqs = (0.32, 0.18, 0.22, 0.28)
    omega = 0.2
    t = 0.3
    rng = np.random.default_rng(0)
    root_dna = "".join(rng.choice(SENSE_CODONS, size=6000))
    sel = mk(kappa=2.0, omega=omega, freqs=freqs)
    neu = mk(kappa=2.0, omega=1.0, freqs=freqs)
    est = []
    for seed in range(6):
        s1, n1 = _count_syn_nonsyn(sel, root_dna, t, seed)
        s0, n0 = _count_syn_nonsyn(neu, root_dna, t, 100 + seed)
        est.append((n1 / s1) / (n0 / s0))
    assert np.mean(est) == pytest.approx(omega, abs=0.06)


@pytest.mark.parametrize("mk", [gy94, mg94])
def test_omega_one_is_neutral(mk):
    m = mk(kappa=2.5, omega=1.0, freqs=(0.3, 0.2, 0.25, 0.25))
    assert expected_dnds(m, m) == pytest.approx(1.0, abs=1e-9)


# --------------------------------------------------------------------------- #
# kappa
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mk", [gy94, mg94])
def test_kappa_is_synonymous_transition_transversion_ratio(mk):
    """Under uniform codon frequencies, every synonymous transition rate is ``kappa`` × every
    synonymous transversion rate (the ω factor cancels: both are synonymous)."""
    kappa = 4.0
    m = mk(kappa=kappa, omega=0.5, freqs=None)     # None -> uniform (equal codon freqs)
    purines = set("AG")
    ts, tv = [], []
    for i, c1 in enumerate(SENSE_CODONS):
        for j, c2 in enumerate(SENSE_CODONS):
            d = _diff_pos(c1, c2)
            if len(d) != 1 or GENETIC_CODE[c1] != GENETIC_CODE[c2]:
                continue                            # single-nt synonymous neighbours only
            a, b = c1[d[0]], c2[d[0]]
            (ts if (a in purines) == (b in purines) else tv).append(m.Q[i, j])
    assert ts and tv
    assert np.ptp(ts) < 1e-12 and np.ptp(tv) < 1e-12   # each class is a single rate
    assert ts[0] / tv[0] == pytest.approx(kappa, rel=1e-9)


# --------------------------------------------------------------------------- #
# GY94 vs MG94, stationarity, engine integration
# --------------------------------------------------------------------------- #
def test_gy94_and_mg94_differ():
    """GY94 weights by target-codon frequency, MG94 by target-nucleotide frequency: same knobs,
    genuinely different matrices (they coincide only under uniform frequencies)."""
    freqs = (0.4, 0.1, 0.2, 0.3)
    g = gy94(kappa=2.0, omega=0.5, freqs=freqs)
    mgg = mg94(kappa=2.0, omega=0.5, freqs=freqs)
    assert not np.allclose(g.Q, mgg.Q)


def test_gy94_stationary_composition_recovered():
    """A long branch drives the base composition to the model's equilibrium.

    The equilibrium is the codon stationary ``π``; its *marginal* base composition is what a long-run
    sequence recovers. That marginal is close to — but not exactly — the input F1×4 base frequencies,
    because excluding the three stop codons and renormalising skews it; so we compare against the
    marginal of ``π`` itself, the honest target.
    """
    base = np.array([0.4, 0.1, 0.2, 0.3])
    m = gy94(kappa=2.0, omega=0.6, freqs=tuple(base))
    order = "ACGT"
    idx = {b: i for i, b in enumerate(order)}
    expected = np.zeros(4)                                          # marginal base comp of pi
    for codon, p in zip(SENSE_CODONS, m.stationary):
        for ch in codon:
            expected[idx[ch]] += p
    expected /= 3.0
    assert not np.allclose(expected, base, atol=1e-4)              # stop-exclusion really does skew it
    a = _N("a")
    root = _N("r", [a])
    rng = np.random.default_rng(1)
    root_dna = "".join(rng.choice(SENSE_CODONS, size=40000))       # arbitrary start
    seqs = evolve_on_tree(root, {a: 30.0}, m, rng, root_seq=root_dna)   # long branch -> equilibrium
    c = Counter(seqs["a"])
    recovered = np.array([c[b] / len(seqs["a"]) for b in order])
    assert np.allclose(recovered, expected, atol=0.01)


def test_zero_branch_copies_parent_and_root_seq_roundtrips():
    a = _N("a")           # t=0 branch: must copy the root exactly
    root = _N("r", [a])
    m = gy94(kappa=2.0, omega=0.5)
    rng = np.random.default_rng(2)
    root_dna = "".join(rng.choice(SENSE_CODONS, size=100))
    seqs = evolve_on_tree(root, {a: 0.0}, m, rng, root_seq=root_dna)
    assert seqs["r"] == root_dna and seqs["a"] == root_dna


@pytest.mark.parametrize("mk", [gy94, mg94])
def test_no_stop_codons_ever_appear(mk):
    m = mk(kappa=3.0, omega=1.5, freqs=(0.3, 0.2, 0.2, 0.3))       # even positive selection
    a, b = _N("a"), _N("b")
    root = _N("r", [a, b])
    rng = np.random.default_rng(4)
    seqs = evolve_on_tree(root, {a: 2.0, b: 2.0}, m, rng, length=500)
    for dna in seqs.values():
        codons = {dna[i:i + 3] for i in range(0, len(dna), 3)}
        assert not (codons & STOP_CODONS)


def test_codon_seq_length_is_codons_via_length():
    m = gy94(kappa=2.0, omega=0.5)
    a = _N("a")
    root = _N("r", [a])
    seqs = evolve_on_tree(root, {a: 0.5}, m, np.random.default_rng(5), length=120)
    assert len(seqs["a"]) == 360                                   # 120 codons -> 360 nt


# --------------------------------------------------------------------------- #
# make_model routing
# --------------------------------------------------------------------------- #
def test_make_model_routes_codon_with_omega():
    assert is_codon_model("gy94") and is_codon_model("MG94") and not is_codon_model("gtr")
    m = make_model("gy94", kappa=2.0, omega=0.25, freqs=(0.25, 0.25, 0.25, 0.25))
    assert m.name == "GY94" and m.k == 61
    # omega recovers through make_model too
    neu = make_model("mg94", omega=1.0)
    sel = make_model("mg94", omega=0.3)
    assert expected_dnds(sel, neu) == pytest.approx(0.3, abs=1e-9)


def test_make_model_warns_on_omega_for_non_codon_model():
    with pytest.warns(UserWarning, match="does not use"):
        make_model("hky85", omega=0.5)
    with pytest.warns(UserWarning, match="does not use"):
        make_model("gy94", rates=(1, 1, 1, 1, 1, 1))              # gy94 has no exchangeabilities
