"""P2 codon-level Halpern-Bruno mutation-selection: the codon kernel + emergent dN/dS.

Exercised with a FixedProfileCritic (a known injected amino-acid preference), so the suite runs
without torch/esm.
"""
from __future__ import annotations

import ast
import inspect

import numpy as np
import pytest

import zombi2.experimental as ex
from zombi2.experimental.codon_selection import (
    GENETIC_CODE, SENSE_CODONS, STOP_CODONS, CodonSelection, _CODON_AA, _codon_mutation,
    _codon_site_model, calibrate_beta, translate,
)
from zombi2.experimental.selection import FixedProfileCritic
from zombi2.sequences.models import AMINO_ACIDS, hky85

_AA = {a: i for i, a in enumerate(AMINO_ACIDS)}
_AA_TO_CODON = {}
for _c in SENSE_CODONS:                       # one representative codon per amino acid
    _AA_TO_CODON.setdefault(GENETIC_CODE[_c], _c)


def _peaked(target: str, hi: float = 0.95) -> np.ndarray:
    p = np.full((len(target), 20), (1.0 - hi) / 19.0)
    for i, a in enumerate(target):
        p[i, _AA[a]] = hi
    return p


def _dna(protein: str) -> str:
    return "".join(_AA_TO_CODON[a] for a in protein)


class _N:
    def __init__(self, gid, children=()):
        self.gid = gid
        self.children = list(children)


# --------------------------------------------------------------------------- #
# genetic code + translation + lifecycle
# --------------------------------------------------------------------------- #
def test_genetic_code_shape():
    assert len(GENETIC_CODE) == 64
    assert len(SENSE_CODONS) == 61 and len(STOP_CODONS) == 3
    assert STOP_CODONS == {"TAA", "TAG", "TGA"}


def test_translate_known_codons():
    assert translate("ATGTGGTTT") == "MWF"           # Met, Trp, Phe
    assert translate("GCTGCCGCAGCG") == "AAAA"        # all alanine (synonymous)
    with pytest.raises(ValueError):
        translate("ATGT")                             # not a multiple of 3


def test_codon_selection_is_experimental():
    ex._warned.discard("CodonSelection")
    with pytest.warns(ex.ExperimentalWarning, match="CodonSelection"):
        CodonSelection(FixedProfileCritic(_peaked("ACDE")))


@pytest.mark.parametrize("bad", [-1.0, float("nan"), float("inf")])
def test_beta_rejects_non_finite_or_negative(bad):
    with pytest.raises(ValueError, match="beta"):
        CodonSelection(FixedProfileCritic(_peaked("ACDE")), beta=bad)


def test_rejects_non_nucleotide_model():
    from zombi2.sequences.models import lg
    with pytest.raises(ValueError, match="nucleotide"):
        CodonSelection(FixedProfileCritic(_peaked("ACDE")), nuc_model=lg())


def test_stop_codon_in_input_is_rejected():
    sel = CodonSelection(FixedProfileCritic(_peaked("MM")))
    root = _N("root", [_N("tip")])
    with pytest.raises(ValueError, match="stop codon"):
        sel.evolve_coding_family(root, {root: 0.0, root.children[0]: 1.0},
                                 "ATGTAA", rng=np.random.default_rng(0))


# --------------------------------------------------------------------------- #
# the codon kernel
# --------------------------------------------------------------------------- #
def test_beta_zero_reduces_to_the_neutral_codon_model():
    nuc = hky85()
    mu, pi_mut = _codon_mutation(nuc)
    Qn = mu.copy()
    np.fill_diagonal(Qn, -Qn.sum(1))
    Qn = Qn / -(pi_mut * np.diag(Qn)).sum()                          # neutral codon Q, mean rate 1
    m0 = _codon_site_model(mu, pi_mut, np.ones(20) / 20.0, 0.0)
    assert np.allclose(m0.Q, Qn, atol=1e-10)
    assert np.allclose(m0.stationary, pi_mut, atol=1e-12)


def test_codon_site_model_is_reversible():
    nuc = hky85()
    mu, pi_mut = _codon_mutation(nuc)
    m = _codon_site_model(mu, pi_mut, _peaked("W", hi=0.9)[0], 4.0)
    flux = m.stationary[:, None] * m.Q
    assert np.allclose(flux, flux.T, atol=1e-12)


# --------------------------------------------------------------------------- #
# the headline: emergent dN/dS
# --------------------------------------------------------------------------- #
def test_dnds_is_one_at_beta_zero_and_decreases_monotonically_toward_zero():
    critic = FixedProfileCritic(_peaked("MQIFVKTLTGKTITLEVE", hi=0.9))
    prot = "MQIFVKTLTGKTITLEVE"
    betas = [0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0]
    omegas = [CodonSelection(critic, beta=b).dnds(prot) for b in betas]
    assert abs(omegas[0] - 1.0) < 1e-9, omegas                       # neutral => omega == 1
    assert all(omegas[i + 1] < omegas[i] for i in range(len(betas) - 1)), omegas   # strictly down
    assert omegas[-1] < 1e-3, omegas                                 # -> 0 under strong selection


def test_calibrate_beta_inverts_dnds():
    # ask for a target dN/dS, get the beta whose model yields it (the usable inverse of dnds)
    critic = FixedProfileCritic(_peaked("MQIFVKTLTGKTITLEVE", hi=0.9))
    prot = "MQIFVKTLTGKTITLEVE"
    for target in (0.7, 0.3, 0.1):
        beta = calibrate_beta(critic, prot, target)
        assert abs(CodonSelection(critic, beta=beta).dnds(prot) - target) < 1e-3, (target, beta)
    with pytest.raises(ValueError, match="target_dnds"):
        calibrate_beta(critic, prot, 1.5)                       # outside (0, 1)
    with pytest.raises(ValueError, match="did not reach tol"):
        calibrate_beta(critic, prot, 0.3, max_iter=2)          # too few iterations -> raise, not a silent miss


def test_synonymous_flux_stays_neutral_dS_is_one():
    critic = FixedProfileCritic(_peaked("MQIFVKTLTG", hi=0.9))
    for beta in (0.0, 2.0, 8.0):
        assert abs(CodonSelection(critic, beta=beta).dnds_syn_check("MQIFVKTLTG") - 1.0) < 1e-9


# --------------------------------------------------------------------------- #
# behaviour: inject -> recover at the DNA level
# --------------------------------------------------------------------------- #
def test_coding_evolution_recovers_the_injected_protein():
    target = "MQIFVKTLTGKTITLE"
    sel = CodonSelection(FixedProfileCritic(_peaked(target, hi=0.95)), beta=6.0)
    root = _N("root", [_N("tip")])
    subst = {root: 0.0, root.children[0]: 50.0}                      # long branch -> reach stationary
    out = sel.evolve_coding_family(root, subst, _dna("A" * len(target)),
                                   rng=np.random.default_rng(0))
    recovered = translate(out["tip"])
    assert len(out["tip"]) == 3 * len(target)                        # stays in frame
    frac = np.mean([a == b for a, b in zip(recovered, target)])
    assert frac > 0.9, f"only recovered {frac:.0%} of the injected protein"


def test_high_beta_near_delta_stays_stable_and_recovers():
    # the near-delta stationary of strong selection (which corrupts a sqrt(pi) eigendecomposition
    # AND a general eig of Q) must still give a valid, converging process under the scipy expm kernel
    target = "WKFED"
    for beta in (20.0, 100.0):
        sel = CodonSelection(FixedProfileCritic(_peaked(target, hi=0.999)), beta=beta)
        root = _N("root", [_N("tip")])
        out = sel.evolve_coding_family(root, {root: 0.0, root.children[0]: 50.0},
                                       _dna("A" * len(target)), rng=np.random.default_rng(0))
        assert translate(out["tip"]) == target


def test_evolution_is_deterministic_and_families_reproducible():
    sel = CodonSelection(FixedProfileCritic(_peaked("ACDEFG")), beta=3.0)
    root = _N("root", [_N("a"), _N("b")])
    subst = {root: 0.0, root.children[0]: 2.0, root.children[1]: 2.0}
    dna = _dna("ACDEFG")
    a = sel.evolve_coding_family(root, subst, dna, rng=np.random.default_rng(5))
    b = sel.evolve_coding_family(root, subst, dna, rng=np.random.default_rng(5))
    assert a == b
    nt = {"f": {"complete": (root, subst), "extant": None}}
    o1 = sel.evolve_coding_families(nt, {"f": dna}, seed=1)
    assert sel.evolve_coding_families(nt, {"f": dna}, seed=1) == o1


def _single_nt_neighbours(codon):
    return [codon[:k] + b + codon[k + 1:] for k in range(3) for b in "ACGT" if b != codon[k]]


def test_genetic_code_matches_the_standard():
    # every amino acid + stop appears; spot-check the tricky degeneracies (a wrong _AAS char is caught)
    assert set(GENETIC_CODE.values()) == set(AMINO_ACIDS) | {"*"}
    ref = {"TTA": "L", "TTG": "L", "CTA": "L", "CTG": "L",        # Leu 6-fold
           "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",        # Arg 6-fold
           "TCA": "S", "AGT": "S", "AGC": "S",                    # Ser 6-fold (split box)
           "ATT": "I", "ATA": "I", "ATG": "M",                    # Ile vs Met
           "TAA": "*", "TAG": "*", "TGA": "*", "TGG": "W"}        # stops vs Trp
    for c, aa in ref.items():
        assert GENETIC_CODE[c] == aa, (c, GENETIC_CODE[c], aa)


def test_mutation_backbone_matches_the_nucleotide_model_under_unequal_frequencies():
    # equal-freq hky85 makes codon mu symmetric + pi_mut uniform, hiding directionality bugs;
    # an unequal-frequency model exposes them, so build mu and pi_mut INDEPENDENTLY and compare.
    nuc = hky85(kappa=2.5, freqs=(0.4, 0.3, 0.2, 0.1))
    mu, pi_mut = _codon_mutation(nuc)
    nb = {b: i for i, b in enumerate("ACGT")}
    exp_pi = np.array([nuc.stationary[nb[c[0]]] * nuc.stationary[nb[c[1]]] * nuc.stationary[nb[c[2]]]
                       for c in SENSE_CODONS])
    assert np.allclose(pi_mut, exp_pi / exp_pi.sum())
    for i, c1 in enumerate(SENSE_CODONS):
        for j, c2 in enumerate(SENSE_CODONS):
            diffs = [k for k in range(3) if c1[k] != c2[k]]
            exp = nuc.Q[nb[c1[diffs[0]]], nb[c2[diffs[0]]]] if len(diffs) == 1 else 0.0
            assert abs(mu[i, j] - exp) < 1e-12, (c1, c2)
    # beta=0 kernel + reversibility must also hold under this asymmetric backbone
    Qn = mu.copy()
    np.fill_diagonal(Qn, -Qn.sum(1))
    Qn = Qn / -(pi_mut * np.diag(Qn)).sum()
    m0 = _codon_site_model(mu, pi_mut, np.ones(20) / 20.0, 0.0)
    assert np.allclose(m0.Q, Qn, atol=1e-10) and np.allclose(m0.stationary, pi_mut, atol=1e-12)
    m = _codon_site_model(mu, pi_mut, _peaked("W", hi=0.9)[0], 4.0)
    flux = m.stationary[:, None] * m.Q
    assert np.allclose(flux, flux.T, atol=1e-12)


def test_expm_kernel_is_a_valid_transition_matrix():
    mu, pi_mut = _codon_mutation(hky85())
    m = _codon_site_model(mu, pi_mut, _peaked("W", hi=0.9)[0], 4.0)
    n = len(SENSE_CODONS)
    P = m.p_matrix(0.3)
    assert np.allclose(P.sum(1), 1.0) and (P >= 0).all()             # a valid stochastic matrix
    assert np.allclose(m.p_matrix(0.0), np.eye(n))                   # P(0) = I
    assert np.allclose(m.stationary @ P, m.stationary)              # stationary is preserved
    assert np.allclose(m.p_matrix(500.0), m.stationary[None, :], atol=1e-6)   # P(inf) -> stationary


def test_nonsense_mutations_to_stop_codons_are_excluded():
    # stop-exclusion is a modelling choice, not just an alphabet omission: a sense codon one nt from
    # a stop has zero flux toward it, so its sense-neighbour count drops by the adjacent stops.
    mu, _ = _codon_mutation(hky85())
    for codon in ("TGG", "TAC"):                                     # Trp and Tyr are stop-adjacent
        i = SENSE_CODONS.index(codon)
        nbrs = _single_nt_neighbours(codon)
        stops = [c for c in nbrs if c in STOP_CODONS]
        assert len(stops) > 0                                        # really stop-adjacent
        assert int((mu[i] > 0).sum()) == len(nbrs) - len(stops)      # stops carry no flux


def test_simulated_dnds_matches_the_analytic_prediction():
    # tie the analytic headline dnds() to the actual expm simulator: evolve many independent codon
    # sites a short branch under BOTH the selected and the neutral (beta=0) model from the SAME roots,
    # count synonymous vs non-synonymous end-state changes, and form
    #   omega = (nonsyn/syn | selected) / (nonsyn/syn | neutral),
    # which should match dnds(). Uses an unequal-frequency backbone (also covers mu directionality).
    rng = np.random.default_rng(0)
    L, beta, t, reps = 500, 1.2, 0.1, 15
    prot = "".join(rng.choice(list(AMINO_ACIDS.replace("W", "")), L))   # skip W (no synonymous site)
    nuc = hky85(kappa=2.5, freqs=(0.4, 0.3, 0.2, 0.1))
    critic = FixedProfileCritic(_peaked(prot, hi=0.8))
    sel = CodonSelection(critic, beta=beta, nuc_model=nuc)
    ms = sel._site_models(prot)
    mn = CodonSelection(critic, beta=0.0, nuc_model=nuc)._site_models(prot)
    ns_s = sy_s = ns_n = sy_n = 0
    for _ in range(reps):
        for msi, mni in zip(ms, mn):
            r = rng.choice(len(SENSE_CODONS), p=msi.stationary)
            cs = rng.choice(len(SENSE_CODONS), p=msi.p_matrix(t)[r])
            cn = rng.choice(len(SENSE_CODONS), p=mni.p_matrix(t)[r])
            if cs != r:
                if _CODON_AA[cs] == _CODON_AA[r]:
                    sy_s += 1
                else:
                    ns_s += 1
            if cn != r:
                if _CODON_AA[cn] == _CODON_AA[r]:
                    sy_n += 1
                else:
                    ns_n += 1
    omega_emp = (ns_s / sy_s) / (ns_n / sy_n)
    assert abs(omega_emp - sel.dnds(prot)) < 0.12, (omega_emp, sel.dnds(prot))


# --------------------------------------------------------------------------- #
# packaging
# --------------------------------------------------------------------------- #
def test_module_has_no_top_level_ml_imports():
    from zombi2.experimental import codon_selection
    tree = ast.parse(inspect.getsource(codon_selection))
    top: set[str] = set()
    for node in tree.body:
        if isinstance(node, ast.Import):
            top.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            top.add(node.module.split(".")[0])
    assert "torch" not in top and "esm" not in top, top


def test_exports_stay_in_the_experimental_namespace():
    import zombi2
    for name in ("CodonSelection", "translate"):
        assert name in ex.__all__ and hasattr(ex, name)
        assert not hasattr(zombi2, name), f"{name} leaked into the top-level zombi2 namespace"
