"""P1 frozen mutation-selection: the Halpern-Bruno per-site kernel + the PLMSelection surface.

The core kernel is exercised with a FixedProfileCritic (a known injected preference), so the whole
suite runs without torch/esm; the ESM2 path is tested only when the optional deps are installed.
"""
from __future__ import annotations

import ast
import inspect

import numpy as np
import pytest

import zombi2.experimental as ex
from zombi2.experimental.selection import (
    FixedProfileCritic, PLMSelection, _site_models, _site_targets,
)
from zombi2.sequences.models import AMINO_ACIDS, lg

_AA = {a: i for i, a in enumerate(AMINO_ACIDS)}


class _N:
    """Minimal duck-typed tree node -- evolve_family only needs ``.gid`` and ``.children``."""

    def __init__(self, gid, children=()):
        self.gid = gid
        self.children = list(children)


def _peaked_profile(target: str, hi: float = 0.97) -> np.ndarray:
    """(L, 20) profile putting mass ``hi`` on ``target[i]`` at site ``i``, rest spread over the others."""
    p = np.full((len(target), 20), (1.0 - hi) / 19.0)
    for i, a in enumerate(target):
        p[i, _AA[a]] = hi
    return p


# --------------------------------------------------------------------------- #
# experimental lifecycle + validation
# --------------------------------------------------------------------------- #
def test_plmselection_is_experimental():
    ex._warned.discard("PLMSelection")
    with pytest.warns(ex.ExperimentalWarning, match="PLMSelection"):
        PLMSelection(FixedProfileCritic(_peaked_profile("ACDE")))


def test_fixedprofilecritic_is_experimental():
    ex._warned.discard("FixedProfileCritic")
    with pytest.warns(ex.ExperimentalWarning, match="FixedProfileCritic"):
        FixedProfileCritic(_peaked_profile("ACDE"))


def test_live_mode_not_yet_implemented():
    with pytest.raises(NotImplementedError):
        PLMSelection(FixedProfileCritic(_peaked_profile("ACDE")), mode="live")


def test_exports_stay_in_the_experimental_namespace():
    import zombi2
    for name in ("Critic", "ESM2Critic", "FixedProfileCritic", "PLMSelection"):
        assert name in ex.__all__ and hasattr(ex, name)
        assert not hasattr(zombi2, name), f"{name} leaked into the top-level zombi2 namespace"


def test_rejects_non_protein_model():
    from zombi2.sequences.models import jc69
    with pytest.raises(ValueError, match="amino-acid"):
        PLMSelection(FixedProfileCritic(_peaked_profile("ACDE")), model=jc69())


@pytest.mark.parametrize("bad", [-1.0, float("nan"), float("inf"), -float("inf")])
def test_beta_rejects_non_finite_or_negative(bad):
    with pytest.raises(ValueError, match="beta"):
        PLMSelection(FixedProfileCritic(_peaked_profile("ACDE")), beta=bad)


# --------------------------------------------------------------------------- #
# the kernel itself (not just its stationary)
# --------------------------------------------------------------------------- #
def test_beta_zero_kernel_reduces_exactly_to_the_base_model():
    # at beta=0 the Halpern-Bruno fixation factor is 1, so each site's Q IS the base model's Q
    base = lg()
    targets = _site_targets(_peaked_profile("MQIFVKTLTG"), base, 0.0)
    for m in _site_models(targets, base.Q, base.stationary, base.alphabet):
        assert np.allclose(m.Q, base.Q, atol=1e-10)
        assert np.allclose(m.stationary, base.stationary, atol=1e-12)


def test_site_model_is_reversible_with_the_target_stationary():
    # detailed balance: pi_a Q_ab == pi_b Q_ba (this is what makes it a valid Halpern-Bruno process)
    base = lg()
    targets = _site_targets(_peaked_profile("MQIFVKTLTG", hi=0.85), base, 3.0)
    m = _site_models(targets, base.Q, base.stationary, base.alphabet)[0]
    flux = targets[0][:, None] * m.Q
    assert np.allclose(flux, flux.T, atol=1e-12)


def test_selection_changes_the_rates_not_only_the_stationary():
    # the whole point of HB over a naive '+F' tilt: beta alters the flux, not just the equilibrium
    base = lg()
    prof = _peaked_profile("MQIFVKTLTG", hi=0.9)
    m0 = _site_models(_site_targets(prof, base, 0.0), base.Q, base.stationary, base.alphabet)[0]
    m5 = _site_models(_site_targets(prof, base, 5.0), base.Q, base.stationary, base.alphabet)[0]
    assert np.allclose(m0.Q, base.Q, atol=1e-10)
    assert not np.allclose(m5.Q, base.Q, atol=1e-3)


# --------------------------------------------------------------------------- #
# behaviour: inject -> recover, monotonic in beta, deterministic
# --------------------------------------------------------------------------- #
def test_frozen_recovers_injected_preference():
    # inject a peaked per-site preference, evolve a long branch under strong selection, and check the
    # tip recovers the injected residue at (almost) every site -- an inject -> recover validation.
    target = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQD"
    sel = PLMSelection(FixedProfileCritic(_peaked_profile(target, hi=0.97)), beta=5.0)
    root = _N("root", [_N("tip")])
    subst = {root: 0.0, root.children[0]: 6.0}
    tips = sel.evolve_family(root, subst, "A" * len(target), rng=np.random.default_rng(0))
    frac = np.mean([a == b for a, b in zip(tips["tip"], target)])
    assert frac > 0.9, f"only recovered {frac:.0%} of the injected preference"


def test_recovery_increases_with_beta():
    target = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQD"
    root = _N("root", [_N("tip")])
    subst = {root: 0.0, root.children[0]: 2.0}
    critic = FixedProfileCritic(_peaked_profile(target, hi=0.9))

    def recovered(beta, seed):
        sel = PLMSelection(critic, beta=beta)
        tips = sel.evolve_family(root, subst, "A" * len(target), rng=np.random.default_rng(seed))
        return np.mean([a == b for a, b in zip(tips["tip"], target)])

    means = [np.mean([recovered(b, s) for s in range(4)]) for b in (0.0, 1.0, 2.0, 4.0, 8.0)]
    assert all(means[i] <= means[i + 1] + 0.05 for i in range(len(means) - 1)), means  # ~monotone
    assert means[-1] > means[0] + 0.3, means                                          # strong effect


def test_frozen_is_deterministic_given_seed():
    sel = PLMSelection(FixedProfileCritic(_peaked_profile("ACDEFGHIKLMNPQRSTVWY")), beta=3.0)
    root = _N("root", [_N("a"), _N("b")])
    subst = {root: 0.0, root.children[0]: 1.0, root.children[1]: 1.0}
    a = sel.evolve_family(root, subst, "A" * 20, rng=np.random.default_rng(7))
    b = sel.evolve_family(root, subst, "A" * 20, rng=np.random.default_rng(7))
    assert a == b


def test_internal_node_is_recorded_and_children_inherit_it():
    # a 3-level tree exercises recursion through a genuine internal node + the split rule
    target = "MQIFVKTLTGKTITLEVEP"
    sel = PLMSelection(FixedProfileCritic(_peaked_profile(target, hi=0.97)), beta=5.0)
    t1, t2 = _N("t1"), _N("t2")
    internal = _N("internal", [t1, t2])
    root = _N("root", [internal])
    subst = {root: 0.0, internal: 6.0, t1: 0.0, t2: 0.0}     # zero-length tips must copy the internal
    out = sel.evolve_family(root, subst, "A" * len(target), rng=np.random.default_rng(0))
    assert "internal" in out
    assert out["t1"] == out["internal"] == out["t2"]
    assert np.mean([a == b for a, b in zip(out["internal"], target)]) > 0.9


def test_high_beta_and_delta_profile_stay_numerically_stable():
    # a hard-zero (delta) profile at very large beta must not produce NaN targets (log-space + floor)
    L = 10
    prof = np.zeros((L, 20))
    prof[:, _AA["W"]] = 1.0                                  # every site prefers W, zeros elsewhere
    for beta in (20.0, 200.0):
        sel = PLMSelection(FixedProfileCritic(prof), beta=beta)
        st = sel.site_targets("A" * L)
        assert np.all(np.isfinite(st)) and st.min() > 0
        root = _N("root", [_N("tip")])
        out = sel.evolve_family(root, {root: 0.0, root.children[0]: 6.0},
                                "A" * L, rng=np.random.default_rng(0))
        assert out["tip"] == "W" * L


# --------------------------------------------------------------------------- #
# batching + packaging
# --------------------------------------------------------------------------- #
def test_evolve_families_batches_and_is_reproducible():
    sel = PLMSelection(FixedProfileCritic(_peaked_profile("ACDEFG")), beta=2.0)
    r1, r2 = _N("r1", [_N("x")]), _N("r2", [_N("y")])
    node_trees = {
        "fam1": {"complete": (r1, {r1: 0.0, r1.children[0]: 1.0}), "extant": None},
        "fam2": {"complete": None, "extant": (r2, {r2: 0.0, r2.children[0]: 1.0})},
        "fam3": {"complete": None, "extant": None},                    # nothing to evolve
    }
    seqs = {"fam1": "AAAAAA", "fam2": "AAAAAA", "fam3": "AAAAAA"}
    out = sel.evolve_families(node_trees, seqs, seed=1)
    assert set(out) == {"fam1", "fam2"}
    assert "x" in out["fam1"] and "y" in out["fam2"]            # every node recorded (root + tip)
    assert sel.evolve_families(node_trees, seqs, seed=1) == out  # reproducible for a fixed seed


def test_selection_module_has_no_top_level_ml_imports():
    # lazy-import contract: torch/esm are imported ONLY inside ESM2Critic, never at module load, so
    # `import zombi2.experimental.selection` and the whole frozen path work without zombi2[selection].
    from zombi2.experimental import selection
    tree = ast.parse(inspect.getsource(selection))
    top: set[str] = set()
    for node in tree.body:
        if isinstance(node, ast.Import):
            top.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            top.add(node.module.split(".")[0])
    assert "torch" not in top and "esm" not in top, f"module-level ML import found: {top}"


# --------------------------------------------------------------------------- #
# the real ESM2 path (only when the optional deps are installed)
# --------------------------------------------------------------------------- #
def test_esm2critic_profile_vocab_order_and_score():
    pytest.importorskip("torch")
    pytest.importorskip("esm")
    from zombi2.experimental.selection import ESM2Critic
    critic = ESM2Critic("esm2_t6_8M_UR50D")
    ubi = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQD"
    prof = critic.profile(ubi)
    assert prof.shape == (len(ubi), 20)
    assert np.allclose(prof.sum(1), 1.0, atol=1e-5)
    # if the ESM2-vocab -> AMINO_ACIDS column mapping is correct, the top-predicted residue matches
    # the native residue at far more than chance (5%); a scrambled mapping would sit near chance.
    top = [AMINO_ACIDS[i] for i in prof.argmax(1)]
    native_match = np.mean([t == a for t, a in zip(top, ubi)])
    assert native_match > 0.3, f"vocab order likely wrong: only {native_match:.0%} native-argmax"
    s = critic.score([ubi, ubi[::-1]])
    assert s.shape == (2,) and np.all(np.isfinite(s))
