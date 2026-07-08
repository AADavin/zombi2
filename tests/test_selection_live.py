"""P3 live/epistatic mode: PLMSelection(mode="live") re-reads the critic on the current sequence
every ``refresh`` substitutions/site, so sites feel each other's changes. Torch-free (FixedProfile /
tiny custom critics)."""
from __future__ import annotations

import numpy as np
import pytest

from zombi2.experimental.selection import Critic, FixedProfileCritic, PLMSelection
from zombi2.sequences.models import AMINO_ACIDS

_AA = {a: i for i, a in enumerate(AMINO_ACIDS)}


def _peaked(target: str, hi: float = 0.95) -> np.ndarray:
    p = np.full((len(target), 20), (1.0 - hi) / 19.0)
    for i, a in enumerate(target):
        p[i, _AA[a]] = hi
    return p


class _N:
    def __init__(self, gid, children=()):
        self.gid = gid
        self.children = list(children)


class _CountingCritic(Critic):
    """Returns a fixed profile but counts how many times it is queried."""

    def __init__(self, profile):
        self.n = 0
        self._p = profile

    def profile(self, seq):
        self.n += 1
        return self._p


class _CoupledCritic(Critic):
    """Epistatic: site 0 always wants W; every other site wants site 0's CURRENT residue. Under live
    the whole sequence follows site 0 to W; under frozen the followers are pinned to the root's site-0
    residue and never move."""

    def profile(self, seq):
        p = np.full((len(seq), 20), 1e-6)
        p[0, _AA["W"]] = 1.0                       # site 0: always W (context-free)
        s0 = _AA[seq[0]]                            # site 0's *current* residue
        for i in range(1, len(seq)):
            p[i, s0] = 1.0
        return p / p.sum(1, keepdims=True)


# --------------------------------------------------------------------------- #
# construction / validation
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), 1e-12, "fast"])
def test_live_rejects_bad_refresh(bad):
    with pytest.raises(ValueError, match="refresh"):
        PLMSelection(FixedProfileCritic(_peaked("ACDE")), mode="live", refresh=bad)


def test_live_constructs_and_infinite_refresh_is_allowed():
    PLMSelection(FixedProfileCritic(_peaked("ACDE")), mode="live", refresh=float("inf"))


# --------------------------------------------------------------------------- #
# the live == frozen limit, and the refresh knob
# --------------------------------------------------------------------------- #
def test_live_with_infinite_refresh_equals_frozen():
    critic = FixedProfileCritic(_peaked("MQIFVKTLTG"))
    root = _N("root", [_N("a"), _N("b")])
    subst = {root: 0.0, root.children[0]: 2.0, root.children[1]: 2.0}
    frozen = PLMSelection(critic, mode="frozen", beta=3.0).evolve_family(
        root, subst, "A" * 10, rng=np.random.default_rng(1))
    live = PLMSelection(critic, mode="live", beta=3.0, refresh=float("inf")).evolve_family(
        root, subst, "A" * 10, rng=np.random.default_rng(1))
    assert frozen == live                           # never refreshing == frozen, exactly


def test_infinite_refresh_equals_frozen_with_ambiguity_and_a_content_critic():
    # a content-dependent critic + an ambiguity code ('X') at the root used to diverge frozen vs
    # live(inf) (frozen scored the raw 'X', live the resampled residue). Both now build from the
    # concrete decoded root, so they match exactly.
    critic = _CoupledCritic()
    root = _N("root", [_N("tip")])
    subst = {root: 0.0, root.children[0]: 3.0}
    frozen = PLMSelection(critic, mode="frozen", beta=6.0).evolve_family(
        root, subst, "XAAAAAAA", rng=np.random.default_rng(1))
    live = PLMSelection(critic, mode="live", beta=6.0, refresh=float("inf")).evolve_family(
        root, subst, "XAAAAAAA", rng=np.random.default_rng(1))
    assert frozen == live


def test_live_accrued_refresh_carries_across_a_branch_boundary():
    # refresh=0.5; the internal branch (0.3) leaves accrued=0.3, so each 0.4 tip crosses 0.5 exactly
    # once => 3 critic reads (root + one per tip). If accrued reset at the node it would be just 1.
    c = _CountingCritic(_peaked("MQIFVKTL"))
    t1, t2 = _N("t1"), _N("t2")
    root = _N("root", [_N("internal", [t1, t2])])
    internal = root.children[0]
    subst = {root: 0.0, internal: 0.3, t1: 0.4, t2: 0.4}
    PLMSelection(c, mode="live", beta=2.0, refresh=0.5).evolve_family(
        root, subst, "A" * 8, rng=np.random.default_rng(0))
    assert c.n == 3, c.n


def test_live_refreshes_more_often_with_a_smaller_interval():
    p = _peaked("MQIFVKTL")
    root = _N("root", [_N("tip")])
    subst = {root: 0.0, root.children[0]: 4.0}

    def calls(mode, refresh=0.25):
        c = _CountingCritic(p)
        PLMSelection(c, mode=mode, beta=2.0, refresh=refresh).evolve_family(
            root, subst, "A" * 8, rng=np.random.default_rng(0))
        return c.n

    assert calls("frozen") == 1                      # frozen reads the critic exactly once
    n_coarse = calls("live", refresh=2.0)
    n_fine = calls("live", refresh=0.5)
    assert n_fine > n_coarse > 1                      # live re-reads; finer interval -> more reads


# --------------------------------------------------------------------------- #
# behaviour: recovery, determinism, and the epistasis payoff
# --------------------------------------------------------------------------- #
def test_live_recovers_injected_preference():
    target = "MQIFVKTLTGKTITLE"
    sel = PLMSelection(FixedProfileCritic(_peaked(target, hi=0.95)),
                       mode="live", beta=5.0, refresh=0.5)
    root = _N("root", [_N("tip")])
    subst = {root: 0.0, root.children[0]: 8.0}
    tip = sel.evolve_family(root, subst, "A" * len(target), rng=np.random.default_rng(0))["tip"]
    assert np.mean([a == b for a, b in zip(tip, target)]) > 0.9


def test_live_is_deterministic_given_seed():
    sel = PLMSelection(FixedProfileCritic(_peaked("ACDEFGHIKL")), mode="live", beta=2.0, refresh=0.5)
    root = _N("root", [_N("t")])
    subst = {root: 0.0, root.children[0]: 3.0}
    a = sel.evolve_family(root, subst, "A" * 10, rng=np.random.default_rng(4))
    b = sel.evolve_family(root, subst, "A" * 10, rng=np.random.default_rng(4))
    assert a == b


def test_live_captures_epistasis_that_frozen_misses():
    # site 0 -> W; under LIVE the followers re-read and chase site 0 to W; under FROZEN they are pinned
    # to the root's site-0 residue (A) forever. This is the whole point of the live mode.
    critic = _CoupledCritic()
    L = 8
    root = _N("root", [_N("tip")])
    subst = {root: 0.0, root.children[0]: 50.0}
    frozen = PLMSelection(critic, mode="frozen", beta=6.0).evolve_family(
        root, subst, "A" * L, rng=np.random.default_rng(0))["tip"]
    live = PLMSelection(critic, mode="live", beta=6.0, refresh=0.5).evolve_family(
        root, subst, "A" * L, rng=np.random.default_rng(0))["tip"]
    assert live.count("W") >= L - 1                  # live: the whole sequence follows site 0 to W
    assert frozen.count("W") <= 2                     # frozen: only site 0 moves
    assert frozen[1:].count("A") >= L - 3             # frozen followers stay at the root residue A


def test_live_works_through_evolve_families():
    sel = PLMSelection(FixedProfileCritic(_peaked("ACDEFG")), mode="live", beta=2.0, refresh=0.5)
    r = _N("r", [_N("x")])
    nt = {"f": {"complete": (r, {r: 0.0, r.children[0]: 2.0}), "extant": None}}
    out = sel.evolve_families(nt, {"f": "A" * 6}, seed=1)
    assert "x" in out["f"]
    assert sel.evolve_families(nt, {"f": "A" * 6}, seed=1) == out    # reproducible for a fixed seed
