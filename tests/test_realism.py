"""P4 realism metric: Fréchet-ESM distance. Contract tests are torch-free; the real discrimination
test runs only with the optional ESM2 deps installed."""
from __future__ import annotations

import ast
import inspect

import numpy as np
import pytest

import zombi2.experimental as ex
from zombi2.experimental.realism import frechet_esm_distance
from zombi2.experimental.selection import FixedProfileCritic


def test_frechet_needs_a_critic_that_embeds():
    # FixedProfileCritic has no embedding space -> a clear NotImplementedError, not a crash
    crit = FixedProfileCritic(np.full((4, 20), 1 / 20.0))
    with pytest.raises(NotImplementedError, match="embed"):
        frechet_esm_distance(["ACDE", "ACDE"], ["ACDE", "ACDE"], crit)


def test_exports_stay_in_the_experimental_namespace():
    import zombi2
    assert "frechet_esm_distance" in ex.__all__ and hasattr(ex, "frechet_esm_distance")
    assert not hasattr(zombi2, "frechet_esm_distance")


def test_module_has_no_top_level_ml_imports():
    from zombi2.experimental import realism
    tree = ast.parse(inspect.getsource(realism))
    top: set[str] = set()
    for node in tree.body:
        if isinstance(node, ast.Import):
            top.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            top.add(node.module.split(".")[0])
    assert "torch" not in top and "esm" not in top and "scipy" not in top, top


class _StubCritic:
    """Duck-typed critic: maps each sequence (used as a lookup key) to a fixed embedding row, so the
    covariance/sqrtm machinery can be tested without torch/esm."""

    def __init__(self, table):
        self.table = table

    def embed(self, seqs):
        return np.array([self.table[s] for s in seqs], dtype=float)


def test_frechet_is_zero_for_a_set_versus_itself():
    # the headline contract: FID(X, X) == 0. Exercises the sqrtm(Sa @ Sa) = Sa path that the
    # discrimination test never touches (there dmu dominates).
    rng = np.random.default_rng(0)
    n, dim = 6, 8
    rows = rng.normal(size=(n, dim))
    crit = _StubCritic({f"s{i}": rows[i] for i in range(n)})
    X = [f"s{i}" for i in range(n)]
    assert frechet_esm_distance(X, X, crit) == pytest.approx(0.0, abs=1e-6)


def test_frechet_covariance_term_drives_the_distance_when_means_match():
    # two sets with the SAME mean but different spread -> the ||dmu||^2 term is ~0, so a nonzero
    # distance can only come from the covariance/sqrtm term. A mean-only metric would score ~0 here.
    rng = np.random.default_rng(1)
    n, dim = 40, 6
    base = rng.normal(size=(n, dim))
    base -= base.mean(0)                                          # centre -> zero mean
    A, B = base, 3.0 * base                                       # same zero mean, 9x the covariance
    table = {}
    ak = [f"a{i}" for i in range(n)]
    bk = [f"b{i}" for i in range(n)]
    for i in range(n):
        table[ak[i]], table[bk[i]] = A[i], B[i]
    crit = _StubCritic(table)
    dmu = A.mean(0) - B.mean(0)
    assert float(dmu @ dmu) < 1e-6                                # means match -> mean-only distance ~0
    assert frechet_esm_distance(ak, bk, crit) > 5.0              # ... so the covariance term drives it


def test_frechet_esm_discriminates_same_family_from_different_protein():
    pytest.importorskip("torch")
    pytest.importorskip("esm")
    from zombi2.experimental.selection import ESM2Critic
    critic = ESM2Critic("esm2_t6_8M_UR50D")
    ubi = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
    other = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR"
    rng = np.random.default_rng(0)
    aa = "ACDEFGHIKLMNPQRSTVWY"

    def variants(seq, n=12):
        out = []
        for _ in range(n):
            s = list(seq)
            for _ in range(3):                                   # a few point mutations each
                s[rng.integers(len(s))] = aa[rng.integers(20)]
            out.append("".join(s))
        return out

    a1, a2, b = variants(ubi), variants(ubi), variants(other)
    d_self = frechet_esm_distance(a1, a2, critic)                # two samples of the same family -> small
    d_diff = frechet_esm_distance(a1, b, critic)                 # a different protein -> larger
    assert np.isfinite(d_self) and np.isfinite(d_diff)
    assert d_diff > d_self                                       # the metric separates real classes
