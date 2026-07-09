"""Unit tests for the RED benchmark's pure core (workflow/scripts/red.py).

These lock the two invariants the benchmark rests on — RED on a time tree equals relative time,
and a strict clock recovers ages exactly — plus the clock factory, the ultrametric guarantee, and
seed determinism. They need only ZOMBI2 (+ numpy), not Snakemake, so they run in CI cheaply.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts"))
import red  # noqa: E402


def _tree(n=120, seed=1):
    return red.simulate_time_tree("yule", n_tips=n, birth=1.0, death=0.0, seed=seed)


def test_red_on_timetree_equals_relative_time():
    # THE core identity: RED with branch lengths == time durations equals node.time / total_age.
    tree = _tree()
    bl = {n: n.branch_length() for n in tree.nodes_preorder()}
    r = red.compute_red(tree, bl)
    dev = max(abs(r[n] - n.time / tree.total_age) for n in tree.internal_nodes())
    assert dev < 1e-9, f"RED-on-timetree deviates from relative time by {dev:.2e}"


def test_red_endpoints():
    tree = _tree()
    bl = {n: n.branch_length() for n in tree.nodes_preorder()}
    r = red.compute_red(tree, bl)
    assert r[tree.root] == 0.0
    assert all(abs(r[lf] - 1.0) < 1e-12 for lf in tree.leaves())


def test_strict_clock_recovers_ages_exactly():
    import zombi2 as z
    tree = _tree()
    scaled = z.StrictClock().scale(tree, seed=3)
    _, m = red.red_recovery(tree, scaled, tree.total_age)
    assert m["pearson_r"] > 0.9999
    assert m["nrmse"] < 1e-6
    assert m["fold_range"] == pytest.approx(1.0)      # a strict clock has no rate spread


def test_ratevar_degrades_but_stays_accurate():
    # Under ~30x bounded rate variation RED still recovers ages well (the Rinke/GTDB result).
    tree = _tree(n=200)
    spec = {"clock": "ratevar", "switch_rate": 1.0, "spread": 5.5, "n_bins": 15}
    rs = [red.red_recovery(tree, red.build_clock(spec).scale(tree, seed=s), tree.total_age)[1]
          for s in range(5)]
    assert np.mean([m["pearson_r"] for m in rs]) > 0.95
    assert np.mean([m["fold_range"] for m in rs]) > 1.0


@pytest.mark.parametrize("spec", [
    {"clock": "strict"},
    {"clock": "ratevar", "switch_rate": 0.5, "spread": 5.5, "n_bins": 15},
    {"clock": "aln", "sigma": 0.5},
    {"clock": "cir", "theta": 1.0, "sigma": 0.5},
    {"clock": "whitenoise", "sigma": 0.4},
    {"clock": "ucln", "sigma": 0.4},
])
def test_build_clock_variants_scale(spec):
    tree = _tree(n=60)
    scaled = red.build_clock(spec).scale(tree, seed=7)
    assert len(scaled.branch_lengths) == len(list(tree.nodes_preorder()))
    _, m = red.red_recovery(tree, scaled, tree.total_age)
    assert 0.0 <= m["pearson_r"] <= 1.0 + 1e-9


def test_build_clock_rejects_unknown():
    with pytest.raises(ValueError):
        red.build_clock({"clock": "nope"})


def test_bd_tree_is_ultrametric():
    tree = red.simulate_time_tree("bd", n_tips=80, birth=1.0, death=0.4, seed=2)
    tips = [lf.time for lf in tree.leaves()]
    assert max(tips) - min(tips) < 1e-6 * tree.total_age
    assert len(tree.leaves()) == 80


def test_bd_tree_calibration_is_correct():
    # The returned bd tree must satisfy the root=0 / total_age==depth contract, so RED on it
    # recovers ages exactly under a strict clock — directly on the object, no Newick round-trip.
    import zombi2 as z
    tree = red.simulate_time_tree("bd", n_tips=100, birth=1.0, death=0.4, seed=1)
    assert tree.root.time == pytest.approx(0.0)
    assert tree.total_age == pytest.approx(max(lf.time for lf in tree.leaves()))
    _, m = red.red_recovery(tree, z.StrictClock().scale(tree, seed=1), tree.total_age)
    assert m["nrmse"] < 1e-6      # calibration, not just correlation


def test_derive_seed_deterministic_and_distinct():
    assert red.derive_seed("a", 1, "x", 0) == red.derive_seed("a", 1, "x", 0)
    assert red.derive_seed("a", 1, "x", 0) != red.derive_seed("a", 1, "x", 1)
