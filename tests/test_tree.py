"""The zombi2.tree toolkit — transforms (with_stem, make_ultrametric, rescale) and analyses
(relative_evolutionary_divergence, red_scaled, distance). read_newick / prune are covered in
test_cli.py."""
import numpy as np
import pytest

from zombi2 import species
from zombi2 import tree as T


def _tree(n=30, seed=1, death=0.0):
    return species.simulate_species_tree(birth=1.0, death=death, n_extant=n, seed=seed).extant_tree


def _depths(t):
    return T._depths(t)


def _stem(t):
    r = t.nodes[t.root]
    return r.end_time - r.birth_time


# ── with_stem ──────────────────────────────────────────────────────────────────────────
def test_with_stem_set_and_add_leave_the_rest_untouched():
    t = _tree()
    s0 = _stem(t)
    assert _stem(T.with_stem(t, 5.0)) == pytest.approx(5.0)
    assert _stem(T.with_stem(t, 2.0, mode="add")) == pytest.approx(s0 + 2.0)
    # a transform returns a copy; the input is unchanged, and every non-root branch is preserved
    out = T.with_stem(t, 5.0)
    assert _stem(t) == pytest.approx(s0)
    for i, n in t.nodes.items():
        if n.parent is not None:
            o = out.nodes[i]
            assert (o.end_time - o.birth_time) == pytest.approx(n.end_time - n.birth_time)


def test_with_stem_rejects_bad_mode():
    with pytest.raises(ValueError):
        T.with_stem(_tree(), 1.0, mode="grow")


# ── make_ultrametric ───────────────────────────────────────────────────────────────────
def test_make_ultrametric_snaps_rounding_noise():
    t = _tree()
    # perturb one tip by a rounding-scale amount
    leaf = next(i for i, n in t.nodes.items() if n.children is None)
    t.nodes[leaf].end_time += 1e-4 * max(_depths(t).values())
    u = T.make_ultrametric(t, tol=1e-2)
    d = [_depths(u)[i] for i, n in u.nodes.items() if n.children is None]
    assert max(d) - min(d) == pytest.approx(0.0, abs=1e-9)


def test_make_ultrametric_refuses_a_real_gap():
    t, _ = T.read_newick("((a:1,b:1):1,c:1.5);",
                         tip_fates={"a": "extant", "b": "extant", "c": "extant"})
    with pytest.raises(ValueError, match="not contemporaneous"):
        T.make_ultrametric(t, tol=1e-3)


# ── rescale ────────────────────────────────────────────────────────────────────────────
def test_rescale_to_height_and_by_factor():
    t = _tree()
    to1 = T.rescale(t, height=1.0)
    d = [_depths(to1)[i] for i, n in to1.nodes.items() if n.children is None]
    assert max(d) == pytest.approx(1.0)
    half = T.rescale(t, factor=0.5)
    assert _stem(half) == pytest.approx(_stem(t) * 0.5)


def test_rescale_needs_exactly_one_of_height_factor():
    t = _tree()
    with pytest.raises(ValueError):
        T.rescale(t)
    with pytest.raises(ValueError):
        T.rescale(t, height=1.0, factor=2.0)


# ── RED ────────────────────────────────────────────────────────────────────────────────
def test_red_is_exact_relative_age_on_an_ultrametric_tree():
    t = _tree(n=40, seed=7)
    red = T.relative_evolutionary_divergence(t)
    root, H = t.nodes[t.root], max(n.end_time for n in t.nodes.values() if n.children is None)
    for i, n in t.nodes.items():
        if n.children is not None and n.parent is not None:
            assert red[i] == pytest.approx((n.end_time - root.end_time) / (H - root.end_time), abs=1e-9)
    assert red[t.root] == 0.0
    assert all(red[i] == pytest.approx(1.0) for i, n in t.nodes.items() if n.children is None)


def test_red_scaled_is_ultrametric_on_unit_interval():
    rs = T.red_scaled(_tree(n=25))
    d = _depths(rs)
    tips = [d[i] for i, n in rs.nodes.items() if n.children is None]
    assert all(x == pytest.approx(1.0) for x in tips)


# ── distance ───────────────────────────────────────────────────────────────────────────
def test_distance_of_a_tree_with_itself_is_zero():
    t = _tree()
    assert T.distance(t, t, metric="rf") == 0.0
    assert T.distance(t, t, metric="rf-normalized") == 0.0
    assert T.distance(t, t, metric="branch-score") == pytest.approx(0.0)


def test_distance_detects_a_topology_difference():
    # same taxa (ids), different topology: build two 4-tip trees over ids {1,2,3,4}
    from zombi2.tree import Node, Tree
    a = Tree({0: Node(0, None, 0.0, 1.0, (5, 6), "speciation"),
              5: Node(5, 0, 1.0, 2.0, (1, 2), "speciation"),
              6: Node(6, 0, 1.0, 2.0, (3, 4), "speciation"),
              1: Node(1, 5, 2.0, 3.0, None, "extant"), 2: Node(2, 5, 2.0, 3.0, None, "extant"),
              3: Node(3, 6, 2.0, 3.0, None, "extant"), 4: Node(4, 6, 2.0, 3.0, None, "extant")}, 0)
    b = Tree({0: Node(0, None, 0.0, 1.0, (5, 6), "speciation"),
              5: Node(5, 0, 1.0, 2.0, (1, 3), "speciation"),          # 1+3 vs 1+2
              6: Node(6, 0, 1.0, 2.0, (2, 4), "speciation"),
              1: Node(1, 5, 2.0, 3.0, None, "extant"), 3: Node(3, 5, 2.0, 3.0, None, "extant"),
              2: Node(2, 6, 2.0, 3.0, None, "extant"), 4: Node(4, 6, 2.0, 3.0, None, "extant")}, 0)
    assert T.distance(a, b, metric="rf") == 4.0        # both non-trivial clades differ, symmetric
    assert 0.0 < T.distance(a, b, metric="rf-normalized") <= 1.0


def test_distance_raises_on_different_leaf_sets():
    with pytest.raises(ValueError, match="different leaf sets"):
        T.distance(_tree(n=10, seed=1), _tree(n=12, seed=2))
