"""Tests for tree distances (``zombi2.tools.treedist``): RF, branch score, quartet."""
from __future__ import annotations

import math

import pytest

from zombi2 import BirthDeath
from zombi2.species.forward import simulate_forward
from zombi2.tools.treedist import (
    robinson_foulds,
    branch_score,
    quartet_distance,
    matching_distance,
    compare_trees,
)

# Two four-leaf trees with the same leaves but incompatible resolutions.
ABCD_1 = "((A,B),(C,D));"
ABCD_2 = "((A,C),(B,D));"


# --------------------------------------------------------------------------- identity

def test_identical_trees_are_zero():
    t = "((A:1,B:1):1,(C:1,D:2):1);"
    assert robinson_foulds(t, t).rf == 0
    assert robinson_foulds(t, t, rooted=False).rf == 0
    assert branch_score(t, t) == 0.0
    assert quartet_distance(t, t).differing == 0


def test_self_compare_on_simulated_tree():
    tree = simulate_forward(BirthDeath(birth=1.0, death=0.0), n_tips=30, seed=7)
    c = compare_trees(tree, tree)
    assert c.rf == 0 and c.rf_unrooted == 0
    assert c.branch_score == pytest.approx(0.0)
    assert c.quartet == 0 and c.quartet_normalized == 0.0


# --------------------------------------------------------------------------- Robinson–Foulds

def test_rf_hand_computed_incompatible_quartet():
    # T1 clusters {A,B},{C,D}; T2 clusters {A,C},{B,D} — fully disjoint.
    r = robinson_foulds(ABCD_1, ABCD_2)               # rooted
    assert (r.rf, r.max_rf, r.normalized) == (4, 4, 1.0)
    u = robinson_foulds(ABCD_1, ABCD_2, rooted=False)  # unrooted
    assert (u.rf, u.max_rf, u.normalized) == (2, 2, 1.0)


def test_rooted_differs_but_unrooted_agrees_under_reroot():
    # Same unrooted topology, different root placement: unrooted RF is 0, rooted RF is not.
    balanced = "((A,B),(C,D));"
    ladder = "(A,(B,(C,D)));"
    assert robinson_foulds(balanced, ladder, rooted=False).rf == 0
    assert robinson_foulds(balanced, ladder, rooted=True).rf == 2


def test_rf_is_symmetric():
    a = robinson_foulds(ABCD_1, ABCD_2).rf
    b = robinson_foulds(ABCD_2, ABCD_1).rf
    assert a == b


# --------------------------------------------------------------------------- branch score

def test_branch_score_hand_computed():
    # Differences live on {D} (1 vs 2) and {A,B} (1 vs 2); everything else matches.
    t1 = "((A:1,B:1):1,(C:1,D:1):1);"
    t2 = "((A:1,B:1):2,(C:1,D:2):1);"
    assert branch_score(t1, t2, order=1) == pytest.approx(2.0)
    assert branch_score(t1, t2, order=2) == pytest.approx(math.sqrt(2.0))


def test_branch_score_is_symmetric_and_nonnegative():
    t1 = "((A:1,B:2):1,(C:3,D:1):2);"
    t2 = "((A:2,B:1):3,(C:1,D:1):1);"
    assert branch_score(t1, t2) == pytest.approx(branch_score(t2, t1))
    assert branch_score(t1, t2) >= 0.0


def test_branch_score_rejects_bad_order():
    with pytest.raises(ValueError):
        branch_score(ABCD_1, ABCD_1, order=3)


# --------------------------------------------------------------------------- quartet

def test_quartet_single_incompatible_quartet():
    q = quartet_distance(ABCD_1, ABCD_2)
    assert (q.differing, q.total, q.normalized) == (1, 1, 1.0)


def test_quartet_caterpillar_topology():
    # (((A,B),C),D) induces the unrooted split AB|CD, same as ((A,B),(C,D)).
    caterpillar = "(((A,B),C),D);"
    assert quartet_distance(caterpillar, ABCD_1).differing == 0
    # but AC|BD disagrees with the caterpillar's AB|CD.
    assert quartet_distance(caterpillar, ABCD_2).differing == 1


def test_quartet_unresolved_star_differs_from_resolved():
    star = "(A,B,C,D);"
    q = quartet_distance(star, ABCD_1)
    assert q.differing == 1                    # star is unresolved, ABCD_1 is AB|CD
    assert quartet_distance(star, star).differing == 0


def test_quartet_guard_on_large_tree():
    tree = simulate_forward(BirthDeath(birth=1.0, death=0.0), n_tips=20, seed=1)
    with pytest.raises(ValueError, match="max_leaves"):
        quartet_distance(tree, tree, max_leaves=10)
    # explicit override runs it
    assert quartet_distance(tree, tree, max_leaves=20).differing == 0


# --------------------------------------------------------------------------- errors & aggregate

def test_mismatched_leaf_sets_raise():
    with pytest.raises(ValueError, match="same leaf set"):
        robinson_foulds("((A,B),C);", "((A,B),D);")


def test_compare_trees_skips_quartet_over_max_leaves():
    tree = simulate_forward(BirthDeath(birth=1.0, death=0.0), n_tips=15, seed=2)
    c = compare_trees(tree, tree, max_leaves=10)   # 15 > 10 ⇒ quartet skipped, not errored
    assert c.quartet is None and c.quartet_normalized is None
    assert c.rf == 0 and c.branch_score == pytest.approx(0.0)


def test_compare_trees_reports_all_metrics():
    c = compare_trees(ABCD_1, ABCD_2)
    assert c.n_leaves == 4
    assert c.rf == 4 and c.rf_unrooted == 2
    assert c.quartet == 1 and c.quartet_normalized == 1.0
    assert c.branch_score == pytest.approx(0.0)    # no branch lengths ⇒ all zero
    assert c.matching == 4 and c.matching_normalized == pytest.approx(0.5)  # rooted (clusters)


# --------------------------------------------------------------------------- matching distance

def test_matching_identical_is_zero():
    assert matching_distance(ABCD_1, ABCD_1).distance == 0
    assert matching_distance(ABCD_1, ABCD_1, rooted=False).distance == 0


def test_matching_hand_computed_incompatible_quartet():
    # rooted clusters {A,B},{C,D} vs {A,C},{B,D}: every pairing costs 2 leaf-moves ⇒ 4.
    r = matching_distance(ABCD_1, ABCD_2, rooted=True)
    assert (r.distance, r.max_distance, r.normalized) == (4, 8, 0.5)
    # unrooted: one split each ({C,D} vs {B,D}), costing 2.
    u = matching_distance(ABCD_1, ABCD_2, rooted=False)
    assert (u.distance, u.max_distance, u.normalized) == (2, 4, 0.5)


def test_matching_padding_against_star():
    # A star has no informative clusters, so both of the resolved tree's clusters go unmatched
    # (each paying its null weight 2) ⇒ distance 4 = max ⇒ normalized 1.
    star = "(A,B,C,D);"
    r = matching_distance(ABCD_1, star, rooted=True)
    assert (r.distance, r.normalized) == (4, 1.0)


def test_matching_is_symmetric():
    assert matching_distance(ABCD_1, ABCD_2).distance == matching_distance(ABCD_2, ABCD_1).distance


def test_matching_guard_on_large_tree():
    tree = simulate_forward(BirthDeath(birth=1.0, death=0.0), n_tips=20, seed=1)
    with pytest.raises(ValueError, match="max_leaves"):
        matching_distance(tree, tree, max_leaves=10)
    assert matching_distance(tree, tree, max_leaves=20).distance == 0    # a tree vs itself
