"""Tests for Relative Evolutionary Divergence (``zombi2.tools.red``)."""
from __future__ import annotations

import pytest

from zombi2 import BirthDeath, StrictClock, RateVariation
from zombi2.species.forward import simulate_forward
from zombi2.tree import read_newick
from zombi2.tools import relative_evolutionary_divergence as red_of


def by_name(tree, red):
    return {n.name: red[n] for n in tree.nodes_preorder()}


def test_root_zero_and_leaves_one():
    tree = read_newick("((A:1,B:1)X:2,C:3)R;")
    red = red_of(tree)
    assert red[tree.root] == 0.0
    for leaf in tree.leaves():
        assert red[leaf] == pytest.approx(1.0)


def test_hand_computed_small_tree():
    # ((A:1,B:1)X:2,C:3)R.  mean-tip-dist(X)=1, so RED(X)=2/(2+1)=2/3; leaves=1; root=0.
    tree = read_newick("((A:1,B:1)X:2,C:3)R;")
    red = by_name(tree, red_of(tree))
    assert red["R"] == 0.0
    assert red["X"] == pytest.approx(2.0 / 3.0)
    assert red["A"] == pytest.approx(1.0)
    assert red["B"] == pytest.approx(1.0)
    assert red["C"] == pytest.approx(1.0)


def test_ultrametric_red_equals_relative_age():
    # On an ultrametric (Yule) tree RED(node) is exactly node.time / total_age.
    tree = simulate_forward(BirthDeath(birth=1.0, death=0.0), n_tips=60, seed=7)
    red = red_of(tree)
    for node in tree.internal_nodes():
        assert red[node] == pytest.approx(node.time / tree.total_age, abs=1e-9)
    for leaf in tree.leaves():
        assert red[leaf] == pytest.approx(1.0)


def test_scale_invariance():
    # Multiplying every branch length by a constant must not change RED.
    tree = simulate_forward(BirthDeath(birth=1.0, death=0.0), n_tips=40, seed=3)
    base = red_of(tree)
    scaled = red_of(tree, branch_length=lambda n: 12.5 * n.branch_length())
    for node in tree.nodes_preorder():
        assert scaled[node] == pytest.approx(base[node], abs=1e-12)


def test_monotone_root_to_tip():
    tree = simulate_forward(BirthDeath(birth=1.0, death=0.0), n_tips=50, seed=11)
    red = red_of(tree)
    for node in tree.nodes_preorder():
        if node.parent is not None:
            assert red[node] >= red[node.parent] - 1e-12
            if node.branch_length() > 0:
                assert red[node] > red[node.parent]


def test_rate_scaled_tree_input_strict_clock():
    # A strict clock rescales every branch by one constant -> phylogram is the timetree x rate,
    # so RED on it must still equal relative age. Exercises the RateScaledTree input path.
    tree = simulate_forward(BirthDeath(birth=1.0, death=0.0), n_tips=45, seed=5)
    phylogram = StrictClock(rate=3.7).scale(tree, seed=0)
    red = red_of(phylogram)
    for node in tree.internal_nodes():
        assert red[node] == pytest.approx(node.time / tree.total_age, abs=1e-9)


def test_rate_variation_phylogram_runs_and_is_bounded():
    # A real rate-varying phylogram: RED stays a valid [0,1] scale (root 0, tips 1, monotone).
    tree = simulate_forward(BirthDeath(birth=1.0, death=0.0), n_tips=80, seed=9)
    phylogram = RateVariation(bins=[0.25, 0.5, 1.0, 2.0, 4.0], switch_rate=0.5).scale(tree, seed=1)
    red = red_of(phylogram)
    assert red[tree.root] == 0.0
    assert all(red[lf] == pytest.approx(1.0) for lf in tree.leaves())
    assert all(0.0 <= red[n] <= 1.0 + 1e-12 for n in tree.nodes_preorder())


def test_rate_scaled_tree_with_accessor_raises():
    tree = simulate_forward(BirthDeath(birth=1.0, death=0.0), n_tips=10, seed=1)
    phylogram = StrictClock(rate=1.0).scale(tree, seed=0)
    with pytest.raises(ValueError, match="not both"):
        red_of(phylogram, branch_length=lambda n: 1.0)


def test_negative_branch_length_raises():
    tree = read_newick("((A:1,B:1)X:2,C:3)R;")
    with pytest.raises(ValueError, match="negative branch length"):
        red_of(tree, branch_length=lambda n: -1.0 if n.parent is not None else 0.0)


def test_zero_length_branch_inherits_parent():
    # A zero-length branch (a+b could still be >0 via subtree) — check the a+b==0 fallback too:
    # a single internal node whose only child is a zero-length leaf inherits the parent RED.
    tree = read_newick("(A:0)R;")
    red = red_of(tree)
    # A sits on a zero-length branch and R is the root (RED 0); with a+b==0, A inherits 0.
    assert red[tree.root] == 0.0
    a = next(n for n in tree.leaves())
    assert red[a] == pytest.approx(0.0)


def test_cli_red_writes_tsv(tmp_path):
    from zombi2.cli import main
    nwk = tmp_path / "tree.nwk"
    nwk.write_text("((A:1,B:1)X:2,C:3)R;\n")
    out = tmp_path / "out"
    rc = main(["tools", "red", "-t", str(nwk), "-o", str(out)])
    assert rc == 0
    text = (out / "RED.tsv").read_text().strip().splitlines()
    header = text[0].split("\t")
    assert header == ["node", "is_leaf", "red"]
    rows = {r.split("\t")[0]: r.split("\t") for r in text[1:]}
    assert float(rows["R"][2]) == pytest.approx(0.0)
    assert float(rows["X"][2]) == pytest.approx(2.0 / 3.0)
    assert float(rows["A"][2]) == pytest.approx(1.0)
    assert rows["X"][1] == "False" and rows["A"][1] == "True"
