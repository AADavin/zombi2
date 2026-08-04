"""The zombi2.tree toolkit — transforms (with_stem, make_ultrametric, rescale) and analyses
(relative_evolutionary_divergence, red_scaled, distance). read_newick / prune are covered in
test_cli.py."""
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


# ── gamma ──────────────────────────────────────────────────────────────────────────────
def test_gamma_is_standard_normal_under_constant_rate_pure_birth():
    """The statistic's whole use is that it has a known null: mean 0, sd 1 on a Yule tree. If that
    drifts, every comparison anyone makes against it is off by the drift."""
    import statistics

    from zombi2.species import simulate_species_tree

    gs = [T.gamma_statistic(simulate_species_tree(birth=1.0, n_extant=100, seed=s).extant_tree)
          for s in range(1, 201)]
    assert statistics.mean(gs) == pytest.approx(0.0, abs=0.25)
    assert statistics.stdev(gs) == pytest.approx(1.0, abs=0.2)


def test_gamma_goes_negative_when_speciation_slows_toward_the_present():
    """The reason to compute it at all: diversity dependence bunches the branching times early."""
    import statistics

    from zombi2.rates import modifiers as mod
    from zombi2.species import simulate_species_tree

    gs = [T.gamma_statistic(simulate_species_tree(birth=1.0 * mod.OnTotalDiversity(cap=110),
                                                  n_extant=100, seed=s).extant_tree)
          for s in range(1, 51)]
    assert statistics.mean(gs) < -3.0


def test_gamma_refuses_a_tree_it_is_not_defined_on():
    from zombi2.species import simulate_species_tree

    with pytest.raises(ValueError, match="ultrametric"):      # extinct lineages left in
        T.gamma_statistic(simulate_species_tree(birth=1.0, death=0.5, n_extant=30,
                                                seed=3).complete_tree)
    with pytest.raises(ValueError, match="at least 4 tips"):
        T.gamma_statistic(simulate_species_tree(birth=1.0, n_extant=3, seed=1).extant_tree)


def test_cli_tree_gamma_prints_one_number(tmp_path, capsys):
    from zombi2.cli.main import main
    from zombi2.species import simulate_species_tree

    t = simulate_species_tree(birth=1.0, n_extant=60, seed=7).extant_tree
    f = str(tmp_path / "t.nwk")
    open(f, "w").write(t.to_newick())
    assert main(["tools", "tree", f, "--gamma"]) == 0
    name, value = capsys.readouterr().out.strip().split("\t")
    assert name == "gamma"
    assert float(value) == pytest.approx(T.gamma_statistic(t), rel=1e-4)


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


# ── CLI: zombi2 tools tree / treedist ────────────────────────────────────────────────────
from zombi2.cli.main import main  # noqa: E402


def _write(tmp_path, name, newick):
    p = tmp_path / name
    p.write_text(newick, encoding="utf-8")
    return str(p)


def test_cli_tree_prune_writes_the_extant_tree_to_stdout(tmp_path, capsys):
    r = species.simulate_species_tree(birth=1.0, death=0.5, n_extant=8, seed=2)
    f = _write(tmp_path, "complete.nwk", r.complete_tree.to_newick())
    assert main(["tools", "tree", f, "--prune"]) == 0
    out = capsys.readouterr().out
    n_tips = out.count(":") - out.count(")")  # rough, but extant < complete tip count
    assert out.strip().endswith(";") and 0 < n_tips


def test_cli_tree_round_makes_a_noisy_tree_ultrametric(tmp_path, capsys):
    f = _write(tmp_path, "dated.nwk", "((a:1.0,b:1.0):0.5,(c:0.80001,d:0.79999):0.7);")
    assert main(["tools", "tree", f, "--round"]) == 0
    out = capsys.readouterr().out
    tree, _ = T.read_newick(out, assume_extant=True)
    d = [T._depths(tree)[i] for i, n in tree.nodes.items() if n.children is None]
    assert max(d) - min(d) == pytest.approx(0.0, abs=1e-9)


def test_cli_tree_red_values_emits_a_table(tmp_path, capsys):
    r = species.simulate_species_tree(birth=1.0, death=0.0, n_extant=6, seed=1)
    f = _write(tmp_path, "t.nwk", r.extant_tree.to_newick())
    assert main(["tools", "tree", f, "--red", "--values"]) == 0
    out = capsys.readouterr().out
    assert out.startswith("node\tRED")
    assert "\t0" in out                                     # the root is RED 0


def test_cli_tree_requires_exactly_one_action(tmp_path):
    f = _write(tmp_path, "t.nwk", "((a:1,b:1):1,(c:1,d:1):1);")
    with pytest.raises(SystemExit):                          # argparse: no action given
        main(["tools", "tree", f])


def test_cli_treedist_self_is_zero_and_mismatch_errors(tmp_path, capsys):
    r = species.simulate_species_tree(birth=1.0, death=0.0, n_extant=8, seed=3)
    f = _write(tmp_path, "t.nwk", r.extant_tree.to_newick())
    assert main(["tools", "treedist", f, f, "--metric", "all"]) == 0
    out = capsys.readouterr().out
    assert "rf\t0" in out and "branch-score\t0" in out
    other = _write(tmp_path, "other.nwk", "((a:1,b:1):1,(c:1,e:1):1);")   # taxon e, not d
    with pytest.raises(SystemExit):
        main(["tools", "treedist", f, other, "--metric", "rf"])


def test_cli_treedist_matches_external_trees_by_label_not_parse_order(tmp_path, capsys):
    # two external trees, same four taxa, DIFFERENT topology — must be rf>0 (matched by label,
    # not by the positionally-minted parse ids, which would falsely report rf 0).
    t1 = _write(tmp_path, "t1.nwk", "((A:1,B:1):1,(C:1,D:1):1);")
    t2 = _write(tmp_path, "t2.nwk", "((A:1,C:1):1,(B:1,D:1):1);")
    assert main(["tools", "treedist", t1, t2, "--metric", "rf"]) == 0
    assert capsys.readouterr().out.strip() == "rf\t4"
    assert main(["tools", "treedist", t1, t1, "--metric", "rf"]) == 0
    assert capsys.readouterr().out.strip() == "rf\t0"


def test_cli_treedist_compares_a_gene_tree_to_a_species_tree_and_says_so(tmp_path, capsys):
    # a gene tree's tips are genes and a species tree's are species, so left alone the two share no
    # labels at all. They ARE comparable on the species each gene sits in — which is the question
    # "does this family recover the species tree?" — but the reader has to be told that is what
    # happened, or a number appears from a comparison that was never like-for-like.
    sp = _write(tmp_path, "sp.nwk", "((n3:1,n4:1)n1:1,(n5:1,n6:1)n2:1)n0;")
    same = _write(tmp_path, "same.nwk",
                  "((n3_g1:1,n4_g2:1)x:1,(n5_g3:1,n6_g4:1)y:1)z;")      # same shape, gene tips
    assert main(["tools", "treedist", same, sp, "--metric", "rf"]) == 0
    cap = capsys.readouterr()
    assert cap.out.strip() == "rf\t0"
    assert "tips are gene copies" in cap.err and "species each gene sits in" in cap.err

    # and it is a real comparison, not one that always answers 0: a transfer-shaped gene tree
    moved = _write(tmp_path, "moved.nwk",
                   "((n3_g1:1,n5_g3:1)x:1,(n4_g2:1,n6_g4:1)y:1)z;")
    assert main(["tools", "treedist", moved, sp, "--metric", "rf"]) == 0
    assert capsys.readouterr().out.strip() == "rf\t4"


def test_cli_treedist_refuses_a_multi_copy_gene_tree_against_a_species_tree(tmp_path, capsys):
    # two copies in one genome means the gene -> species mapping is not one-to-one, so there is no
    # well-defined distance. A plausible number would be worse than a refusal.
    sp = _write(tmp_path, "sp.nwk", "((n3:1,n4:1)n1:1,(n5:1,n6:1)n2:1)n0;")
    dup = _write(tmp_path, "dup.nwk",
                 "(((n3_g1:1,n3_g9:1)d:1,n4_g2:1)x:1,(n5_g3:1,n6_g4:1)y:1)z;")
    with pytest.raises(SystemExit):
        main(["tools", "treedist", dup, sp, "--metric", "rf"])
    assert "n3 carry several copies" in capsys.readouterr().err


def test_cli_treedist_leaves_two_trees_of_the_same_kind_alone(tmp_path, capsys):
    # both gene trees: nothing to detect, nothing to say
    a = _write(tmp_path, "a.nwk", "((n3_g1:1,n4_g2:1)x:1,(n5_g3:1,n6_g4:1)y:1)z;")
    b = _write(tmp_path, "b.nwk", "((n3_g1:1,n5_g3:1)x:1,(n4_g2:1,n6_g4:1)y:1)z;")
    assert main(["tools", "treedist", a, b, "--metric", "rf"]) == 0
    cap = capsys.readouterr()
    assert cap.out.strip() == "rf\t4"
    assert "gene copies" not in cap.err


def test_a_tree_with_no_branch_lengths_is_refused():
    """A topology-only tree spans no time, so every rate — all of them per unit time — fires zero
    times. The run used to succeed having simulated nothing: a genome per node, a gene tree per
    family, and an event log holding only the originations and the speciations the topology forced.
    That is invisible from the outside, so it is refused at the door."""
    import pytest

    from zombi2.tree import read_newick

    with pytest.raises(ValueError, match="no branch lengths"):
        read_newick("(((A,B),C),D);")

    # ...but a geometric caller (treedist, the tree transforms) may still load one: a bare topology
    # is a perfectly good input to a tree comparison
    tree, _ = read_newick("(((A,B),C),D);", assume_extant=True)
    assert len(tree.nodes) == 7


def test_a_written_tree_is_still_ultrametric_when_it_is_read_back():
    """An extant tree from a dated run is ultrametric to ~1e-16 in memory, and the writer used to
    undo it: a tip's depth is a **sum** of branch lengths, so at the old 7 significant digits the
    rounding accumulated down the path and two tips came out ~1e-6 apart — far above the ~1e-8
    ``ape::is.ultrametric()`` allows. The very first thing anyone does with `species_extant.nwk` in R
    therefore failed, on a tree that was never not ultrametric.

    Checked after a round trip, because reading back is what was broken, and on the plain tree as
    well as the snapped one: ``--round`` was never the only thing affected."""
    from zombi2 import species
    from zombi2.tree import make_ultrametric, read_newick

    result = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=40, seed=1)

    def relative_spread(newick):
        tree, _ = read_newick(newick)
        depth = {}
        for i in sorted(tree.nodes):
            node = tree.nodes[i]
            depth[i] = (0.0 if node.parent is None else depth[node.parent]) + \
                       (node.end_time - node.birth_time)
        tips = [depth[i] for i, n in tree.nodes.items() if n.children is None]
        return (max(tips) - min(tips)) / max(tips)

    ape_tolerance = 1e-8
    assert relative_spread(result.extant_tree.to_newick()) < ape_tolerance
    assert relative_spread(make_ultrametric(result.extant_tree, tol=1e-3).to_newick()) < ape_tolerance
    # the old default is still reachable, and still not good enough — which is why it is not the default
    assert relative_spread(result.extant_tree.to_newick(precision=7)) > ape_tolerance


def test_pruning_keeps_the_complete_tree_s_clade_order():
    """The extant tree must draw its clades on the same side the complete tree does.

    Children were rebuilt in ascending node-id order. Ids are assigned in birth order, so that
    coincides with the original order until a child is pruned away — the surviving descendant that
    replaces it can have a far larger id than its sibling, and the pair comes out swapped. A figure
    showing the complete tree beside the extant one then put the same clade top-left in one and
    bottom-left in the other, and any reader joining the two by position disagreed with itself.

    Checked on tip order, which is what a reader actually sees, over enough seeds to catch it: the
    extant tips of the complete tree, in Newick order, are the extant tree's tips in that same order.
    """
    import re

    from zombi2 import species

    for seed in range(1, 25):
        result = species.simulate_species_tree(birth=1.0, death=0.45, n_extant=12, seed=seed)
        tips = lambda nwk: re.findall(r"[(,]([ne]\d+):", nwk)          # noqa: E731
        survivors = [t for t in tips(result.complete_tree.to_newick()) if t.startswith("n")]
        assert survivors == tips(result.extant_tree.to_newick()), f"seed {seed}"
