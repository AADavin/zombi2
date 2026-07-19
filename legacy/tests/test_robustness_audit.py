"""Regression tests for the robustness/consistency fixes (2026-07 audit).

Each of these was a crash, silent data corruption, or an internal inconsistency on a realistic but
untested edge input. Grouped by subsystem.
"""


import numpy as np
import pytest

from zombi2.tree import read_newick


# --------------------------------------------------------------- Newick parsing

def test_read_newick_strips_whitespace_in_names():
    # a space after a comma / a line-wrapped file must not leak into the following node name
    assert [n.name for n in read_newick("(A:1, B:1)root:0;").leaves()] == ["A", "B"]
    assert [n.name for n in read_newick("(A:1,\n  B:1);").leaves()] == ["A", "B"]


def test_read_newick_normal_tree_unchanged():
    t = read_newick("((a:1,b:1)i:1,(c:1,d:1)j:1)R:0;")
    assert {n.name for n in t.leaves()} == {"a", "b", "c", "d"}


def test_read_newick_quoted_label_with_space_stays_one_leaf():
    # a quoted label keeps its embedded space (and the quotes are stripped); the whitespace-skip
    # must not split it into two phantom leaves
    leaves = read_newick("('Homo sapiens':1,B:2);").leaves()
    assert [(n.name, n.branch_length()) for n in leaves] == [("Homo sapiens", 1.0), ("B", 2.0)]


# --------------------------------------------------------------- sequence models / FASTA

def test_read_fasta_blank_header(tmp_path):
    p = tmp_path / "x.fasta"
    p.write_text(">\nACGT\n>s2\nTTTT\n")            # blank first header
    from zombi2.sequences.models import read_fasta
    recs = read_fasta(str(p))
    assert recs == {"": "ACGT", "s2": "TTTT"}


def test_zero_stationary_frequency_rejected_clearly():
    from zombi2.sequences.models import make_model
    with pytest.raises(ValueError, match="strictly positive"):
        make_model("hky85", freqs=(0.5, 0.5, 0.0, 0.0))
    make_model("hky85", freqs=(0.25, 0.25, 0.25, 0.25))   # normal still builds


# --------------------------------------------------------------- profile / rate I/O

def test_from_coo_tsv_keeps_family_named_family():
    from zombi2.genomes.profiles import ProfileMatrix
    pm = ProfileMatrix.from_coo_tsv("family\tspecies\tcopies\nfamily\tA\t3\ng2\tA\t2\n")
    assert sorted(pm.copies_per_family().tolist()) == [2, 3]


def test_duplicate_coo_cells_dense_matches_reductions():
    from zombi2.genomes.profiles import ProfileMatrix
    pm = ProfileMatrix.from_coo_tsv("family\tspecies\tcopies\ng1\tA\t2\ng1\tA\t5\n")
    assert pm.copies_per_family().tolist() == [7]          # 2 + 5 summed (COO semantics)
    assert pm.matrix.sum(axis=1).tolist() == [7]           # dense view agrees (was last-write-wins)


def test_read_family_rates_requires_family_column(tmp_path):
    from zombi2.genomes.read_rates import read_family_rates
    bad = tmp_path / "bad.tsv"
    bad.write_text("duplication\ttransfer\tloss\n3\t2\t1\n")
    with pytest.raises(ValueError, match="family"):
        read_family_rates(str(bad))
    good = tmp_path / "good.tsv"
    good.write_text("family\tduplication\ttransfer\tloss\nfamA\t3\t2\t1\n")
    assert read_family_rates(str(good)) == {"famA": (3.0, 2.0, 1.0)}


# --------------------------------------------------------------- reconciliation consistency

def _species_tree():
    from zombi2.tools.reconciliation import SpeciesTree
    return SpeciesTree.from_tree(read_newick("((A:1,B:1)i1:1,(C:1,D:1)i2:1)root;"))


def test_dated_unknown_species_raises_like_undated():
    from zombi2.tools.reconciliation import GeneTree, dated_loglik, undated_loglik, DatedDTL, UndatedDTL
    sp = _species_tree()
    gt = GeneTree.from_newick("(A|1,Z|2)r;")               # Z is not a species-tree leaf
    with pytest.raises(KeyError):
        undated_loglik(gt, sp, UndatedDTL(0.1, 0.1, 0.1))
    with pytest.raises(KeyError):                          # dated used to silently return -inf
        dated_loglik(gt, sp, DatedDTL(0.1, 0.1, 0.1), n_steps=20)


def test_dated_normal_still_scores():
    from zombi2.tools.reconciliation import GeneTree, dated_loglik, DatedDTL
    ll = dated_loglik(GeneTree.from_newick("(A|1,B|2)r;"), _species_tree(),
                      DatedDTL(0.1, 0.1, 0.1), n_steps=20)
    assert np.isfinite(ll)


def test_dated_single_leaf_tree_clear_error():
    from zombi2.tools.reconciliation import SpeciesTree, GeneTree, dated_loglik, DatedDTL
    sp = SpeciesTree.from_tree(read_newick("A;"))
    with pytest.raises(ValueError, match="at least two leaves"):
        dated_loglik(GeneTree.from_newick("A|1;"), sp, DatedDTL(0.1, 0.05, 0.1))


def test_write_scores_tsv_derives_columns_from_rows(tmp_path):
    import types
    from zombi2.tools.reconciliation.scoring import write_scores_tsv
    rows = [types.SimpleNamespace(family="f1", extant_tips=3, logliks={"reldated": -5.0})]
    p = tmp_path / "scores.tsv"
    write_scores_tsv(rows, str(p))                         # default no longer hard-codes dated/undated
    assert p.read_text().splitlines()[0].endswith("reldated_loglik")


# --------------------------------------------------------------- clocks / traits

def test_rate_variation_root_rate_consistent():
    from zombi2.sequences.clocks import RateVariation
    rv = RateVariation(bins=[0.5, 2.0, 8.0], switch_rate=0.0, start=2)
    tree = read_newick("((A:1,B:1):1,C:2);")
    scaled = rv.scale(tree, rng=np.random.default_rng(0))
    assert rv.root_rate == scaled.branch_rate[tree.root] == 8.0


def test_pagel_delta_ultrametric_holds_tips_at_present():
    from zombi2.traits.models import pagel_delta
    t = read_newick("((A:1,B:1):1,(C:1,D:1):1);")         # ultrametric, tips at T=2
    for delta in (0.2, 1.0, 2.5):
        assert all(abs(n.time - 2.0) < 1e-9 for n in pagel_delta(t, delta).leaves())


def test_pagel_delta_keeps_all_branches_positive():
    # the depth rescaling is monotonic, so no branch may invert — including on a non-ultrametric
    # (fossil) tree and for delta < 1, which pushes internal splits toward the present
    from zombi2.traits.models import pagel_delta
    for nwk in ("((A:1,B:1):1,(C:1,D:1):1);", "((A:1,B:1,C:0.1)P:1)R:0;"):
        t = read_newick(nwk)
        for delta in (0.1, 0.5, 2.0):
            out = pagel_delta(t, delta)
            assert all(n.branch_length() >= -1e-12 for n in out.nodes_preorder()
                       if n.parent is not None), f"negative branch: {nwk} delta={delta}"


def test_musse_frozen_character_stationary_is_uniform():
    from zombi2.coevolve.sse import BiSSE
    pi = BiSSE(1, 3, 0.2, 0.2, 0.0, 0.0).stationary_distribution()   # frozen character (q=0)
    assert np.allclose(pi, [0.5, 0.5])
