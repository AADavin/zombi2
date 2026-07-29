"""Tests for :mod:`zombi2.tools.reconciliation` — the history projected onto what a dataset holds.

The complete reconciliation is the truth and is not what a method is scored against, because a method
sees only survivors. These check the projection: that it keeps exactly the observable events, that
what it does with the unobservable ones is defensible, and — the part most likely to go wrong — that
the result is still a *valid* reconciliation once nodes have been moved.
"""

import collections
import xml.etree.ElementTree as ET

import pytest

from zombi2.genomes import simulate_genomes_family
from zombi2.species import simulate_species_tree
from zombi2.tools.reconciliation import (ORIGINATION, TRANSFER, extant_reconciliation,
                                         origins_tsv, visible_branches)
from zombi2.tree import prune

#: regimes that differ in what is *invisible*: none, dead lineages, and unobserved survivors
_REGIMES = {
    "no extinction": (dict(birth=1.0, death=0.0, n_extant=25, seed=3),
                      dict(duplication=0.2, loss=0.3, origination=0.4, initial_families=30, seed=5)),
    "heavy extinction": (dict(birth=1.2, death=0.6, n_extant=40, seed=11),
                         dict(duplication=0.15, transfer=0.3, loss=0.2, origination=0.3,
                              initial_families=30, max_family_size=None, seed=11)),
    "unsampled survivors": (dict(birth=1.2, death=0.6, n_extant=30, sampling=0.6, seed=5),
                            dict(duplication=0.2, transfer=0.4, loss=0.3, initial_families=25,
                                 max_family_size=None, seed=5)),
}


def _run(name):
    sp_kw, g_kw = _REGIMES[name]
    sp = simulate_species_tree(**sp_kw)
    return sp, simulate_genomes_family(sp, **g_kw)


def _walk(root):
    stack, out = [root], []
    while stack:
        n = stack.pop()
        out.append(n)
        stack.extend(n.children)
    return out


def _extant_leaves(root):
    return sorted((n.species, n.copy) for n in _walk(root) if n.kind == "extant")


@pytest.mark.parametrize("regime", list(_REGIMES))
@pytest.mark.parametrize("scope", ["true", "recoverable"])
def test_the_projection_is_a_valid_reconciliation(regime, scope):
    # the part most likely to be wrong: nodes MOVE when their own branch is not one the extant tree
    # names, and a moved node can easily end up somewhere its parent cannot reach. Checked here
    # against the trees themselves rather than against the builder's own bookkeeping.
    sp, g = _run(regime)
    tree = sp.complete_tree
    et = prune(tree, keep="extant")
    vis = visible_branches(tree)
    ancestors, tips = {}, {i for i, n in et.nodes.items() if n.children is None}
    for i in sorted(et.nodes):
        p = et.nodes[i].parent
        ancestors[i] = (ancestors[p] | {p}) if p is not None else set()

    checked = 0
    for fam, gt in g.gene_trees.items():
        rec = extant_reconciliation(gt, tree, scope=scope, visible=vis)
        if rec is None:
            continue
        for n in _walk(rec.root):
            checked += 1
            assert n.species in et.nodes, "every branch named must be one the extant tree has"
            if n.is_leaf:
                assert n.kind in ("extant", "loss")
                if n.kind == "extant":
                    assert n.species in tips
            elif n.kind in ("speciation", "duplication"):
                # vertical descent only: a child sits on its parent's branch or below it. (A transfer
                # is the exception by definition — that is what makes it a transfer.)
                for c in n.children:
                    assert c.species == n.species or n.species in ancestors[c.species]
    assert checked > 100


@pytest.mark.parametrize("regime", list(_REGIMES))
def test_the_reconciliation_is_of_the_gene_tree_a_dataset_holds(regime):
    # whatever the projection does to internal nodes, the observable genes must come through
    # untouched — otherwise it is a reconciliation of some other tree
    sp, g = _run(regime)
    vis = visible_branches(sp.complete_tree)
    seen = 0
    for fam, gt in g.gene_trees.items():
        if gt.extant is None:
            continue
        seen += 1
        for scope in ("true", "recoverable"):
            rec = extant_reconciliation(gt, sp.complete_tree, scope=scope, visible=vis)
            assert _extant_leaves(rec.root) == _extant_leaves(gt.extant)
    assert seen > 5


def test_true_contains_recoverable():
    # the whole reason both are written: 'true' carries the family's ancestral presence and the
    # losses that narrowed it, 'recoverable' is that trimmed back to the surviving copies' ancestor.
    # Trimming is possible, inventing is not, so the richer one is the one to ship.
    sp, g = _run("heavy extinction")
    vis = visible_branches(sp.complete_tree)
    extra, pairs = 0, 0
    for fam, gt in g.gene_trees.items():
        t = extant_reconciliation(gt, sp.complete_tree, scope="true", visible=vis)
        r = extant_reconciliation(gt, sp.complete_tree, scope="recoverable", visible=vis)
        if t is None:
            continue
        pairs += 1
        assert t.losses >= r.losses                       # only ever adds
        assert len(_walk(t.root)) >= len(_walk(r.root))
        extra += t.losses - r.losses
    assert pairs > 5 and extra > 0, "some family must be ancestrally present above its survivors"


def test_a_family_that_arrived_from_a_dead_lineage_says_so():
    # the marker exists because the arrival is INVISIBLE in the reconciliation itself: with the donor
    # gone the copy simply appears, which a reader would call an origination and be wrong about
    sp, g = _run("heavy extinction")
    vis = visible_branches(sp.complete_tree)
    kinds = collections.Counter()
    for fam, gt in g.gene_trees.items():
        rec = extant_reconciliation(gt, sp.complete_tree, visible=vis)
        if rec is not None:
            kinds[rec.entered_by] += 1
    assert kinds[ORIGINATION] > 0 and kinds[TRANSFER] > 0
    assert set(kinds) == {ORIGINATION, TRANSFER}


def test_with_no_extinction_there_is_nothing_to_move():
    # the control. Every lineage survives, so every branch is one the extant tree names and each is
    # its own image: no node can be relocated and no transfer synthesised, which is what isolates
    # those behaviours to the regimes that actually have invisible branches.
    #
    # The two scopes still differ here, and that is the point worth pinning: 'true' vs 'recoverable'
    # is about ancestral presence above the surviving copies, NOT about extinction. A family present
    # at the root and lost outside one clade has history no method can reach even when every lineage
    # is sampled.
    sp, g = _run("no extinction")
    tree = sp.complete_tree
    kept, image = visible_branches(tree)
    assert kept == set(tree.nodes)
    assert image == {i: i for i in tree.nodes}
    for fam, gt in g.gene_trees.items():
        rec = extant_reconciliation(gt, tree, visible=(kept, image))
        if rec is not None:
            assert rec.entered_by == ORIGINATION      # no dead donor to arrive from


def test_scope_is_validated():
    sp, g = _run("no extinction")
    gt = next(iter(g.gene_trees.values()))
    with pytest.raises(ValueError, match="'true' or 'recoverable'"):
        extant_reconciliation(gt, sp.complete_tree, scope="best")


def test_origins_tsv_names_every_family_once():
    sp, g = _run("heavy extinction")
    vis = visible_branches(sp.complete_tree)
    recs = {f: r for f, gt in g.gene_trees.items()
            if (r := extant_reconciliation(gt, sp.complete_tree, visible=vis)) is not None}
    rows = origins_tsv(recs).splitlines()
    assert rows[0] == "family\tentered_by\tbranch\tlosses"
    assert len(rows) == len(recs) + 1
    assert all(r.split("\t")[1] in (ORIGINATION, TRANSFER) for r in rows[1:])
    assert all(r.split("\t")[2].startswith("n") for r in rows[1:])   # a branch of the EXTANT tree


def test_the_cli_writes_both_scopes_and_the_origins_table(tmp_path):
    from zombi2.cli.main import main

    run = tmp_path / "r"
    assert main(["species", str(run), "--birth", "1.2", "--death", "0.6", "--n-extant", "20",
                 "--seed", "11", "--quiet"]) == 0
    assert main(["genomes", str(run), "--duplication", "0.15", "--transfer", "0.3", "--loss", "0.2",
                 "--origination", "0.3", "--initial-families", "10", "--seed", "11", "--quiet"]) == 0
    assert main(["tools", "format", str(run), "--format", "recphylo", "--recphylo", "extant",
                 "--quiet"]) == 0
    out = run / "genomes" / "recphylo"
    assert (out / "family_origins.tsv").exists()
    assert not list(out.glob("recphylo_fam*[0-9].xml")), "extant scope writes no complete files"
    pairs = sorted(p.name for p in out.glob("*_true.xml"))
    assert pairs and len(pairs) == len(list(out.glob("*_recoverable.xml")))

    for p in out.glob("*_true.xml"):
        root = ET.fromstring(p.read_text(encoding="utf-8"))       # it is well-formed XML...
        names = {c.text for c in root.find("spTree").iter("name")}
        assert names and not any(n.startswith("e") for n in names)   # ...against the EXTANT tree
